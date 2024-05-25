#include <vtkActor.h>
#include <vtkCallbackCommand.h>
#include <vtkCollection.h>
#include <vtkCommand.h>
#include <vtkObjectFactory.h>
#include <vtkPropCollection.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkWidgetCallbackMapper.h>
#include <vtkWidgetEvent.h>

#include "simpleTransformRepresentation.hpp"
#include "transformWidget.hpp"

vtkStandardNewMacro(transformWidget);

transformWidget::transformWidget()
{
    this->state_ = WIDGETSTATE::start;

    // define widget events
    {
        this->CallbackMapper->SetCallbackMethod(
            vtkCommand::LeftButtonPressEvent, vtkWidgetEvent::Select, this,
            transformWidget::SelectAction);

        this->CallbackMapper->SetCallbackMethod(
            vtkCommand::LeftButtonReleaseEvent, vtkWidgetEvent::EndSelect, this,
            transformWidget::EndSelectAction);

        this->CallbackMapper->SetCallbackMethod(vtkCommand::MouseMoveEvent,
                                                vtkWidgetEvent::Move, this,
                                                transformWidget::MoveAction);

        this->keyEventCallbackCommand_ = vtkCallbackCommand::New();
        this->keyEventCallbackCommand_->SetClientData(this);
        this->keyEventCallbackCommand_->SetCallback(
            transformWidget::ProcessKeyEvents);
    }
}

transformWidget::~transformWidget()
{
    this->keyEventCallbackCommand_->Delete();
}

void transformWidget::SetEnabled(int enabling)
{
    const int enabled = this->Enabled;

    if (enabling && !enabled)
    {
        // We do this step first because it sets the CurrentRenderer
        vtkAbstractWidget::SetEnabled(enabling);

        this->CreateDefaultRepresentation();

        if (this->Parent)
        {
            this->Parent->AddObserver(vtkCommand::KeyPressEvent,
                                      this->keyEventCallbackCommand_,
                                      this->Priority);
            this->Parent->AddObserver(vtkCommand::KeyReleaseEvent,
                                      this->keyEventCallbackCommand_,
                                      this->Priority);
        }
        else
        {
            this->Interactor->AddObserver(vtkCommand::KeyPressEvent,
                                          this->keyEventCallbackCommand_,
                                          this->Priority);
            this->Interactor->AddObserver(vtkCommand::KeyReleaseEvent,
                                          this->keyEventCallbackCommand_,
                                          this->Priority);
        }

        // Add the actor
        {
            vtkNew<vtkPropCollection> propsCollection;
            this->WidgetRep->GetActors(propsCollection);

            vtkCollectionSimpleIterator sIt;
            propsCollection->InitTraversal(sIt);

            const int nProps = propsCollection->GetNumberOfItems();
            for (int i = 0; i < nProps; i++)
            {
                vtkProp *actor =
                    vtkProp::SafeDownCast(propsCollection->GetNextProp(sIt));
                if (actor != nullptr)
                {
                    this->CurrentRenderer->AddActor(actor);
                }
            }
        }
    }
    else if (!enabling && enabled)
    {
        if (this->Parent)
        {
            this->Parent->RemoveObserver(this->keyEventCallbackCommand_);
        }
        else
        {
            this->Interactor->RemoveObserver(this->keyEventCallbackCommand_);
        }

        // Remove the actor
        {
            vtkNew<vtkPropCollection> propsCollection;
            this->WidgetRep->GetActors(propsCollection);

            vtkCollectionSimpleIterator sIt;
            propsCollection->InitTraversal(sIt);

            const int nProps = propsCollection->GetNumberOfItems();
            for (int i = 0; i < nProps; i++)
            {
                vtkProp *actor =
                    vtkProp::SafeDownCast(propsCollection->GetNextProp(sIt));
                if (actor != nullptr)
                {
                    this->CurrentRenderer->RemoveActor(actor);
                }
            }
        }

        // We do this step last because the CurrentRenderer is still useful
        vtkAbstractWidget::SetEnabled(enabling);
    }

    this->Interactor->Render();
}

void transformWidget::CreateDefaultRepresentation()
{
    if (!this->WidgetRep)
    {
        this->WidgetRep = simpleTransformRepresentation::New();
    }
}

void transformWidget::SetRepresentation(transformRepresentation *rep)
{
    this->Superclass::SetWidgetRepresentation(
        vtkWidgetRepresentation::SafeDownCast(rep));
}

void transformWidget::PrintSelf(ostream &os, vtkIndent indent)
{
    vtkAbstractWidget::PrintSelf(os, indent);
}

void transformWidget::SelectAction(vtkAbstractWidget *w)
{
    auto self = reinterpret_cast<transformWidget *>(w);
    if (self->WidgetRep->GetInteractionState() == INTERACTIONSTATE::outside)
    {
        return;
    }

    const int x = self->Interactor->GetEventPosition()[0];
    const int y = self->Interactor->GetEventPosition()[1];

    self->state_ = transformWidget::WIDGETSTATE::active;

    self->GrabFocus(self->EventCallbackCommand);

    double eventPos[2];
    eventPos[0] = static_cast<double>(x);
    eventPos[1] = static_cast<double>(y);
    self->WidgetRep->StartWidgetInteraction(eventPos);

    // start the interaction
    self->EventCallbackCommand->SetAbortFlag(1);
    self->StartInteraction();
    self->InvokeEvent(vtkCommand::StartInteractionEvent, nullptr);
    self->Render();
}

void transformWidget::EndSelectAction(vtkAbstractWidget *w)
{
    auto self = reinterpret_cast<transformWidget *>(w);
    if (self->state_ == transformWidget::WIDGETSTATE::start)
    {
        return;
    }

    // Return state to not active
    self->state_ = transformWidget::WIDGETSTATE::start;
    self->ReleaseFocus();
    self->InvokeEvent(vtkCommand::LeftButtonReleaseEvent,
                      nullptr);  // handles observe this
    self->EventCallbackCommand->SetAbortFlag(1);
    self->InvokeEvent(vtkCommand::EndInteractionEvent, nullptr);
    self->EndInteraction();
    self->Render();
}

void transformWidget::MoveAction(vtkAbstractWidget *w)
{
    auto self = reinterpret_cast<transformWidget *>(w);

    const int x = self->Interactor->GetEventPosition()[0];
    const int y = self->Interactor->GetEventPosition()[1];

    if (self->state_ ==
        transformWidget::WIDGETSTATE::start)  // if is not active
    {
        self->Interactor->Disable();  // avoid extra renders

        const int prevState = self->WidgetRep->GetInteractionState();
        const int currState = self->WidgetRep->ComputeInteractionState(x, y);
        int cursorChanged = 0;

        {
            if (currState == INTERACTIONSTATE::outside)
            {
                cursorChanged = self->RequestCursorShape(VTK_CURSOR_DEFAULT);
            }
            else
            {
                cursorChanged = self->RequestCursorShape(VTK_CURSOR_HAND);
            }
        }

        self->WidgetRep->Highlight(1);

        self->Interactor->Enable();  // avoid extra renders

        if (cursorChanged || (prevState != currState))
        {
            self->Render();
        }
    }
    else  // if is active
    {
        double e[2];
        e[0] = static_cast<double>(x);
        e[1] = static_cast<double>(y);

        self->InvokeEvent(vtkCommand::MouseMoveEvent,
                          nullptr);  // handles observe this

        self->WidgetRep->WidgetInteraction(e);

        self->EventCallbackCommand->SetAbortFlag(1);
        self->InvokeEvent(vtkCommand::InteractionEvent, nullptr);
        self->Render();
    }
}

void transformWidget::ProcessKeyEvents(vtkObject *, unsigned long, void *,
                                       void *)
{
}
