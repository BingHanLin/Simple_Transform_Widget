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

#include "movableAxesRepresentation.hpp"
#include "movableAxesWidget.hpp"

vtkStandardNewMacro(movableAxesWidget);

movableAxesWidget::movableAxesWidget()
{
    this->state_ = WIDGETSTATE::start;

    // define widget events
    {
        this->CallbackMapper->SetCallbackMethod(
            vtkCommand::LeftButtonPressEvent, vtkWidgetEvent::Select, this,
            movableAxesWidget::SelectAction);

        this->CallbackMapper->SetCallbackMethod(vtkCommand::MouseMoveEvent,
                                                vtkWidgetEvent::Move, this,
                                                movableAxesWidget::MoveAction);

        this->keyEventCallbackCommand_ = vtkCallbackCommand::New();
        this->keyEventCallbackCommand_->SetClientData(this);
        // this->keyEventCallbackCommand_->SetCallback(movableAxesWidget::ProcessKeyEvents);
    }
}

movableAxesWidget::~movableAxesWidget()
{
    this->keyEventCallbackCommand_->Delete();
}

void movableAxesWidget::SetEnabled(int enabling)
{
    const int enabled = this->Enabled;

    // We do this step first because it sets the CurrentRenderer
    vtkAbstractWidget::SetEnabled(enabling);

    if (enabling && !enabled)
    {
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

        // Add the ring actor
        {
            vtkNew<vtkPropCollection> propsCollection;
            this->WidgetRep->GetActors(propsCollection);

            vtkCollectionSimpleIterator sIt;
            propsCollection->InitTraversal(sIt);

            const int nProps = propsCollection->GetNumberOfItems();
            for (int i = 0; i < nProps; i++)
            {
                vtkActor *actor =
                    vtkActor::SafeDownCast(propsCollection->GetNextProp(sIt));
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
    }
}

void movableAxesWidget::CreateDefaultRepresentation()
{
    if (!this->WidgetRep)
    {
        this->WidgetRep = movableAxesRepresentation::New();
    }
}

void movableAxesWidget::SetProcessEvents(vtkTypeBool enabled)
{
}

void movableAxesWidget::PrintSelf(ostream &os, vtkIndent indent)
{
    vtkAbstractWidget::PrintSelf(os, indent);
}

void movableAxesWidget::SelectAction(vtkAbstractWidget *w)
{
    std::cout << "SelectAction......." << std::endl;
    auto self = reinterpret_cast<movableAxesWidget *>(w);
    if (self->WidgetRep->GetInteractionState() ==
        movableAxesRepresentation::INTERACTIONSTATE::outside)
    {
        return;
    }
}
void movableAxesWidget::MoveAction(vtkAbstractWidget *w)
{
    std::cout << "MoveAction......." << std::endl;
    auto self = reinterpret_cast<movableAxesWidget *>(w);

    // compute some info we need for all actions
    const int x = self->Interactor->GetEventPosition()[0];
    const int y = self->Interactor->GetEventPosition()[1];

    if (self->state_ == movableAxesWidget::WIDGETSTATE::start)
    {
        self->Interactor->Disable();  // avoid extra renders

        const int prevState = self->WidgetRep->GetInteractionState();
        const int currState = self->WidgetRep->ComputeInteractionState(x, y);
        int cursorChanged = 0;
        std::cout << "currState: " << currState << std::endl;

        if (currState == movableAxesRepresentation::INTERACTIONSTATE::outside)
        {
            cursorChanged = self->RequestCursorShape(VTK_CURSOR_DEFAULT);
        }
        else
        {
            cursorChanged = self->RequestCursorShape(VTK_CURSOR_HAND);

            // if (currState ==
            //         movableAxesRepresentation::INTERACTIONSTATE::onXRing ||
            //     currState ==
            //         movableAxesRepresentation::INTERACTIONSTATE::onYRing ||
            //     currState ==
            //         movableAxesRepresentation::INTERACTIONSTATE::onZRing)
            // {
            //     self->Point1Widget->SetEnabled(1);
            // }
        }
        self->Interactor->Enable();

        if (cursorChanged || prevState != currState)
        {
            self->Render();
        }
    }
    else  // if ( self->WidgetState == vtkLineWidget2::Active )
    {
    }
}