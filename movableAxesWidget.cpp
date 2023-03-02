#include <vtkActor.h>
#include <vtkCallbackCommand.h>
#include <vtkCollection.h>
#include <vtkCommand.h>
#include <vtkObjectFactory.h>
#include <vtkPropCollection.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>

#include "movableAxesRepresentation.hpp"
#include "movableAxesWidget.hpp"

vtkStandardNewMacro(movableAxesWidget);

movableAxesWidget::movableAxesWidget() {
  this->state_ = WIDGETSTATE::Start;

  this->keyEventCallbackCommand_ = vtkCallbackCommand::New();
  this->keyEventCallbackCommand_->SetClientData(this);
  // this->keyEventCallbackCommand_->SetCallback(movableAxesWidget::ProcessKeyEvents);
}

movableAxesWidget::~movableAxesWidget() {
  this->keyEventCallbackCommand_->Delete();
}

void movableAxesWidget::SetEnabled(int enabling) {

  const int enabled = this->Enabled;

  // We do this step first because it sets the CurrentRenderer
  vtkAbstractWidget::SetEnabled(enabling);

  if (enabling && !enabled) {
    this->CreateDefaultRepresentation();

    if (this->Parent) {
      this->Parent->AddObserver(vtkCommand::KeyPressEvent,
                                this->keyEventCallbackCommand_, this->Priority);
      this->Parent->AddObserver(vtkCommand::KeyReleaseEvent,
                                this->keyEventCallbackCommand_, this->Priority);
    } else {
      this->Interactor->AddObserver(vtkCommand::KeyPressEvent,
                                    this->keyEventCallbackCommand_,
                                    this->Priority);
      this->Interactor->AddObserver(vtkCommand::KeyReleaseEvent,
                                    this->keyEventCallbackCommand_,
                                    this->Priority);
    }

    // Add the ring actor
    vtkNew<vtkPropCollection> propsCollection;
    this->WidgetRep->GetActors(propsCollection);

    vtkCollectionSimpleIterator sIt;
    propsCollection->InitTraversal(sIt);

    const int nProps = propsCollection->GetNumberOfItems();
    for (int i = 0; i < nProps; i++) {
      vtkActor *actor =
          vtkActor::SafeDownCast(propsCollection->GetNextProp(sIt));
      if (actor != nullptr) {
        this->CurrentRenderer->AddActor(actor);
      }
    }

  } else if (!enabling && enabled) {
    if (this->Parent) {
      this->Parent->RemoveObserver(this->keyEventCallbackCommand_);
    } else {
      this->Interactor->RemoveObserver(this->keyEventCallbackCommand_);
    }
  }
}

void movableAxesWidget::CreateDefaultRepresentation() {
  if (!this->WidgetRep) {
    this->WidgetRep = movableAxesRepresentation::New();
  }
}

void movableAxesWidget::SetProcessEvents(vtkTypeBool enabled) {}

void movableAxesWidget::PrintSelf(ostream &os, vtkIndent indent) {
  vtkAbstractWidget::PrintSelf(os, indent);
}
