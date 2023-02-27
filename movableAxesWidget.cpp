#include <vtkObjectFactory.h>

#include "movableAxesWidget.hpp"
#include "myLineRepresentation.h"

vtkStandardNewMacro(movableAxesWidget);

movableAxesWidget::movableAxesWidget() {}

movableAxesWidget::~movableAxesWidget() {}

void movableAxesWidget::SetEnabled(int enabled) {}

void movableAxesWidget::CreateDefaultRepresentation() {
  if (!this->WidgetRep) {
    this->WidgetRep = myLineRepresentation::New();
  }
}

void movableAxesWidget::SetProcessEvents(vtkTypeBool enabled) {}

void movableAxesWidget::PrintSelf(ostream &os, vtkIndent indent) {
  this->Superclass::PrintSelf(os, indent);
}
