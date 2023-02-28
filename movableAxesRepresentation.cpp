#include <vtkActor.h>
#include <vtkBox.h>
#include <vtkCallbackCommand.h>
#include <vtkCamera.h>
#include <vtkCellPicker.h>
#include <vtkConeSource.h>
#include <vtkFollower.h>
#include <vtkInteractorObserver.h>
#include <vtkLine.h>
#include <vtkLineSource.h>
#include <vtkMath.h>
#include <vtkObjectFactory.h>
#include <vtkPointHandleRepresentation3D.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSphereSource.h>
#include <vtkVector.h>
#include <vtkVectorOperators.h>
#include <vtkVectorText.h>
#include <vtkWindow.h>

#include "movableAxesRepresentation.hpp"

vtkStandardNewMacro(movableAxesRepresentation);

movableAxesRepresentation::movableAxesRepresentation() {}

movableAxesRepresentation::~movableAxesRepresentation() {}

void movableAxesRepresentation::SetRenderer(vtkRenderer *ren) {
  this->Superclass::SetRenderer(ren);
}

void movableAxesRepresentation::StartWidgetInteraction(double e[2]) {}

void movableAxesRepresentation::WidgetInteraction(double e[2]) {}

void movableAxesRepresentation::PlaceWidget(double bds[6]) {}

int movableAxesRepresentation::ComputeInteractionState(int x, int y,
                                                       int vtkNotUsed(modify)) {
  return 0;
}

double *movableAxesRepresentation::GetBounds() {
  double a;
  return &a;
}

void movableAxesRepresentation::BuildRepresentation() {}

//------------------------------------------------------------------------------
vtkMTimeType movableAxesRepresentation::GetMTime() {
  vtkMTimeType mTime = this->Superclass::GetMTime();
  return mTime;
}

//------------------------------------------------------------------------------
void movableAxesRepresentation::PrintSelf(ostream &os, vtkIndent indent) {
  this->Superclass::PrintSelf(os, indent);
}
