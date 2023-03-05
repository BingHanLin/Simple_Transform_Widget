#include <vtkActor.h>
#include <vtkAssemblyPath.h>
#include <vtkCallbackCommand.h>
#include <vtkCellPicker.h>
#include <vtkInteractorObserver.h>
#include <vtkObjectFactory.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRegularPolygonSource.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkTransform.h>
#include <vtkWindow.h>

#include "movableAxesRepresentation.hpp"

vtkStandardNewMacro(movableAxesRepresentation);

movableAxesRepresentation::movableAxesRepresentation()
{
    for (auto i = 0; i < 3; i++)
    {
        axisRingActors_[i] = vtkSmartPointer<vtkActor>::New();

        vtkNew<vtkRegularPolygonSource> polygonSource;
        polygonSource->SetNumberOfSides(100);
        polygonSource->SetRadius(0.5);
        polygonSource->SetCenter(0, 0, 0);

        vtkNew<vtkPolyDataMapper> mapper;
        mapper->SetInputConnection(polygonSource->GetOutputPort());

        axisRingActors_[i]->SetMapper(mapper);
        axisRingActors_[i]->GetProperty()->SetOpacity(1.0);
        axisRingActors_[i]->GetProperty()->SetLineWidth(5);
        axisRingActors_[i]->GetProperty()->SetRepresentationToWireframe();
        axisRingActors_[i]->GetProperty()->SetColor(axisNormalColor_[i][0],
                                                    axisNormalColor_[i][1],
                                                    axisNormalColor_[i][2]);
    }

    {
        vtkSmartPointer<vtkTransform> trans =
            vtkSmartPointer<vtkTransform>::New();
        trans->RotateY(90);
        axisRingActors_[0]->SetUserTransform(trans);
    }

    {
        vtkSmartPointer<vtkTransform> trans =
            vtkSmartPointer<vtkTransform>::New();
        trans->RotateX(90);
        axisRingActors_[1]->SetUserTransform(trans);
    }

    {
        vtkSmartPointer<vtkTransform> trans =
            vtkSmartPointer<vtkTransform>::New();
        axisRingActors_[2]->SetUserTransform(trans);
    }

    {
        picker_ = vtkSmartPointer<vtkCellPicker>::New();
        picker_->SetTolerance(0.001);
        for (auto i = 0; i < 3; i++)
        {
            picker_->AddPickList(axisRingActors_[i]);
        }
        picker_->PickFromListOn();
    }

    // Define the point coordinates
    double bounds[6];
    bounds[0] = -0.5;
    bounds[1] = 0.5;
    bounds[2] = -0.5;
    bounds[3] = 0.5;
    bounds[4] = -0.5;
    bounds[5] = 0.5;
    this->PlaceWidget(bounds);
}

movableAxesRepresentation::~movableAxesRepresentation()
{
}

void movableAxesRepresentation::SetRenderer(vtkRenderer *ren)
{
    this->Superclass::SetRenderer(ren);
}

void movableAxesRepresentation::StartWidgetInteraction(double e[2])
{
}

void movableAxesRepresentation::WidgetInteraction(double e[2])
{
}

void movableAxesRepresentation::PlaceWidget(double bds[6])
{
}

int movableAxesRepresentation::ComputeInteractionState(int x, int y,
                                                       int vtkNotUsed(modify))
{
    if (!this->Renderer || !this->Renderer->IsInViewport(x, y))
    {
        this->InteractionState = INTERACTIONSTATE::outside;
        return this->InteractionState;
    }

    auto path = this->GetAssemblyPath(x, y, 0., picker_);
    if (path != nullptr)
    {
        currActor_ =
            reinterpret_cast<vtkActor *>(path->GetFirstNode()->GetViewProp());

        if (currActor_ == axisRingActors_[0])
        {
            this->InteractionState = INTERACTIONSTATE::onXRing;
        }
        else if (currActor_ == axisRingActors_[1])
        {
            this->InteractionState = INTERACTIONSTATE::onYRing;
        }
        else if (currActor_ == axisRingActors_[2])
        {
            this->InteractionState = INTERACTIONSTATE::onZRing;
        }
    }
    else
    {
        this->InteractionState = INTERACTIONSTATE::outside;
    }

    return this->InteractionState;
}

double *movableAxesRepresentation::GetBounds()
{
    double bounds[6];
    bounds[0] = -0.5;
    bounds[1] = 0.5;
    bounds[2] = -0.5;
    bounds[3] = 0.5;
    bounds[4] = -0.5;
    bounds[5] = 0.5;
    return bounds;
}

void movableAxesRepresentation::BuildRepresentation()
{
}

vtkMTimeType movableAxesRepresentation::GetMTime()
{
    vtkMTimeType mTime = this->Superclass::GetMTime();
    return mTime;
}

void movableAxesRepresentation::GetActors(vtkPropCollection *pc)
{
    for (auto i = 0; i < 3; i++)
    {
        axisRingActors_[i]->GetActors(pc);
    }
}

void movableAxesRepresentation::PrintSelf(ostream &os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);
}
