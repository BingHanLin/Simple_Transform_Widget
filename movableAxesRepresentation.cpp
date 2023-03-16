#include <vtkActor.h>
#include <vtkAssemblyPath.h>
#include <vtkCallbackCommand.h>
#include <vtkCellPicker.h>
#include <vtkConeSource.h>
#include <vtkInteractorObserver.h>
#include <vtkObjectFactory.h>
#include <vtkParametricFunctionSource.h>
#include <vtkParametricTorus.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkTransform.h>
#include <vtkWindow.h>

#include "movableAxesRepresentation.hpp"

vtkStandardNewMacro(movableAxesRepresentation);

movableAxesRepresentation::movableAxesRepresentation()
{
    const std::array<std::array<double, 3>, 3> coneDirection = {
        {{0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0}}};

    const double ringRadius = 0.5;
    const double ringSectionRadius = 0.025;

    for (auto i = 0; i < 3; i++)
    {
        {
            axisRingActors_[i] = vtkSmartPointer<vtkActor>::New();

            vtkNew<vtkParametricTorus> paramTorus;
            paramTorus->SetRingRadius(ringRadius);
            paramTorus->SetCrossSectionRadius(ringSectionRadius);

            vtkNew<vtkParametricFunctionSource> paramSource;
            paramSource->SetParametricFunction(paramTorus);
            paramSource->SetUResolution(60);
            paramSource->SetVResolution(15);
            paramSource->Update();

            vtkNew<vtkPolyDataMapper> mapper;
            mapper->SetInputConnection(paramSource->GetOutputPort());

            axisRingActors_[i]->SetMapper(mapper);
            axisRingActors_[i]->GetProperty()->SetDiffuse(0.8);
            axisRingActors_[i]->GetProperty()->SetSpecular(0.5);
            axisRingActors_[i]->GetProperty()->SetSpecularPower(30.0);
            axisRingActors_[i]->GetProperty()->SetOpacity(1.0);
            axisRingActors_[i]->GetProperty()->SetLineWidth(5);
            axisRingActors_[i]->GetProperty()->SetRepresentationToWireframe();
            axisRingActors_[i]->GetProperty()->SetColor(axisNormalColor_[i][0],
                                                        axisNormalColor_[i][1],
                                                        axisNormalColor_[i][2]);
        }

        {
            axisRingConeActors_[i] = vtkSmartPointer<vtkActor>::New();

            vtkNew<vtkConeSource> coneSource;
            coneSource->SetResolution(20);
            coneSource->SetHeight(0.25);
            coneSource->SetAngle(20);
            coneSource->SetDirection(coneDirection[i].data());
            coneSource->Update();

            vtkNew<vtkPolyDataMapper> mapper;
            mapper->SetInputConnection(coneSource->GetOutputPort());

            axisRingConeActors_[i]->SetMapper(mapper);
            axisRingConeActors_[i]->GetProperty()->SetDiffuse(0.8);
            axisRingConeActors_[i]->GetProperty()->SetSpecular(0.5);
            axisRingConeActors_[i]->GetProperty()->SetSpecularPower(30.0);
            axisRingConeActors_[i]->GetProperty()->SetColor(
                axisNormalColor_[i][0], axisNormalColor_[i][1],
                axisNormalColor_[i][2]);
        }
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
        vtkSmartPointer<vtkTransform> trans =
            vtkSmartPointer<vtkTransform>::New();
        trans->Translate(0.0, 0.0, -ringRadius);
        axisRingConeActors_[0]->SetUserTransform(trans);
    }

    {
        vtkSmartPointer<vtkTransform> trans =
            vtkSmartPointer<vtkTransform>::New();
        trans->Translate(-ringRadius, 0.0, 0.0);
        axisRingConeActors_[1]->SetUserTransform(trans);
    }

    {
        vtkSmartPointer<vtkTransform> trans =
            vtkSmartPointer<vtkTransform>::New();
        trans->Translate(0.0, -ringRadius, 0.0);
        axisRingConeActors_[2]->SetUserTransform(trans);
    }

    {
        picker_ = vtkSmartPointer<vtkCellPicker>::New();
        picker_->SetTolerance(0.001);
        for (auto i = 0; i < 3; i++)
        {
            picker_->AddPickList(axisRingActors_[i]);
            picker_->AddPickList(axisRingConeActors_[i]);
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
    startEventPosition_ = {e[0], e[1], 0.0};
    lastEventPosition_ = {e[0], e[1], 0.0};

    this->ComputeInteractionState(static_cast<int>(e[0]),
                                  static_cast<int>(e[1]), 0);
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
    for (auto i = 0; i < 3; i++)
    {
        axisRingActors_[i]->GetProperty()->SetDiffuse(0.8);
        axisRingActors_[i]->GetProperty()->SetSpecular(0.5);
        axisRingConeActors_[i]->GetProperty()->SetDiffuse(0.8);
        axisRingConeActors_[i]->GetProperty()->SetSpecular(0.5);
    }

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

        if (currActor_ == axisRingActors_[0] ||
            currActor_ == axisRingConeActors_[0])
        {
            this->InteractionState = INTERACTIONSTATE::onXRing;

            axisRingActors_[0]->GetProperty()->SetDiffuse(1.0);
            axisRingActors_[0]->GetProperty()->SetSpecular(0.0);
            axisRingConeActors_[0]->GetProperty()->SetDiffuse(1.0);
            axisRingConeActors_[0]->GetProperty()->SetSpecular(0.0);
        }
        else if (currActor_ == axisRingActors_[1] ||
                 currActor_ == axisRingConeActors_[1])
        {
            this->InteractionState = INTERACTIONSTATE::onYRing;

            axisRingActors_[1]->GetProperty()->SetDiffuse(1.0);
            axisRingActors_[1]->GetProperty()->SetSpecular(0.0);
            axisRingConeActors_[1]->GetProperty()->SetDiffuse(1.0);
            axisRingConeActors_[1]->GetProperty()->SetSpecular(0.0);
        }
        else if (currActor_ == axisRingActors_[2] ||
                 currActor_ == axisRingConeActors_[2])
        {
            this->InteractionState = INTERACTIONSTATE::onZRing;

            axisRingActors_[2]->GetProperty()->SetDiffuse(1.0);
            axisRingActors_[2]->GetProperty()->SetSpecular(0.0);
            axisRingConeActors_[2]->GetProperty()->SetDiffuse(1.0);
            axisRingConeActors_[2]->GetProperty()->SetSpecular(0.0);
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
        axisRingConeActors_[i]->GetActors(pc);
    }
}

void movableAxesRepresentation::PrintSelf(ostream &os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);
}
