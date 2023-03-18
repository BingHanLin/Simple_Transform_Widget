#include <vtkActor.h>
#include <vtkAssembly.h>
#include <vtkAssemblyPath.h>
#include <vtkCallbackCommand.h>
#include <vtkCellPicker.h>
#include <vtkConeSource.h>
#include <vtkCylinderSource.h>
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
    const std::array<std::array<double, 3>, 3> lineDirection = {
        {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};

    const std::array<std::array<double, 3>, 3> coneDirection = {
        {{0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0}}};

    const double ringRadius = 0.5;
    const double ringSectionRadius = 0.025;

    std::array<vtkSmartPointer<vtkActor>, 3> axisRingCircleActors;
    std::array<vtkSmartPointer<vtkActor>, 3> axisRingConeActors;
    for (auto i = 0; i < 3; i++)
    {
        {
            vtkNew<vtkParametricTorus> paramTorus;
            paramTorus->SetRingRadius(ringRadius);
            paramTorus->SetCrossSectionRadius(ringSectionRadius);

            vtkNew<vtkParametricFunctionSource> paramSource;
            paramSource->SetParametricFunction(paramTorus);
            paramSource->SetUResolution(60);
            paramSource->SetVResolution(20);
            paramSource->Update();

            vtkNew<vtkPolyDataMapper> mapper;
            mapper->SetInputConnection(paramSource->GetOutputPort());

            axisRingCircleActors[i] = vtkSmartPointer<vtkActor>::New();
            axisRingCircleActors[i]->SetMapper(mapper);
            axisRingCircleActors[i]->GetProperty()->SetDiffuse(0.8);
            axisRingCircleActors[i]->GetProperty()->SetSpecular(0.5);
            axisRingCircleActors[i]->GetProperty()->SetSpecularPower(30.0);
            axisRingCircleActors[i]->GetProperty()->SetOpacity(1.0);
            // axisRingCircleActors[i]->GetProperty()->SetLineWidth(10);
            // axisRingCircleActors[i]
            // ->GetProperty()
            // ->SetRepresentationToWireframe();
            axisRingCircleActors[i]->GetProperty()->SetColor(
                axisNormalColor_[i][0], axisNormalColor_[i][1],
                axisNormalColor_[i][2]);
        }

        {
            vtkNew<vtkConeSource> coneSource;
            coneSource->SetResolution(30);
            coneSource->SetHeight(ringSectionRadius * 10);
            coneSource->SetAngle(20);
            coneSource->SetDirection(coneDirection[i].data());
            coneSource->Update();

            vtkNew<vtkPolyDataMapper> mapper;
            mapper->SetInputConnection(coneSource->GetOutputPort());

            axisRingConeActors[i] = vtkSmartPointer<vtkActor>::New();
            axisRingConeActors[i]->SetMapper(mapper);
            axisRingConeActors[i]->GetProperty()->SetDiffuse(0.8);
            axisRingConeActors[i]->GetProperty()->SetSpecular(0.5);
            axisRingConeActors[i]->GetProperty()->SetSpecularPower(30.0);
            axisRingConeActors[i]->GetProperty()->SetColor(
                axisNormalColor_[i][0], axisNormalColor_[i][1],
                axisNormalColor_[i][2]);
        }
    }

    std::array<vtkSmartPointer<vtkActor>, 3> axisLineActors;
    for (auto i = 0; i < 3; i++)
    {
        {
            vtkNew<vtkCylinderSource> cylinderSource;
            cylinderSource->SetResolution(20);
            cylinderSource->SetHeight(2 * ringRadius * 1.50);
            cylinderSource->SetRadius(ringSectionRadius);

            vtkNew<vtkPolyDataMapper> mapper;
            mapper->SetInputConnection(cylinderSource->GetOutputPort());

            axisLineActors[i] = vtkSmartPointer<vtkActor>::New();
            axisLineActors[i]->SetMapper(mapper);
            axisLineActors[i]->GetProperty()->SetDiffuse(0.8);
            axisLineActors[i]->GetProperty()->SetSpecular(0.5);
            axisLineActors[i]->GetProperty()->SetSpecularPower(30.0);
            axisLineActors[i]->GetProperty()->SetColor(axisNormalColor_[i][0],
                                                       axisNormalColor_[i][1],
                                                       axisNormalColor_[i][2]);
        }
    }

    {
        vtkSmartPointer<vtkTransform> trans =
            vtkSmartPointer<vtkTransform>::New();
        trans->RotateY(90);
        axisRingCircleActors[0]->SetUserTransform(trans);
    }

    {
        vtkSmartPointer<vtkTransform> trans =
            vtkSmartPointer<vtkTransform>::New();
        trans->RotateX(90);
        axisRingCircleActors[1]->SetUserTransform(trans);
    }

    {
        vtkSmartPointer<vtkTransform> trans =
            vtkSmartPointer<vtkTransform>::New();
        axisRingCircleActors[2]->SetUserTransform(trans);
    }

    {
        vtkSmartPointer<vtkTransform> trans =
            vtkSmartPointer<vtkTransform>::New();
        trans->Translate(0.0, 0.0, -ringRadius);
        axisRingConeActors[0]->SetUserTransform(trans);
    }

    {
        vtkSmartPointer<vtkTransform> trans =
            vtkSmartPointer<vtkTransform>::New();
        trans->Translate(-ringRadius, 0.0, 0.0);
        axisRingConeActors[1]->SetUserTransform(trans);
    }

    {
        vtkSmartPointer<vtkTransform> trans =
            vtkSmartPointer<vtkTransform>::New();
        trans->Translate(0.0, -ringRadius, 0.0);
        axisRingConeActors[2]->SetUserTransform(trans);
    }

    {
        vtkSmartPointer<vtkTransform> trans =
            vtkSmartPointer<vtkTransform>::New();
        trans->RotateZ(90);
        axisLineActors[0]->SetUserTransform(trans);
    }

    {
        vtkSmartPointer<vtkTransform> trans =
            vtkSmartPointer<vtkTransform>::New();
        axisLineActors[1]->SetUserTransform(trans);
    }

    {
        vtkSmartPointer<vtkTransform> trans =
            vtkSmartPointer<vtkTransform>::New();
        trans->RotateX(90);
        axisLineActors[2]->SetUserTransform(trans);
    }

    for (auto i = 0; i < 3; i++)
    {
        axisRingActors_[i] = vtkSmartPointer<vtkAssembly>::New();
        axisRingActors_[i]->AddPart(axisRingCircleActors[i]);
        axisRingActors_[i]->AddPart(axisRingConeActors[i]);

        axisRingActors_[i]->SetOrigin(0.0, 0.0, 0.0);

        vtkSmartPointer<vtkTransform> trans =
            vtkSmartPointer<vtkTransform>::New();

        axisRingActors_[i]->SetUserTransform(trans);
    }

    for (auto i = 0; i < 3; i++)
    {
        axisArrowActors_[i] = vtkSmartPointer<vtkAssembly>::New();
        axisArrowActors_[i]->AddPart(axisLineActors[i]);

        axisArrowActors_[i]->SetOrigin(0.0, 0.0, 0.0);

        vtkSmartPointer<vtkTransform> trans =
            vtkSmartPointer<vtkTransform>::New();

        axisArrowActors_[i]->SetUserTransform(trans);
    }

    {
        picker_ = vtkSmartPointer<vtkCellPicker>::New();
        picker_->SetTolerance(0.001);
        for (auto i = 0; i < 3; i++)
        {
            picker_->AddPickList(axisRingActors_[i]);
            picker_->AddPickList(axisArrowActors_[i]);
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
    std::cout << "WidgetInteraction: " << e[0] << ", " << e[1] << std::endl;

    double prevPickedWorldPoint[4];
    vtkInteractorObserver::ComputeDisplayToWorld(
        this->Renderer, lastEventPosition_[0], lastEventPosition_[1], 0,
        prevPickedWorldPoint);

    double currPickedWorldPoint[4];
    vtkInteractorObserver::ComputeDisplayToWorld(this->Renderer, e[0], e[1], 0,
                                                 currPickedWorldPoint);

    std::cout << "prevPickedWorldPoint: " << prevPickedWorldPoint[0] << ", "
              << prevPickedWorldPoint[1] << ", " << prevPickedWorldPoint[2]
              << std::endl;

    std::cout << "currPickedWorldPoint: " << currPickedWorldPoint[0] << ", "
              << currPickedWorldPoint[1] << ", " << currPickedWorldPoint[2]
              << std::endl;

    if (this->InteractionState == INTERACTIONSTATE::onXRing ||
        this->InteractionState == INTERACTIONSTATE::onYRing ||
        this->InteractionState == INTERACTIONSTATE::onZRing)
    {
        double normal[3];
        if (this->InteractionState == INTERACTIONSTATE::onXRing)
        {
            normal[0] = 1.0;
            normal[1] = 0.0;
            normal[2] = 0.0;
        }
        else if (this->InteractionState == INTERACTIONSTATE::onYRing)
        {
            normal[0] = 0.0;
            normal[1] = 1.0;
            normal[2] = 0.0;
        }
        else
        {
            normal[0] = 0.0;
            normal[1] = 0.0;
            normal[2] = 1.0;
        }

        for (int i = 0; i < 3; i++)
        {
            vtkMatrix4x4 *originMatrix =
                this->axisRingActors_[i]->GetUserMatrix();

            vtkNew<vtkMatrix4x4> newMatrix;
            {
                vtkNew<vtkMatrix4x4> invertedMatrix;
                vtkMatrix4x4::Invert(originMatrix, invertedMatrix);

                auto pos =
                    invertedMatrix->MultiplyDoublePoint(prevPickedWorldPoint);
                double originPrevPickedWorldPoint[3] = {pos[0], pos[1], pos[2]};

                pos = invertedMatrix->MultiplyDoublePoint(currPickedWorldPoint);
                double originCurrPickedWorldPoint[3] = {pos[0], pos[1], pos[2]};

                double projectedPrevPickedPoint[3];
                {
                    const double dist =
                        vtkMath::Dot(originPrevPickedWorldPoint, normal);
                    projectedPrevPickedPoint[0] =
                        originPrevPickedWorldPoint[0] - dist * normal[0];
                    projectedPrevPickedPoint[1] =
                        originPrevPickedWorldPoint[1] - dist * normal[1];
                    projectedPrevPickedPoint[2] =
                        originPrevPickedWorldPoint[2] - dist * normal[2];
                }

                double projectedCurrPickedPoint[3];
                {
                    const double dist =
                        vtkMath::Dot(originCurrPickedWorldPoint, normal);
                    projectedCurrPickedPoint[0] =
                        originCurrPickedWorldPoint[0] - dist * normal[0];
                    projectedCurrPickedPoint[1] =
                        originCurrPickedWorldPoint[1] - dist * normal[1];
                    projectedCurrPickedPoint[2] =
                        originCurrPickedWorldPoint[2] - dist * normal[2];
                }

                vtkMath::Normalize(projectedPrevPickedPoint);
                vtkMath::Normalize(projectedCurrPickedPoint);

                double axis[3];
                vtkMath::Cross(projectedPrevPickedPoint,
                               projectedCurrPickedPoint, axis);

                const double angleRadians = std::acos(vtkMath::Dot(
                    projectedPrevPickedPoint, projectedCurrPickedPoint));

                const double angleDegrees =
                    vtkMath::DegreesFromRadians(angleRadians);

                vtkNew<vtkTransform> trans;
                trans->SetMatrix(originMatrix);
                trans->RotateWXYZ(angleDegrees, axis);
                newMatrix->DeepCopy(trans->GetMatrix());
            }
            this->axisRingActors_[i]->SetUserMatrix(newMatrix);
        }

        for (int i = 0; i < 3; i++)
        {
            vtkMatrix4x4 *originMatrix =
                this->axisArrowActors_[i]->GetUserMatrix();

            vtkNew<vtkMatrix4x4> newMatrix;
            {
                vtkNew<vtkMatrix4x4> invertedMatrix;
                vtkMatrix4x4::Invert(originMatrix, invertedMatrix);

                auto pos =
                    invertedMatrix->MultiplyDoublePoint(prevPickedWorldPoint);
                double originPrevPickedWorldPoint[3] = {pos[0], pos[1], pos[2]};

                pos = invertedMatrix->MultiplyDoublePoint(currPickedWorldPoint);
                double originCurrPickedWorldPoint[3] = {pos[0], pos[1], pos[2]};

                double projectedPrevPickedPoint[3];
                {
                    const double dist =
                        vtkMath::Dot(originPrevPickedWorldPoint, normal);
                    projectedPrevPickedPoint[0] =
                        originPrevPickedWorldPoint[0] - dist * normal[0];
                    projectedPrevPickedPoint[1] =
                        originPrevPickedWorldPoint[1] - dist * normal[1];
                    projectedPrevPickedPoint[2] =
                        originPrevPickedWorldPoint[2] - dist * normal[2];
                }

                double projectedCurrPickedPoint[3];
                {
                    const double dist =
                        vtkMath::Dot(originCurrPickedWorldPoint, normal);
                    projectedCurrPickedPoint[0] =
                        originCurrPickedWorldPoint[0] - dist * normal[0];
                    projectedCurrPickedPoint[1] =
                        originCurrPickedWorldPoint[1] - dist * normal[1];
                    projectedCurrPickedPoint[2] =
                        originCurrPickedWorldPoint[2] - dist * normal[2];
                }

                vtkMath::Normalize(projectedPrevPickedPoint);
                vtkMath::Normalize(projectedCurrPickedPoint);

                double axis[3];
                vtkMath::Cross(projectedPrevPickedPoint,
                               projectedCurrPickedPoint, axis);

                const double angleRadians = std::acos(vtkMath::Dot(
                    projectedPrevPickedPoint, projectedCurrPickedPoint));

                const double angleDegrees =
                    vtkMath::DegreesFromRadians(angleRadians);

                vtkNew<vtkTransform> trans;
                trans->SetMatrix(originMatrix);
                trans->RotateWXYZ(angleDegrees, axis);
                newMatrix->DeepCopy(trans->GetMatrix());
            }
            this->axisArrowActors_[i]->SetUserMatrix(newMatrix);
        }
    }
    else if (this->InteractionState == INTERACTIONSTATE::onXArrow ||
             this->InteractionState == INTERACTIONSTATE::onYArrow ||
             this->InteractionState == INTERACTIONSTATE::onZArrow)
    {
        double direction[3];
        if (this->InteractionState == INTERACTIONSTATE::onXArrow)
        {
            direction[0] = 1.0;
            direction[1] = 0.0;
            direction[2] = 0.0;
        }
        else if (this->InteractionState == INTERACTIONSTATE::onYArrow)
        {
            direction[0] = 0.0;
            direction[1] = 1.0;
            direction[2] = 0.0;
        }
        else
        {
            direction[0] = 0.0;
            direction[1] = 0.0;
            direction[2] = 1.0;
        }

        for (int i = 0; i < 3; i++)
        {
            vtkMatrix4x4 *originMatrix =
                this->axisRingActors_[i]->GetUserMatrix();

            vtkNew<vtkMatrix4x4> newMatrix;
            {
                vtkNew<vtkMatrix4x4> invertedMatrix;
                vtkMatrix4x4::Invert(originMatrix, invertedMatrix);

                auto pos =
                    invertedMatrix->MultiplyDoublePoint(prevPickedWorldPoint);
                double originPrevPickedWorldPoint[3] = {pos[0], pos[1], pos[2]};

                pos = invertedMatrix->MultiplyDoublePoint(currPickedWorldPoint);
                double originCurrPickedWorldPoint[3] = {pos[0], pos[1], pos[2]};

                double projectedPrevPickedPoint[3];
                {
                    const double dist =
                        vtkMath::Dot(originPrevPickedWorldPoint, direction);
                    projectedPrevPickedPoint[0] = dist * direction[0];
                    projectedPrevPickedPoint[1] = dist * direction[1];
                    projectedPrevPickedPoint[2] = dist * direction[2];
                }

                double projectedCurrPickedPoint[3];
                {
                    const double dist =
                        vtkMath::Dot(originCurrPickedWorldPoint, direction);
                    projectedCurrPickedPoint[0] = dist * direction[0];
                    projectedCurrPickedPoint[1] = dist * direction[1];
                    projectedCurrPickedPoint[2] = dist * direction[2];
                }

                double projectedDiff[3];
                vtkMath::Subtract(projectedCurrPickedPoint,
                                  projectedPrevPickedPoint, projectedDiff);

                vtkNew<vtkTransform> trans;
                trans->SetMatrix(originMatrix);
                trans->Translate(direction[0] * projectedDiff[0],
                                 direction[1] * projectedDiff[1],
                                 direction[2] * projectedDiff[2]);
                newMatrix->DeepCopy(trans->GetMatrix());
            }
            this->axisRingActors_[i]->SetUserMatrix(newMatrix);
        }

        for (int i = 0; i < 3; i++)
        {
            vtkMatrix4x4 *originMatrix =
                this->axisArrowActors_[i]->GetUserMatrix();

            vtkNew<vtkMatrix4x4> newMatrix;
            {
                vtkNew<vtkMatrix4x4> invertedMatrix;
                vtkMatrix4x4::Invert(originMatrix, invertedMatrix);

                auto pos =
                    invertedMatrix->MultiplyDoublePoint(prevPickedWorldPoint);
                double originPrevPickedWorldPoint[3] = {pos[0], pos[1], pos[2]};

                pos = invertedMatrix->MultiplyDoublePoint(currPickedWorldPoint);
                double originCurrPickedWorldPoint[3] = {pos[0], pos[1], pos[2]};

                double projectedPrevPickedPoint[3];
                {
                    const double dist =
                        vtkMath::Dot(originPrevPickedWorldPoint, direction);
                    projectedPrevPickedPoint[0] = dist * direction[0];
                    projectedPrevPickedPoint[1] = dist * direction[1];
                    projectedPrevPickedPoint[2] = dist * direction[2];
                }

                double projectedCurrPickedPoint[3];
                {
                    const double dist =
                        vtkMath::Dot(originCurrPickedWorldPoint, direction);
                    projectedCurrPickedPoint[0] = dist * direction[0];
                    projectedCurrPickedPoint[1] = dist * direction[1];
                    projectedCurrPickedPoint[2] = dist * direction[2];
                }

                double projectedDiff[3];
                vtkMath::Subtract(projectedCurrPickedPoint,
                                  projectedPrevPickedPoint, projectedDiff);

                vtkNew<vtkTransform> trans;
                trans->SetMatrix(originMatrix);
                trans->Translate(direction[0] * projectedDiff[0],
                                 direction[1] * projectedDiff[1],
                                 direction[2] * projectedDiff[2]);
                newMatrix->DeepCopy(trans->GetMatrix());
            }
            this->axisArrowActors_[i]->SetUserMatrix(newMatrix);
        }
    }

    lastEventPosition_ = {e[0], e[1], 0.0};

    this->Modified();
    this->BuildRepresentation();
}

void movableAxesRepresentation::PlaceWidget(double bds[6])
{
}

int movableAxesRepresentation::ComputeInteractionState(int x, int y,
                                                       int vtkNotUsed(modify))
{
    for (auto i = 0; i < 3; i++)
    {
        vtkNew<vtkPropCollection> propsCollection;
        axisRingActors_[i]->GetActors(propsCollection);

        vtkCollectionSimpleIterator sIt;
        propsCollection->InitTraversal(sIt);
        const int nProps = propsCollection->GetNumberOfItems();
        for (int i = 0; i < nProps; i++)
        {
            vtkActor *actor =
                vtkActor::SafeDownCast(propsCollection->GetNextProp(sIt));
            if (actor != nullptr)
            {
                actor->GetProperty()->SetDiffuse(0.8);
                actor->GetProperty()->SetSpecular(0.5);
            }
        }
    }

    for (auto i = 0; i < 3; i++)
    {
        vtkNew<vtkPropCollection> propsCollection;
        axisArrowActors_[i]->GetActors(propsCollection);

        vtkCollectionSimpleIterator sIt;
        propsCollection->InitTraversal(sIt);
        const int nProps = propsCollection->GetNumberOfItems();
        for (int i = 0; i < nProps; i++)
        {
            vtkActor *actor =
                vtkActor::SafeDownCast(propsCollection->GetNextProp(sIt));
            if (actor != nullptr)
            {
                actor->GetProperty()->SetDiffuse(0.8);
                actor->GetProperty()->SetSpecular(0.5);
            }
        }
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
            reinterpret_cast<vtkProp3D *>(path->GetFirstNode()->GetViewProp());

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
        else if (currActor_ == axisArrowActors_[0])
        {
            this->InteractionState = INTERACTIONSTATE::onXArrow;
        }
        else if (currActor_ == axisArrowActors_[1])
        {
            this->InteractionState = INTERACTIONSTATE::onYArrow;
        }
        else if (currActor_ == axisArrowActors_[2])
        {
            this->InteractionState = INTERACTIONSTATE::onZArrow;
        }

        vtkNew<vtkPropCollection> propsCollection;
        currActor_->GetActors(propsCollection);

        vtkCollectionSimpleIterator sIt;
        propsCollection->InitTraversal(sIt);
        const int nProps = propsCollection->GetNumberOfItems();
        for (int i = 0; i < nProps; i++)
        {
            vtkActor *actor =
                vtkActor::SafeDownCast(propsCollection->GetNextProp(sIt));
            if (actor != nullptr)
            {
                actor->GetProperty()->SetDiffuse(1.0);
                actor->GetProperty()->SetSpecular(0.0);
            }
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
        pc->AddItem(axisRingActors_[i]);
        pc->AddItem(axisArrowActors_[i]);
    }
}

void movableAxesRepresentation::PrintSelf(ostream &os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);
}
