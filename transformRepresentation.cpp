#include <vtkActor.h>
#include <vtkAssembly.h>
#include <vtkAssemblyPath.h>
#include <vtkCallbackCommand.h>
#include <vtkCamera.h>
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

#include "transformRepresentation.hpp"

vtkStandardNewMacro(transformRepresentation);

transformRepresentation::transformRepresentation()
{
    vtkSmartPointer<vtkTransform> trans = vtkSmartPointer<vtkTransform>::New();
    trans->Identity();

    dummyActor_ = vtkSmartPointer<vtkAssembly>::New();
    dummyActor_->SetUserTransform(trans);

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
    std::array<vtkSmartPointer<vtkActor>, 3> startConeActors;
    std::array<vtkSmartPointer<vtkActor>, 3> endConeActors;
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
        {
            vtkNew<vtkConeSource> startConeSource;
            startConeSource->SetResolution(20);
            startConeSource->SetHeight(ringSectionRadius * 10);
            startConeSource->SetAngle(20);
            startConeSource->SetCenter(lineDirection[i][0] * ringRadius * 1.50,
                                       lineDirection[i][1] * ringRadius * 1.50,
                                       lineDirection[i][2] * ringRadius * 1.50);
            startConeSource->SetDirection(
                lineDirection[i][0], lineDirection[i][1], lineDirection[i][2]);
            startConeSource->Update();

            vtkNew<vtkPolyDataMapper> mapper;
            mapper->SetInputConnection(startConeSource->GetOutputPort());

            startConeActors[i] = vtkSmartPointer<vtkActor>::New();
            startConeActors[i]->SetMapper(mapper);
            startConeActors[i]->GetProperty()->SetDiffuse(0.8);
            startConeActors[i]->GetProperty()->SetSpecular(0.5);
            startConeActors[i]->GetProperty()->SetSpecularPower(30.0);
            startConeActors[i]->GetProperty()->SetColor(axisNormalColor_[i][0],
                                                        axisNormalColor_[i][1],
                                                        axisNormalColor_[i][2]);
        }

        {
            vtkNew<vtkConeSource> endConeSource;
            endConeSource->SetResolution(20);
            endConeSource->SetHeight(ringSectionRadius * 10);
            endConeSource->SetAngle(20);
            endConeSource->SetCenter(-lineDirection[i][0] * ringRadius * 1.50,
                                     -lineDirection[i][1] * ringRadius * 1.50,
                                     -lineDirection[i][2] * ringRadius * 1.50);
            endConeSource->SetDirection(-lineDirection[i][0],
                                        -lineDirection[i][1],
                                        -lineDirection[i][2]);
            endConeSource->Update();

            vtkNew<vtkPolyDataMapper> mapper;
            mapper->SetInputConnection(endConeSource->GetOutputPort());

            endConeActors[i] = vtkSmartPointer<vtkActor>::New();
            endConeActors[i]->SetMapper(mapper);
            endConeActors[i]->GetProperty()->SetDiffuse(0.8);
            endConeActors[i]->GetProperty()->SetSpecular(0.5);
            endConeActors[i]->GetProperty()->SetSpecularPower(30.0);
            endConeActors[i]->GetProperty()->SetColor(axisNormalColor_[i][0],
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
        rotateActors_[i] = vtkSmartPointer<vtkAssembly>::New();
        rotateActors_[i]->AddPart(axisRingCircleActors[i]);
        rotateActors_[i]->AddPart(axisRingConeActors[i]);

        rotateActors_[i]->SetOrigin(0.0, 0.0, 0.0);

        vtkSmartPointer<vtkTransform> trans =
            vtkSmartPointer<vtkTransform>::New();

        rotateActors_[i]->SetUserTransform(trans);
    }

    for (auto i = 0; i < 3; i++)
    {
        translateActors_[i] = vtkSmartPointer<vtkAssembly>::New();
        translateActors_[i]->AddPart(axisLineActors[i]);
        translateActors_[i]->AddPart(startConeActors[i]);
        translateActors_[i]->AddPart(endConeActors[i]);

        translateActors_[i]->SetOrigin(0.0, 0.0, 0.0);

        vtkSmartPointer<vtkTransform> trans =
            vtkSmartPointer<vtkTransform>::New();

        translateActors_[i]->SetUserTransform(trans);
    }

    {
        picker_ = vtkSmartPointer<vtkCellPicker>::New();
        picker_->SetTolerance(0.001);
        for (auto i = 0; i < 3; i++)
        {
            picker_->AddPickList(rotateActors_[i]);
            picker_->AddPickList(translateActors_[i]);
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

transformRepresentation::~transformRepresentation()
{
}

void transformRepresentation::StartWidgetInteraction(double e[2])
{
    auto path = this->GetAssemblyPath(e[0], e[1], 0., picker_);
    if (path != nullptr)
    {
        picker_->GetPickPosition(prevEventWorldPosition_.data());
    }

    prevEventPosition_ = {e[0], e[1], 0.0};

    this->ComputeInteractionState(static_cast<int>(e[0]),
                                  static_cast<int>(e[1]), 0);
}

void transformRepresentation::WidgetInteraction(double e[2])
{
    // refer to vtkBoxWidget.cxx
    double focalPoint[4];
    vtkInteractorObserver::ComputeWorldToDisplay(
        this->Renderer, prevEventWorldPosition_[0], prevEventWorldPosition_[1],
        prevEventWorldPosition_[2], focalPoint);

    vtkInteractorObserver::ComputeDisplayToWorld(
        this->Renderer, prevEventPosition_[0], prevEventPosition_[1],
        focalPoint[2], prevEventWorldPosition_.data());

    vtkInteractorObserver::ComputeDisplayToWorld(
        this->Renderer, e[0], e[1], focalPoint[2],
        currEventWorldPosition_.data());

    prevEventPosition_ = {e[0], e[1], 0.0};

    this->BuildRepresentation();

    this->Modified();
}

void transformRepresentation::PlaceWidget(double bounds[6])
{
    // xmin,xmax,ymin,ymax,zmin,
    const std::array<double, 3> minPoint = {bounds[0], bounds[2], bounds[4]};
    const std::array<double, 3> maxPoint = {bounds[1], bounds[3], bounds[5]};

    const auto diagonalLength =
        vtkMath::Distance2BetweenPoints(minPoint.data(), maxPoint.data()) * 0.5;

    for (auto actors :
         {rotateActors_[0], rotateActors_[1], rotateActors_[2],
          translateActors_[0], translateActors_[1], translateActors_[2]})
    {
        vtkNew<vtkTransform> trans;
        trans->Scale(diagonalLength / 1.0, diagonalLength / 1.0,
                     diagonalLength / 1.0);

        vtkNew<vtkMatrix4x4> newMatrix;
        newMatrix->DeepCopy(trans->GetMatrix());

        actors->SetUserMatrix(newMatrix);
    }
}

int transformRepresentation::ComputeInteractionState(int x, int y,
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
            reinterpret_cast<vtkProp3D *>(path->GetFirstNode()->GetViewProp());

        if (currActor_ == rotateActors_[0])
        {
            this->InteractionState = INTERACTIONSTATE::onXRing;
        }
        else if (currActor_ == rotateActors_[1])
        {
            this->InteractionState = INTERACTIONSTATE::onYRing;
        }
        else if (currActor_ == rotateActors_[2])
        {
            this->InteractionState = INTERACTIONSTATE::onZRing;
        }
        else if (currActor_ == translateActors_[0])
        {
            this->InteractionState = INTERACTIONSTATE::onXArrow;
        }
        else if (currActor_ == translateActors_[1])
        {
            this->InteractionState = INTERACTIONSTATE::onYArrow;
        }
        else if (currActor_ == translateActors_[2])
        {
            this->InteractionState = INTERACTIONSTATE::onZArrow;
        }
    }
    else
    {
        this->InteractionState = INTERACTIONSTATE::outside;
    }

    return this->InteractionState;
}

void transformRepresentation::GetTransform(vtkTransform *t)
{
    t->SetMatrix(dummyActor_->GetUserMatrix());
}

void transformRepresentation::Highlight(int highlight)
{
    for (auto actors :
         {rotateActors_[0], rotateActors_[1], rotateActors_[2],
          translateActors_[0], translateActors_[1], translateActors_[2]})
    {
        vtkNew<vtkPropCollection> propsCollection;
        actors->GetActors(propsCollection);

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

    if (highlight)
    {
        const auto state = this->InteractionState;
        vtkNew<vtkPropCollection> propsCollection;

        if (state == INTERACTIONSTATE::onXRing)
        {
            rotateActors_[0]->GetActors(propsCollection);
        }
        else if (state == INTERACTIONSTATE::onYRing)
        {
            rotateActors_[1]->GetActors(propsCollection);
        }
        else if (state == INTERACTIONSTATE::onZRing)
        {
            rotateActors_[2]->GetActors(propsCollection);
        }
        else if (state == INTERACTIONSTATE::onXArrow)
        {
            translateActors_[0]->GetActors(propsCollection);
        }
        else if (state == INTERACTIONSTATE::onYArrow)
        {
            translateActors_[1]->GetActors(propsCollection);
        }
        else if (state == INTERACTIONSTATE::onZArrow)
        {
            translateActors_[2]->GetActors(propsCollection);
        }

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
}

void transformRepresentation::BuildRepresentation()
{
    if (this->GetMTime() > this->BuildTime ||
        (this->Renderer && this->Renderer->GetVTKWindow() &&
         (this->Renderer->GetVTKWindow()->GetMTime() > this->BuildTime ||
          this->Renderer->GetActiveCamera()->GetMTime() > this->BuildTime)))
    {
        if (this->InteractionState == INTERACTIONSTATE::onXRing ||
            this->InteractionState == INTERACTIONSTATE::onYRing ||
            this->InteractionState == INTERACTIONSTATE::onZRing)
        {
            double originPoint[4] = {0.0, 0.0, 0.0, 1.0};

            double normalAxisEndPoint[4] = {0.0, 0.0, 0.0, 1.0};
            if (this->InteractionState == INTERACTIONSTATE::onXRing)
            {
                normalAxisEndPoint[0] = 1.0;
                normalAxisEndPoint[1] = 0.0;
                normalAxisEndPoint[2] = 0.0;
            }
            else if (this->InteractionState == INTERACTIONSTATE::onYRing)
            {
                normalAxisEndPoint[0] = 0.0;
                normalAxisEndPoint[1] = 1.0;
                normalAxisEndPoint[2] = 0.0;
            }
            else
            {
                normalAxisEndPoint[0] = 0.0;
                normalAxisEndPoint[1] = 0.0;
                normalAxisEndPoint[2] = 1.0;
            }

            for (auto actors :
                 {rotateActors_[0], rotateActors_[1], rotateActors_[2],
                  translateActors_[0], translateActors_[1], translateActors_[2],
                  dummyActor_})
            {
                vtkMatrix4x4 *actorTransformMatrix = actors->GetUserMatrix();

                vtkNew<vtkMatrix4x4> newMatrix;
                {
                    double tranformedOriginPoint[4];
                    actorTransformMatrix->MultiplyPoint(originPoint,
                                                        tranformedOriginPoint);

                    double tranformedNormalAxisEndPoint[4];
                    actorTransformMatrix->MultiplyPoint(
                        normalAxisEndPoint, tranformedNormalAxisEndPoint);

                    double projectedPrevPickedPoint[3];
                    {
                        const double vec[3] = {prevEventWorldPosition_[0] -
                                                   tranformedOriginPoint[0],
                                               prevEventWorldPosition_[1] -
                                                   tranformedOriginPoint[1],
                                               prevEventWorldPosition_[2] -
                                                   tranformedOriginPoint[2]};

                        double unitNormal[3] = {
                            tranformedNormalAxisEndPoint[0] -
                                tranformedOriginPoint[0],
                            tranformedNormalAxisEndPoint[1] -
                                tranformedOriginPoint[1],
                            tranformedNormalAxisEndPoint[2] -
                                tranformedOriginPoint[2]};

                        vtkMath::Normalize(unitNormal);

                        const double dist = vtkMath::Dot(vec, unitNormal);

                        projectedPrevPickedPoint[0] =
                            prevEventWorldPosition_[0] - dist * unitNormal[0];
                        projectedPrevPickedPoint[1] =
                            prevEventWorldPosition_[1] - dist * unitNormal[1];
                        projectedPrevPickedPoint[2] =
                            prevEventWorldPosition_[2] - dist * unitNormal[2];
                    }

                    double projectedCurrPickedPoint[3];
                    {
                        const double vec[3] = {currEventWorldPosition_[0] -
                                                   tranformedOriginPoint[0],
                                               currEventWorldPosition_[1] -
                                                   tranformedOriginPoint[1],
                                               currEventWorldPosition_[2] -
                                                   tranformedOriginPoint[2]};

                        double unitNormal[3] = {
                            tranformedNormalAxisEndPoint[0] -
                                tranformedOriginPoint[0],
                            tranformedNormalAxisEndPoint[1] -
                                tranformedOriginPoint[1],
                            tranformedNormalAxisEndPoint[2] -
                                tranformedOriginPoint[2]};

                        vtkMath::Normalize(unitNormal);

                        const double dist = vtkMath::Dot(vec, unitNormal);
                        projectedCurrPickedPoint[0] =
                            currEventWorldPosition_[0] - dist * unitNormal[0];
                        projectedCurrPickedPoint[1] =
                            currEventWorldPosition_[1] - dist * unitNormal[1];
                        projectedCurrPickedPoint[2] =
                            currEventWorldPosition_[2] - dist * unitNormal[2];
                    }

                    double projectedPrevPickedVec[3] = {
                        projectedPrevPickedPoint[0] - tranformedOriginPoint[0],
                        projectedPrevPickedPoint[1] - tranformedOriginPoint[1],
                        projectedPrevPickedPoint[2] - tranformedOriginPoint[2]};

                    double projectedCurrPickedVec[3] = {
                        projectedCurrPickedPoint[0] - tranformedOriginPoint[0],
                        projectedCurrPickedPoint[1] - tranformedOriginPoint[1],
                        projectedCurrPickedPoint[2] - tranformedOriginPoint[2]};

                    vtkMath::Normalize(projectedPrevPickedVec);
                    vtkMath::Normalize(projectedCurrPickedVec);

                    double rotateAxis[3];
                    vtkMath::Cross(projectedPrevPickedVec,
                                   projectedCurrPickedVec, rotateAxis);

                    const auto dot = vtkMath::Dot(projectedPrevPickedVec,
                                                  projectedCurrPickedVec);
                    const double angleRadians = std::acos(dot);

                    const double angleDegrees =
                        vtkMath::DegreesFromRadians(angleRadians);

                    {
                        const double translationX =
                            actorTransformMatrix->GetElement(0, 3);
                        const double translationY =
                            actorTransformMatrix->GetElement(1, 3);
                        const double translationZ =
                            actorTransformMatrix->GetElement(2, 3);

                        vtkNew<vtkTransform> trans;
                        trans->PostMultiply();
                        trans->SetMatrix(actorTransformMatrix);
                        trans->Translate(-translationX, -translationY,
                                         -translationZ);
                        trans->RotateWXYZ(angleDegrees, rotateAxis);
                        trans->Translate(translationX, translationY,
                                         translationZ);

                        newMatrix->DeepCopy(trans->GetMatrix());
                    }
                }
                actors->SetUserMatrix(newMatrix);
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

            for (vtkProp3D *actors :
                 {rotateActors_[0], rotateActors_[1], rotateActors_[2],
                  translateActors_[0], translateActors_[1], translateActors_[2],
                  dummyActor_})
            {
                vtkMatrix4x4 *originMatrix = actors->GetUserMatrix();

                vtkNew<vtkMatrix4x4> newMatrix;
                {
                    vtkNew<vtkMatrix4x4> invertedMatrix;
                    vtkMatrix4x4::Invert(originMatrix, invertedMatrix);

                    auto pos = invertedMatrix->MultiplyDoublePoint(
                        prevEventWorldPosition_.data());
                    double originPrevPickedWorldPoint[3] = {pos[0], pos[1],
                                                            pos[2]};

                    pos = invertedMatrix->MultiplyDoublePoint(
                        currEventWorldPosition_.data());
                    double originCurrPickedWorldPoint[3] = {pos[0], pos[1],
                                                            pos[2]};

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

                    {
                        vtkNew<vtkTransform> trans;
                        trans->SetMatrix(originMatrix);
                        trans->Translate(direction[0] * projectedDiff[0],
                                         direction[1] * projectedDiff[1],
                                         direction[2] * projectedDiff[2]);
                        newMatrix->DeepCopy(trans->GetMatrix());
                    }
                }
                actors->SetUserMatrix(newMatrix);
            }
        }

        this->BuildTime.Modified();
    }
}

void transformRepresentation::GetActors(vtkPropCollection *pc)
{
    for (auto i = 0; i < 3; i++)
    {
        pc->AddItem(rotateActors_[i]);
        pc->AddItem(translateActors_[i]);
    }
}
