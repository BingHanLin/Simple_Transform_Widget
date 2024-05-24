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
#include <vtkPolygon.h>
#include <vtkProperty.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkTransform.h>
#include <vtkWindow.h>

#include "modernTransformRepresentation.hpp"

const static double DEFAULT_SIZE = 1.0;
const static double RING_RADIUS = 0.5;
const static double RING_CROSS_SECTION_RADIUS = 0.025;
const static double SCALE_INDICATOR_POS = RING_RADIUS * 1.50;
const static int NUMBER_OF_ARROW_POINTS = 7;

enum class ARROWDIRECTION
{
    X = 0,
    Y = 1,
    Z = 2
};

vtkSmartPointer<vtkActor> generateArrowActor(const ARROWDIRECTION &direction)
{
    /**
    Shape of the arrow:
        /\
       /  \
      /    \
     /_    _\
       |  |
       |__|
    **/

    const std::array<double, NUMBER_OF_ARROW_POINTS> arrowShapeVOffsets{
        0.0, 0.5, 0.5, 1.0, 0.5, 0.5, 0.0};

    const std::array<double, NUMBER_OF_ARROW_POINTS> arrowShapeHOffsets{
        -0.25, -0.25, -0.5, 0.0, 0.5, 0.25, 0.25};

    // Step 1: create the points for the arrow
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(NUMBER_OF_ARROW_POINTS);
    for (int i = 0; i < NUMBER_OF_ARROW_POINTS; i++)
    {
        if (direction == ARROWDIRECTION::X)
        {
            points->SetPoint(i, arrowShapeVOffsets[i], arrowShapeHOffsets[i],
                             0.0);
        }
        else if (direction == ARROWDIRECTION::Y)
        {
            points->SetPoint(i, arrowShapeHOffsets[i], arrowShapeVOffsets[i],
                             0.0);
        }
        else
        {
            points->SetPoint(i, arrowShapeHOffsets[i], 0.0,
                             arrowShapeVOffsets[i]);
        }
    }

    // Step 2: create the polygon for the arrow
    vtkSmartPointer<vtkPolygon> polygon = vtkSmartPointer<vtkPolygon>::New();
    polygon->GetPointIds()->SetNumberOfIds(NUMBER_OF_ARROW_POINTS);
    for (int i = 0; i < NUMBER_OF_ARROW_POINTS; i++)
    {
        polygon->GetPointIds()->SetId(i, i);
    }

    // Step 3: Create a cell array to store the polygon
    vtkSmartPointer<vtkCellArray> polygons =
        vtkSmartPointer<vtkCellArray>::New();
    polygons->InsertNextCell(polygon);

    // Step 4: Create a polydata to store the points and the polygon
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetPolys(polygons);

    // Step 5: Create a mapper and actor
    vtkSmartPointer<vtkPolyDataMapper> mapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(polyData);

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    return actor;
}

vtkStandardNewMacro(modernTransformRepresentation);

modernTransformRepresentation::modernTransformRepresentation()
{
    auto xShapeActor = generateArrowActor(ARROWDIRECTION::X);
    auto yShapeActor = vtkSmartPointer<vtkActor>::New();
    auto zShapeActor = vtkSmartPointer<vtkActor>::New();

    const std::array<vtkSmartPointer<vtkActor>, 3> arrowShapeActors{
        xShapeActor, yShapeActor, zShapeActor};

    auto scaleConeActor = vtkSmartPointer<vtkActor>::New();
    {
        const std::array<double, 3> direction = {1.0, 1.0, 1.0};

        vtkNew<vtkConeSource> coneSource;
        coneSource->SetResolution(30);
        coneSource->SetHeight(RING_CROSS_SECTION_RADIUS * 10);
        coneSource->SetAngle(20);
        coneSource->SetDirection(direction.data());
        coneSource->Update();

        vtkNew<vtkPolyDataMapper> mapper;
        mapper->SetInputConnection(coneSource->GetOutputPort());

        scaleConeActor = vtkSmartPointer<vtkActor>::New();
        scaleConeActor->SetMapper(mapper);
        scaleConeActor->GetProperty()->SetDiffuse(0.8);
        scaleConeActor->GetProperty()->SetSpecular(0.5);
        scaleConeActor->GetProperty()->SetSpecularPower(30.0);
        // coneActor->GetProperty()->SetColor(1.0, 1.0, 1.0);

        vtkSmartPointer<vtkTransform> trans =
            vtkSmartPointer<vtkTransform>::New();
        trans->Translate(SCALE_INDICATOR_POS, SCALE_INDICATOR_POS,
                         SCALE_INDICATOR_POS);
        scaleConeActor->SetUserTransform(trans);
    }

    for (auto i = 0; i < 3; i++)
    {
        rotateActors_[i] = vtkSmartPointer<vtkAssembly>::New();
        // rotateActors_[i]->AddPart(axisRingCircleActors[i]);
        // rotateActors_[i]->AddPart(axisRingConeActors[i]);

        // rotateActors_[i]->SetOrigin(0.0, 0.0, 0.0);

        // vtkSmartPointer<vtkTransform> trans =
        //     vtkSmartPointer<vtkTransform>::New();

        // rotateActors_[i]->SetUserTransform(trans);
    }

    for (auto i = 0; i < 3; i++)
    {
        translateActors_[i] = vtkSmartPointer<vtkAssembly>::New();
        translateActors_[i]->AddPart(arrowShapeActors[i]);

        translateActors_[i]->SetOrigin(0.0, 0.0, 0.0);

        vtkSmartPointer<vtkTransform> trans =
            vtkSmartPointer<vtkTransform>::New();

        translateActors_[i]->SetUserTransform(trans);
    }

    {
        scaleActor_ = vtkSmartPointer<vtkAssembly>::New();
        scaleActor_->AddPart(scaleConeActor);

        scaleActor_->SetOrigin(0.0, 0.0, 0.0);

        vtkSmartPointer<vtkTransform> trans =
            vtkSmartPointer<vtkTransform>::New();

        scaleActor_->SetUserTransform(trans);
    }

    {
        assembleActor_ = vtkSmartPointer<vtkAssembly>::New();

        for (auto i = 0; i < 3; i++)
        {
            // assembleActor_->AddPart(rotateActors_[i]);
            assembleActor_->AddPart(translateActors_[i]);
        }
        assembleActor_->AddPart(scaleActor_);

        assembleActor_->SetOrigin(0.0, 0.0, 0.0);

        vtkSmartPointer<vtkTransform> trans =
            vtkSmartPointer<vtkTransform>::New();

        assembleActor_->SetUserTransform(trans);
    }

    {
        picker_ = vtkSmartPointer<vtkCellPicker>::New();
        picker_->SetTolerance(0.001);
        picker_->AddPickList(assembleActor_);

        picker_->PickFromListOn();
    }

    // Define the point coordinates
    double bounds[6];
    bounds[0] = -DEFAULT_SIZE;
    bounds[1] = DEFAULT_SIZE;
    bounds[2] = -DEFAULT_SIZE;
    bounds[3] = DEFAULT_SIZE;
    bounds[4] = -DEFAULT_SIZE;
    bounds[5] = DEFAULT_SIZE;
    this->PlaceWidget(bounds);
}

modernTransformRepresentation::~modernTransformRepresentation()
{
}

void modernTransformRepresentation::StartWidgetInteraction(double e[2])
{
    auto path = this->GetAssemblyPath(e[0], e[1], 0., picker_);
    if (path != nullptr)
    {
        picker_->GetPickPosition(prevEventWorldPosition_.data());
    }

    prevEventPosition_ = {e[0], e[1], 0.0};
}

void modernTransformRepresentation::WidgetInteraction(double e[2])
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

void modernTransformRepresentation::PlaceWidget(double bounds[6])
{
    // bounds[6]: xmin, xmax, ymin, ymax, zmin
    const std::array<double, 3> minPoint = {bounds[0], bounds[2], bounds[4]};
    const std::array<double, 3> maxPoint = {bounds[1], bounds[3], bounds[5]};

    const auto newDiagonalLengthSquare =
        vtkMath::Distance2BetweenPoints(minPoint.data(), maxPoint.data());
    const auto newDiagonalLength = sqrt(newDiagonalLengthSquare);

    const auto diagonalLength = (RING_RADIUS * 2.0) * sqrt(3);

    const std::array<double, 3> placeScale = {
        newDiagonalLength / diagonalLength, newDiagonalLength / diagonalLength,
        newDiagonalLength / diagonalLength};

    const std::array<double, 3> placeCenter = {
        (minPoint[0] + maxPoint[0]) * 0.5, (minPoint[1] + maxPoint[1]) * 0.5,
        (minPoint[2] + maxPoint[2]) * 0.5};

    placeCenter_ = placeCenter;
    placeScale_ = placeScale;

    vtkNew<vtkTransform> trans;
    trans->PostMultiply();

    trans->Scale(placeScale.data());
    trans->Translate(placeCenter.data());

    assembleActor_->SetUserTransform(trans);
}

int modernTransformRepresentation::ComputeInteractionState(
    int x, int y, int vtkNotUsed(modify))
{
    if (!this->Renderer || !this->Renderer->IsInViewport(x, y))
    {
        this->InteractionState = INTERACTIONSTATE::outside;
        return this->InteractionState;
    }

    this->InteractionState = INTERACTIONSTATE::outside;

    auto path = this->GetAssemblyPath(x, y, 0., picker_);
    if (path != nullptr)
    {
        vtkCollectionSimpleIterator simpleIt;
        path->InitTraversal(simpleIt);
        for (int i = 0; i < path->GetNumberOfItems(); i++)
        {
            auto pickedActor = path->GetNextNode(simpleIt)->GetViewProp();

            if (pickedActor == rotateActors_[0])
            {
                this->InteractionState = INTERACTIONSTATE::onXRing;
            }
            else if (pickedActor == rotateActors_[1])
            {
                this->InteractionState = INTERACTIONSTATE::onYRing;
            }
            else if (pickedActor == rotateActors_[2])
            {
                this->InteractionState = INTERACTIONSTATE::onZRing;
            }
            else if (pickedActor == translateActors_[0])
            {
                this->InteractionState = INTERACTIONSTATE::onXArrow;
            }
            else if (pickedActor == translateActors_[1])
            {
                this->InteractionState = INTERACTIONSTATE::onYArrow;
            }
            else if (pickedActor == translateActors_[2])
            {
                this->InteractionState = INTERACTIONSTATE::onZArrow;
            }
            else if (pickedActor == scaleActor_)
            {
                this->InteractionState = INTERACTIONSTATE::onScale;
            }
        }
    }

    return this->InteractionState;
}

void modernTransformRepresentation::GetTransform(vtkTransform *t)
{
    t->Identity();
    t->PostMultiply();

    t->Translate(-placeCenter_[0], -placeCenter_[1], -placeCenter_[2]);
    t->Scale(1.0 / placeScale_[0], 1.0 / placeScale_[1], 1.0 / placeScale_[2]);

    t->Concatenate(assembleActor_->GetUserMatrix());
}

void modernTransformRepresentation::Highlight(int highlight)
{
    vtkNew<vtkPropCollection> propsCollection;
    assembleActor_->GetActors(propsCollection);

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
        else if (state == INTERACTIONSTATE::onScale)
        {
            scaleActor_->GetActors(propsCollection);
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

void modernTransformRepresentation::BuildRepresentation()
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
            const double originPoint[4] = {0.0, 0.0, 0.0, 1.0};

            double normalAxisEndPoint[4] = {0.0, 0.0, 0.0, 1.0};
            {
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
            }

            vtkMatrix4x4 *originMatrix = assembleActor_->GetUserMatrix();

            double tranformedOriginPoint[4];
            originMatrix->MultiplyPoint(originPoint, tranformedOriginPoint);

            double tranformedNormalAxisEndPoint[4];
            originMatrix->MultiplyPoint(normalAxisEndPoint,
                                        tranformedNormalAxisEndPoint);

            double projectedPrevPickedPoint[3];
            {
                const double vec[3] = {
                    prevEventWorldPosition_[0] - tranformedOriginPoint[0],
                    prevEventWorldPosition_[1] - tranformedOriginPoint[1],
                    prevEventWorldPosition_[2] - tranformedOriginPoint[2]};

                double unitNormal[3] = {
                    tranformedNormalAxisEndPoint[0] - tranformedOriginPoint[0],
                    tranformedNormalAxisEndPoint[1] - tranformedOriginPoint[1],
                    tranformedNormalAxisEndPoint[2] - tranformedOriginPoint[2]};

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
                const double vec[3] = {
                    currEventWorldPosition_[0] - tranformedOriginPoint[0],
                    currEventWorldPosition_[1] - tranformedOriginPoint[1],
                    currEventWorldPosition_[2] - tranformedOriginPoint[2]};

                double unitNormal[3] = {
                    tranformedNormalAxisEndPoint[0] - tranformedOriginPoint[0],
                    tranformedNormalAxisEndPoint[1] - tranformedOriginPoint[1],
                    tranformedNormalAxisEndPoint[2] - tranformedOriginPoint[2]};

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
            vtkMath::Cross(projectedPrevPickedVec, projectedCurrPickedVec,
                           rotateAxis);

            const auto dot =
                vtkMath::Dot(projectedPrevPickedVec, projectedCurrPickedVec);
            const double angleRadians = std::acos(dot);

            const double angleDegrees =
                vtkMath::DegreesFromRadians(angleRadians);

            {
                vtkNew<vtkTransform> trans;
                trans->PostMultiply();
                trans->SetMatrix(originMatrix);

                std::array<double, 3> positionValues;
                trans->GetPosition(positionValues.data());

                trans->Translate(-positionValues[0], -positionValues[1],
                                 -positionValues[2]);

                trans->RotateWXYZ(angleDegrees, rotateAxis);

                trans->Translate(positionValues[0], positionValues[1],
                                 positionValues[2]);

                assembleActor_->SetUserTransform(trans);
            }
        }
        else if (this->InteractionState == INTERACTIONSTATE::onXArrow ||
                 this->InteractionState == INTERACTIONSTATE::onYArrow ||
                 this->InteractionState == INTERACTIONSTATE::onZArrow)
        {
            double direction[3];
            {
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
            }

            vtkMatrix4x4 *originMatrix = assembleActor_->GetUserMatrix();
            vtkNew<vtkMatrix4x4> invertedMatrix;
            vtkMatrix4x4::Invert(originMatrix, invertedMatrix);

            auto pos = invertedMatrix->MultiplyDoublePoint(
                prevEventWorldPosition_.data());
            double originPrevPickedWorldPoint[3] = {pos[0], pos[1], pos[2]};

            pos = invertedMatrix->MultiplyDoublePoint(
                currEventWorldPosition_.data());
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

            {
                vtkNew<vtkTransform> trans;
                trans->SetMatrix(originMatrix);
                trans->Translate(direction[0] * projectedDiff[0],
                                 direction[1] * projectedDiff[1],
                                 direction[2] * projectedDiff[2]);

                assembleActor_->SetUserTransform(trans);
            }
        }
        else if (this->InteractionState == INTERACTIONSTATE::onScale)
        {
            const double direction[3] = {1.0, 1.0, 1.0};

            vtkNew<vtkMatrix4x4> originMatrix;
            originMatrix->DeepCopy(assembleActor_->GetUserMatrix());

            {
                vtkNew<vtkMatrix4x4> invertedMatrix;
                vtkMatrix4x4::Invert(originMatrix, invertedMatrix);

                auto pos = invertedMatrix->MultiplyDoublePoint(
                    prevEventWorldPosition_.data());
                const double originPrevPickedWorldPoint[3] = {pos[0], pos[1],
                                                              pos[2]};

                pos = invertedMatrix->MultiplyDoublePoint(
                    currEventWorldPosition_.data());
                const double originCurrPickedWorldPoint[3] = {pos[0], pos[1],
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
                };

                std::array<double, 3> originScaleValues;
                vtkNew<vtkTransform> trans;
                trans->SetMatrix(originMatrix);
                trans->GetScale(originScaleValues.data());

                const double origin[3] = {0.0, 0.0, 0.0};
                const double initialPoint[3] = {
                    SCALE_INDICATOR_POS * originScaleValues[0],
                    SCALE_INDICATOR_POS * originScaleValues[1],
                    SCALE_INDICATOR_POS * originScaleValues[2]};

                const double projectedCurrDiffSquare =
                    vtkMath::Distance2BetweenPoints(origin,
                                                    projectedCurrPickedPoint);
                const double projectedCurrDiff = sqrt(projectedCurrDiffSquare);

                const double projectedPreviffSquare =
                    vtkMath::Distance2BetweenPoints(origin,
                                                    projectedPrevPickedPoint);
                const double projectedPrevDiff = sqrt(projectedPreviffSquare);

                const double ratio = projectedCurrDiff / projectedPrevDiff;

                {
                    vtkNew<vtkTransform> trans;
                    trans->PostMultiply();
                    trans->SetMatrix(originMatrix);

                    std::array<double, 3> positionValues;
                    trans->GetPosition(positionValues.data());

                    trans->Translate(-positionValues[0], -positionValues[1],
                                     -positionValues[2]);

                    trans->Scale(ratio, ratio, ratio);

                    trans->Translate(positionValues[0], positionValues[1],
                                     positionValues[2]);

                    assembleActor_->SetUserTransform(trans);
                }
            }
        }

        this->BuildTime.Modified();
    }
}

void modernTransformRepresentation::GetActors(vtkPropCollection *pc)
{
    pc->AddItem(assembleActor_);
}