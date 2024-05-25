#include <vector>
#include <vtkActor.h>
#include <vtkAssembly.h>
#include <vtkAssemblyPath.h>
#include <vtkCallbackCommand.h>
#include <vtkCamera.h>
#include <vtkCellPicker.h>
#include <vtkConeSource.h>
#include <vtkCylinderSource.h>
#include <vtkInteractorObserver.h>
#include <vtkNamedColors.h>
#include <vtkObjectFactory.h>
#include <vtkParametricFunctionSource.h>
#include <vtkParametricTorus.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyLine.h>
#include <vtkPolygon.h>
#include <vtkProperty.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkTransform.h>
#include <vtkTriangleFilter.h>
#include <vtkWindow.h>

#include "modernTransformRepresentation.hpp"

vtkStandardNewMacro(modernTransformRepresentation);

enum class ARROWDIRECTION
{
    X = 0,
    Y = 1,
    Z = 2
};

const double ROTATE_ARROW_LEFT_ANGLE = 135.0;
const double ROTATE_ARROW_RIGHT_ANGLE = 45.0;
const int NUMBER_OF_ROTATE_ARROW_ARC_POINTS = 10;
const int NUMBER_OF_ROTATE_ARROW_POINTS =
    NUMBER_OF_ROTATE_ARROW_ARC_POINTS * 2 + 6;

const static double ARROW_SHAFT_WIDTH = 0.1;
const static double ARROW_HEAD_WIDTH = 0.2;
const static double ARROW_HEAD_LENGTH = 0.1;

const static double ROTATE_ARROW_RADIUS = 0.3;

const static double DEFAULT_SIZE = 1.0;
const static double RING_RADIUS = 0.5;
const static double RING_CROSS_SECTION_RADIUS = 0.025;
const static double SCALE_INDICATOR_POS = RING_RADIUS * 1.50;
const static int NUMBER_OF_ARROW_POINTS = 7;

/**
Shape of the arrow:
    /\
   /  \
  /    \
 /_    _\
   |  |
   |__|
**/

const static std::vector<double> TRANSLATE_ARROW_HORIZONTAL_OFFSETS{
    -0.1, -0.1, -0.2, 0.0, 0.2, 0.1, 0.1};

const static std::vector<double> TRANSLATE_ARROW_VERTICAL_OFFSETS{
    0.5, 0.75, 0.75, 1.0, 0.75, 0.75, 0.5};

const static std::vector<double> SCALE_ARROW_HORIZONTAL_OFFSETS{
    0.0, -0.1, -0.05, -0.05, -0.1, 0.0, 0.1, 0.05, 0.05, 0.1};

const static std::vector<double> SCALE_ARROW_VERTICAL_OFFSETS{
    0.6, 0.7, 0.7, 0.9, 0.9, 1.0, 0.9, 0.9, 0.7, 0.7};

void getRotateArrowPoints(std::vector<double> &horizontalOffsets,
                          std::vector<double> &verticalOffsets)
{
    horizontalOffsets.resize(NUMBER_OF_ROTATE_ARROW_POINTS);
    verticalOffsets.resize(NUMBER_OF_ROTATE_ARROW_POINTS);

    const auto headRadians = ARROW_HEAD_LENGTH / ROTATE_ARROW_RADIUS;

    const double resolutionAngle =
        (ROTATE_ARROW_LEFT_ANGLE - ROTATE_ARROW_RIGHT_ANGLE) /
        (NUMBER_OF_ROTATE_ARROW_ARC_POINTS - 1);

    int indexCounter = 0;

    {
        const double angleRadians =
            vtkMath::RadiansFromDegrees(ROTATE_ARROW_LEFT_ANGLE) + headRadians;
        const double h = std::cos(angleRadians) * ROTATE_ARROW_RADIUS;
        const double v = std::sin(angleRadians) * ROTATE_ARROW_RADIUS;
        horizontalOffsets[indexCounter] = h;
        verticalOffsets[indexCounter] = v;
        indexCounter++;
    }

    {
        const double angleRadians =
            vtkMath::RadiansFromDegrees(ROTATE_ARROW_LEFT_ANGLE);
        const double h = std::cos(angleRadians) *
                         (ROTATE_ARROW_RADIUS + ARROW_HEAD_WIDTH * 0.5);
        const double v = std::sin(angleRadians) *
                         (ROTATE_ARROW_RADIUS + ARROW_HEAD_WIDTH * 0.5);

        horizontalOffsets[indexCounter] = h;
        verticalOffsets[indexCounter] = v;
        indexCounter++;
    }

    for (double angle = ROTATE_ARROW_LEFT_ANGLE;
         angle >= ROTATE_ARROW_RIGHT_ANGLE; angle -= resolutionAngle)
    {
        const double angleRadians = vtkMath::RadiansFromDegrees(angle);
        const double h = std::cos(angleRadians) *
                         (ROTATE_ARROW_RADIUS + ARROW_SHAFT_WIDTH * 0.5);
        const double v = std::sin(angleRadians) *
                         (ROTATE_ARROW_RADIUS + ARROW_SHAFT_WIDTH * 0.5);

        horizontalOffsets[indexCounter] = h;
        verticalOffsets[indexCounter] = v;
        indexCounter++;
    }

    {
        const double angleRadians =
            vtkMath::RadiansFromDegrees(ROTATE_ARROW_RIGHT_ANGLE);
        const double h = std::cos(angleRadians) *
                         (ROTATE_ARROW_RADIUS + ARROW_HEAD_WIDTH * 0.5);
        const double v = std::sin(angleRadians) *
                         (ROTATE_ARROW_RADIUS + ARROW_HEAD_WIDTH * 0.5);

        horizontalOffsets[indexCounter] = h;
        verticalOffsets[indexCounter] = v;
        indexCounter++;
    }

    {
        const double angleRadians =
            vtkMath::RadiansFromDegrees(ROTATE_ARROW_RIGHT_ANGLE) - headRadians;
        const double h = std::cos(angleRadians) * ROTATE_ARROW_RADIUS;
        const double v = std::sin(angleRadians) * ROTATE_ARROW_RADIUS;
        horizontalOffsets[indexCounter] = h;
        verticalOffsets[indexCounter] = v;
        indexCounter++;
    }

    {
        const double angleRadians =
            vtkMath::RadiansFromDegrees(ROTATE_ARROW_RIGHT_ANGLE);
        const double h = std::cos(angleRadians) *
                         (ROTATE_ARROW_RADIUS - ARROW_HEAD_WIDTH * 0.5);
        const double v = std::sin(angleRadians) *
                         (ROTATE_ARROW_RADIUS - ARROW_HEAD_WIDTH * 0.5);

        horizontalOffsets[indexCounter] = h;
        verticalOffsets[indexCounter] = v;
        indexCounter++;
    }

    for (double angle = ROTATE_ARROW_RIGHT_ANGLE;
         angle <= ROTATE_ARROW_LEFT_ANGLE; angle += resolutionAngle)
    {
        const double angleRadians = vtkMath::RadiansFromDegrees(angle);
        const double h = std::cos(angleRadians) *
                         (ROTATE_ARROW_RADIUS - ARROW_SHAFT_WIDTH * 0.5);
        const double v = std::sin(angleRadians) *
                         (ROTATE_ARROW_RADIUS - ARROW_SHAFT_WIDTH * 0.5);

        horizontalOffsets[indexCounter] = h;
        verticalOffsets[indexCounter] = v;
        indexCounter++;
    }

    {
        const double angleRadians =
            vtkMath::RadiansFromDegrees(ROTATE_ARROW_LEFT_ANGLE);
        const double h = std::cos(angleRadians) *
                         (ROTATE_ARROW_RADIUS - ARROW_HEAD_WIDTH * 0.5);
        const double v = std::sin(angleRadians) *
                         (ROTATE_ARROW_RADIUS - ARROW_HEAD_WIDTH * 0.5);

        horizontalOffsets[indexCounter] = h;
        verticalOffsets[indexCounter] = v;
    }
}

vtkSmartPointer<vtkActor> getShapeActor(
    const std::vector<double> &horizontalOffsets,
    const std::vector<double> &verticalOffsets, const ARROWDIRECTION &direction)
{
    if (horizontalOffsets.size() != verticalOffsets.size())
    {
        return nullptr;
    }

    // Step 1: Create the points for the arrow
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(horizontalOffsets.size());
    for (int i = 0; i < points->GetNumberOfPoints(); i++)
    {
        if (direction == ARROWDIRECTION::X)
        {
            points->SetPoint(i, verticalOffsets[i], horizontalOffsets[i], 0.0);
        }
        else if (direction == ARROWDIRECTION::Y)
        {
            points->SetPoint(i, 0.0, verticalOffsets[i], horizontalOffsets[i]);
        }
        else
        {
            points->SetPoint(i, horizontalOffsets[i], 0.0, verticalOffsets[i]);
        }
    }

    // Step 2: Create the polygon for the arrow
    vtkSmartPointer<vtkPolygon> polygon = vtkSmartPointer<vtkPolygon>::New();
    polygon->GetPointIds()->SetNumberOfIds(points->GetNumberOfPoints());
    for (int i = 0; i < polygon->GetNumberOfPoints(); i++)
    {
        polygon->GetPointIds()->SetId(i, i);
    }

    // Step 3: Create a cell array to store the polygon
    vtkSmartPointer<vtkCellArray> polygonCell =
        vtkSmartPointer<vtkCellArray>::New();
    polygonCell->InsertNextCell(polygon);

    // Step 4: Create a polydata to store the points and the polygon
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetPolys(polygonCell);

    // Step 5: Convert input polygonCell and strips to triangles
    vtkSmartPointer<vtkTriangleFilter> geometryFilter =
        vtkSmartPointer<vtkTriangleFilter>::New();
    geometryFilter->SetInputData(polyData);
    geometryFilter->Update();

    // Step 6: Create a mapper and actor
    vtkSmartPointer<vtkPolyDataMapper> mapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(geometryFilter->GetOutput());

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    vtkNew<vtkNamedColors> colors;
    actor->GetProperty()->SetColor(colors->GetColor3d("White").GetData());

    return actor;
}

vtkSmartPointer<vtkActor> getOutlineActor(
    const std::vector<double> &horizontalOffsets,
    const std::vector<double> &verticalOffsets, const ARROWDIRECTION &direction)
{
    if (horizontalOffsets.size() != verticalOffsets.size())
    {
        return nullptr;
    }

    // Step 1: Create the points for the arrow
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(horizontalOffsets.size() + 1);
    for (int i = 0; i < points->GetNumberOfPoints() - 1; i++)
    {
        if (direction == ARROWDIRECTION::X)
        {
            points->SetPoint(i, verticalOffsets[i], horizontalOffsets[i], 0.0);
        }
        else if (direction == ARROWDIRECTION::Y)
        {
            points->SetPoint(i, 0.0, verticalOffsets[i], horizontalOffsets[i]);
        }
        else
        {
            points->SetPoint(i, horizontalOffsets[i], 0.0, verticalOffsets[i]);
        }
    }

    if (direction == ARROWDIRECTION::X)
    {
        points->SetPoint(points->GetNumberOfPoints() - 1, verticalOffsets[0],
                         horizontalOffsets[0], 0.0);
    }
    else if (direction == ARROWDIRECTION::Y)
    {
        points->SetPoint(points->GetNumberOfPoints() - 1, 0.0,
                         verticalOffsets[0], horizontalOffsets[0]);
    }
    else
    {
        points->SetPoint(points->GetNumberOfPoints() - 1, horizontalOffsets[0],
                         0.0, verticalOffsets[0]);
    }

    // Step 2: Create the polygon for the arrow
    vtkSmartPointer<vtkPolyLine> polyLine = vtkSmartPointer<vtkPolyLine>::New();
    polyLine->GetPointIds()->SetNumberOfIds(points->GetNumberOfPoints());
    for (int i = 0; i < polyLine->GetNumberOfPoints(); i++)
    {
        polyLine->GetPointIds()->SetId(i, i);
    }

    // Step 3: Create a cell array to store the polygon
    vtkSmartPointer<vtkCellArray> polygonCell =
        vtkSmartPointer<vtkCellArray>::New();
    polygonCell->InsertNextCell(polyLine);

    // Step 4: Create a polydata to store the points and the polygon
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetLines(polygonCell);

    // Step 6: Create a mapper and actor
    vtkSmartPointer<vtkPolyDataMapper> mapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(polyData);

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    vtkNew<vtkNamedColors> colors;
    actor->GetProperty()->SetColor(colors->GetColor3d("Black").GetData());
    actor->GetProperty()->SetLineWidth(3.0);

    return actor;
}

vtkSmartPointer<vtkActor> getTranslateArrowShapeActor(
    const ARROWDIRECTION &direction)
{
    return getShapeActor(TRANSLATE_ARROW_HORIZONTAL_OFFSETS,
                         TRANSLATE_ARROW_VERTICAL_OFFSETS, direction);
}

vtkSmartPointer<vtkActor> getTranslateArrowOutlineActor(
    const ARROWDIRECTION &direction)
{
    return getOutlineActor(TRANSLATE_ARROW_HORIZONTAL_OFFSETS,
                           TRANSLATE_ARROW_VERTICAL_OFFSETS, direction);
}

vtkSmartPointer<vtkActor> getRotateArrowShapeActor(
    const ARROWDIRECTION &direction)
{
    std::vector<double> horizontalOffsets;
    std::vector<double> verticalOffsets;
    getRotateArrowPoints(horizontalOffsets, verticalOffsets);

    return getShapeActor(horizontalOffsets, verticalOffsets, direction);
}

vtkSmartPointer<vtkActor> getRotateArrowOutlineActor(
    const ARROWDIRECTION &direction)
{
    std::vector<double> horizontalOffsets;
    std::vector<double> verticalOffsets;
    getRotateArrowPoints(horizontalOffsets, verticalOffsets);

    return getOutlineActor(horizontalOffsets, verticalOffsets, direction);
}

vtkSmartPointer<vtkActor> getScaleArrowShapeActor()
{
    auto actor = getShapeActor(SCALE_ARROW_HORIZONTAL_OFFSETS,
                               SCALE_ARROW_VERTICAL_OFFSETS, ARROWDIRECTION::X);

    vtkSmartPointer<vtkTransform> trans = vtkSmartPointer<vtkTransform>::New();
    trans->RotateX(45.0);
    trans->RotateZ(45.0);
    trans->RotateY(45.0);
    actor->SetUserTransform(trans);

    return actor;
}

vtkSmartPointer<vtkActor> getScaleArrowOutlineActor()
{
    auto actor =
        getOutlineActor(SCALE_ARROW_HORIZONTAL_OFFSETS,
                        SCALE_ARROW_VERTICAL_OFFSETS, ARROWDIRECTION::X);

    vtkSmartPointer<vtkTransform> trans = vtkSmartPointer<vtkTransform>::New();
    trans->RotateX(45.0);
    trans->RotateZ(45.0);
    trans->RotateY(45.0);
    actor->SetUserTransform(trans);

    return actor;
}

modernTransformRepresentation::modernTransformRepresentation()
{
    std::array<vtkSmartPointer<vtkActor>, 3> translateShapeActors_ = {
        getTranslateArrowShapeActor(ARROWDIRECTION::X),
        getTranslateArrowShapeActor(ARROWDIRECTION::Y),
        getTranslateArrowShapeActor(ARROWDIRECTION::Z)};

    std::array<vtkSmartPointer<vtkActor>, 3> translateOutlineActors_ = {
        getTranslateArrowOutlineActor(ARROWDIRECTION::X),
        getTranslateArrowOutlineActor(ARROWDIRECTION::Y),
        getTranslateArrowOutlineActor(ARROWDIRECTION::Z)};

    for (auto i = 0; i < 3; i++)
    {
        translateActors_[i] = vtkSmartPointer<vtkAssembly>::New();
        translateActors_[i]->AddPart(translateShapeActors_[i]);
        translateActors_[i]->AddPart(translateOutlineActors_[i]);

        translateActors_[i]->SetOrigin(0.0, 0.0, 0.0);

        vtkSmartPointer<vtkTransform> trans =
            vtkSmartPointer<vtkTransform>::New();

        translateActors_[i]->SetUserTransform(trans);
    }

    std::array<vtkSmartPointer<vtkActor>, 3> rotateShapeActors_ = {
        getRotateArrowShapeActor(ARROWDIRECTION::Y),
        getRotateArrowShapeActor(ARROWDIRECTION::Z),
        getRotateArrowShapeActor(ARROWDIRECTION::X)};

    std::array<vtkSmartPointer<vtkActor>, 3> rotateOutlineActors_ = {
        getRotateArrowOutlineActor(ARROWDIRECTION::Y),
        getRotateArrowOutlineActor(ARROWDIRECTION::Z),
        getRotateArrowOutlineActor(ARROWDIRECTION::X)};

    for (auto i = 0; i < 3; i++)
    {
        rotateActors_[i] = vtkSmartPointer<vtkAssembly>::New();
        rotateActors_[i]->AddPart(rotateShapeActors_[i]);
        rotateActors_[i]->AddPart(rotateOutlineActors_[i]);

        rotateActors_[i]->SetOrigin(0.0, 0.0, 0.0);

        vtkSmartPointer<vtkTransform> trans =
            vtkSmartPointer<vtkTransform>::New();

        rotateActors_[i]->SetUserTransform(trans);
    }

    auto scaleShapeActor = getScaleArrowShapeActor();
    auto scaleOutlineActor = getScaleArrowOutlineActor();
    {
        scaleActor_ = vtkSmartPointer<vtkAssembly>::New();
        scaleActor_->AddPart(scaleShapeActor);
        scaleActor_->AddPart(scaleOutlineActor);

        scaleActor_->SetOrigin(0.0, 0.0, 0.0);

        vtkSmartPointer<vtkTransform> trans =
            vtkSmartPointer<vtkTransform>::New();

        scaleActor_->SetUserTransform(trans);
    }

    {
        assembleActor_ = vtkSmartPointer<vtkAssembly>::New();

        for (auto i = 0; i < 3; i++)
        {
            assembleActor_->AddPart(rotateActors_[i]);
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

    this->Highlight(0);
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
            actor->GetProperty()->SetOpacity(0.7);
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
                actor->GetProperty()->SetOpacity(1.0);
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
