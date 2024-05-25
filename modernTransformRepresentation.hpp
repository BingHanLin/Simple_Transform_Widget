#pragma once
#include <array>

#include <vtkActor.h>
#include <vtkMatrix4x4.h>
#include <vtkSmartPointer.h>
#include <vtkWidgetRepresentation.h>

#include "transformRepresentation.hpp"

class vtkAssembly;
class vtkCellPicker;
class modernTransformRepresentation : public transformRepresentation
{
   public:
    static modernTransformRepresentation *New();
    vtkTypeMacro(modernTransformRepresentation, transformRepresentation);

    void PlaceWidget(double bounds[6]) override;
    void StartWidgetInteraction(double e[2]) override;
    void WidgetInteraction(double e[2]) override;
    int ComputeInteractionState(int X, int Y, int modify = 0) override;
    void Highlight(int highlight) override;

    void BuildRepresentation() override;

    void GetActors(vtkPropCollection *pc) override;

    void GetTransform(vtkTransform *t);

   protected:
    modernTransformRepresentation();

   private:
    const std::array<std::array<double, 3>, 3> axisNormalColor_ = {
        {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};

    std::array<vtkSmartPointer<vtkActor>, 3> rotateShapeActors_;
    std::array<vtkSmartPointer<vtkActor>, 3> rotateOutlineActors_;
    std::array<vtkSmartPointer<vtkActor>, 3> translateShapeActors_;
    std::array<vtkSmartPointer<vtkActor>, 3> translateOutlineActors_;
    vtkSmartPointer<vtkActor> scaleShapeActor_;
    vtkSmartPointer<vtkActor> scaleOutlineActor_;

    std::array<vtkSmartPointer<vtkAssembly>, 3> rotateActors_;
    std::array<vtkSmartPointer<vtkAssembly>, 3> translateActors_;
    vtkSmartPointer<vtkAssembly> scaleActor_;
    vtkSmartPointer<vtkAssembly> assembleActor_;

    vtkSmartPointer<vtkCellPicker> picker_;

    std::array<double, 3> prevEventPosition_;
    std::array<double, 4> prevEventWorldPosition_;
    std::array<double, 4> currEventWorldPosition_;

    std::array<double, 3> placeCenter_;
    std::array<double, 3> placeScale_;

    modernTransformRepresentation(const modernTransformRepresentation &) =
        delete;
    void operator=(const modernTransformRepresentation &) = delete;
};
