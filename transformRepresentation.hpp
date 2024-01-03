#pragma once
#include <array>

#include <vtkMatrix4x4.h>
#include <vtkSmartPointer.h>
#include <vtkWidgetRepresentation.h>

class vtkAssembly;
class vtkCellPicker;
class transformRepresentation : public vtkWidgetRepresentation
{
   public:
    enum INTERACTIONSTATE
    {
        outside = 0,
        onXRing,
        onYRing,
        onZRing,
        onXArrow,
        onYArrow,
        onZArrow,
        onScale
    };

    static transformRepresentation *New();
    vtkTypeMacro(transformRepresentation, vtkWidgetRepresentation);

    void PlaceWidget(double bounds[6]) override;
    void StartWidgetInteraction(double e[2]) override;
    void WidgetInteraction(double e[2]) override;
    int ComputeInteractionState(int X, int Y, int modify = 0) override;
    void Highlight(int highlight) override;

    void BuildRepresentation() override;

    void GetActors(vtkPropCollection *pc) override;

    void GetTransform(vtkTransform *t);

   protected:
    transformRepresentation();
    ~transformRepresentation() override;

   private:
    const std::array<std::array<double, 3>, 3> axisNormalColor_ = {
        {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};

    std::array<vtkSmartPointer<vtkAssembly>, 3> rotateActors_;
    std::array<vtkSmartPointer<vtkAssembly>, 3> translateActors_;
    vtkSmartPointer<vtkAssembly> scaleActor_;

    vtkSmartPointer<vtkAssembly> operationActor_;

    vtkSmartPointer<vtkCellPicker> picker_;

    std::array<double, 3> prevEventPosition_;
    std::array<double, 4> prevEventWorldPosition_;
    std::array<double, 4> currEventWorldPosition_;

    transformRepresentation(const transformRepresentation &) = delete;
    void operator=(const transformRepresentation &) = delete;
};
