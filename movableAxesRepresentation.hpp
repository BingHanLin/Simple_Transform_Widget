#pragma once
#include <array>

#include <vtkMatrix4x4.h>
#include <vtkSmartPointer.h>
#include <vtkWidgetRepresentation.h>

class vtkAssembly;
class vtkCellPicker;
class movableAxesRepresentation : public vtkWidgetRepresentation
{
   public:
    /**
     * Instantiate the class.
     */
    static movableAxesRepresentation *New();

    ///@{
    /**
     * Standard methods for the class.
     */
    vtkTypeMacro(movableAxesRepresentation, vtkWidgetRepresentation);
    void PrintSelf(ostream &os, vtkIndent indent) override;
    ///@}

    ///@{
    /**
     * These are methods that satisfy vtkWidgetRepresentation's API.
     */
    void PlaceWidget(double bounds[6]) override;
    void BuildRepresentation() override;
    int ComputeInteractionState(int X, int Y, int modify = 0) override;
    void StartWidgetInteraction(double e[2]) override;
    void WidgetInteraction(double e[2]) override;
    double *GetBounds() VTK_SIZEHINT(6) override;

    /**
     * Overload the superclasses' GetMTime() because internal classes
     * are used to keep the state of the representation.
     */
    vtkMTimeType GetMTime() override;

    /**
     * Overridden to set the rendererer on the internal representations.
     */
    void SetRenderer(vtkRenderer *ren) override;

    void GetActors(vtkPropCollection *pc) override;

    enum INTERACTIONSTATE
    {
        outside = 0,
        onXRing,
        onYRing,
        onZRing
    };

   protected:
    movableAxesRepresentation();
    ~movableAxesRepresentation() override;

   private:
    const std::array<std::array<double, 3>, 3> axisNormalColor_ = {
        {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};

    std::array<vtkSmartPointer<vtkAssembly>, 3> axisRingActors_;
    std::array<vtkSmartPointer<vtkMatrix4x4>, 3> axisRingInitMatrix_;

    vtkSmartPointer<vtkCellPicker> picker_;
    vtkProp3D *currActor_ = nullptr;

    std::array<double, 3> startEventPosition_;
    std::array<double, 3> lastEventPosition_;

    movableAxesRepresentation(const movableAxesRepresentation &) = delete;
    void operator=(const movableAxesRepresentation &) = delete;
};
