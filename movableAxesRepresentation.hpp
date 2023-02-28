#pragma once

#include "vtkDeprecation.h"              // For VTK_DEPRECATED_IN_9_0_0
#include "vtkInteractionWidgetsModule.h" // For export macro
#include "vtkWidgetRepresentation.h"

class vtkActor;
class vtkConeSource;
class vtkPolyDataMapper;
class vtkLineSource;
class vtkProperty;
class vtkPolyData;
class vtkPolyDataAlgorithm;
class vtkPointHandleRepresentation3D;
class vtkBox;
class vtkFollower;
class vtkVectorText;
class vtkPolyDataMapper;
class vtkCellPicker;

class movableAxesRepresentation : public vtkWidgetRepresentation {
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

protected:
  movableAxesRepresentation();
  ~movableAxesRepresentation() override;

private:
  movableAxesRepresentation(const movableAxesRepresentation &) = delete;
  void operator=(const movableAxesRepresentation &) = delete;
};
