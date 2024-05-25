#pragma once

#include <vtkMatrix4x4.h>
#include <vtkSmartPointer.h>
#include <vtkWidgetRepresentation.h>

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

class vtkAssembly;
class vtkCellPicker;
class transformRepresentation : public vtkWidgetRepresentation
{
   public:
    vtkTypeMacro(transformRepresentation, vtkWidgetRepresentation);

    virtual void GetTransform(vtkTransform *t) = 0;

   protected:
    transformRepresentation(){};
    virtual ~transformRepresentation() override{};

   private:
    transformRepresentation(const transformRepresentation &) = delete;
    void operator=(const transformRepresentation &) = delete;
};
