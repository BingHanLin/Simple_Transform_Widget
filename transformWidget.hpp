#pragma once

#include "vtkAbstractWidget.h"

class vtkCallbackCommand;

class movableAxesRepresentation;

class transformWidget : public vtkAbstractWidget
{
   public:
    static transformWidget *New();

    vtkTypeMacro(transformWidget, transformWidget);
    void PrintSelf(ostream &os, vtkIndent indent) override;

    void SetEnabled(int enabling) override;

    void SetRepresentation(movableAxesRepresentation *r)
    {
        this->Superclass::SetWidgetRepresentation(
            reinterpret_cast<vtkWidgetRepresentation *>(r));
    }

    void CreateDefaultRepresentation() override;

   protected:
    transformWidget();
    ~transformWidget() override;

    static void SelectAction(vtkAbstractWidget *w);
    static void EndSelectAction(vtkAbstractWidget *w);
    static void MoveAction(vtkAbstractWidget *w);
    static void ProcessKeyEvents(vtkObject *, unsigned long, void *, void *);

   private:
    enum WIDGETSTATE
    {
        start = 0,
        active
    };

    int state_;
    vtkCallbackCommand *keyEventCallbackCommand_;

    transformWidget(const transformWidget &) = delete;
    void operator=(const transformWidget &) = delete;
};
