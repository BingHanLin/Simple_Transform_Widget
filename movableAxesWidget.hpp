#pragma once

#include "vtkAbstractWidget.h"

class vtkCallbackCommand;

class movableAxesRepresentation;

class movableAxesWidget : public vtkAbstractWidget {
public:
  /**
   * Instantiate the object.
   */
  static movableAxesWidget *New();

  /**
   * Standard vtkObject methods
   */
  vtkTypeMacro(movableAxesWidget, movableAxesWidget);
  void PrintSelf(ostream &os, vtkIndent indent) override;

  /**
   * Override superclasses' SetEnabled() method because the line
   * widget must enable its internal handle widgets.
   */
  void SetEnabled(int enabling) override;

  /**
   * Specify an instance of vtkWidgetRepresentation used to represent this
   * widget in the scene. Note that the representation is a subclass of vtkProp
   * so it can be added to the renderer independent of the widget.
   */
  void SetRepresentation(movableAxesRepresentation *r) {
    this->Superclass::SetWidgetRepresentation(
        reinterpret_cast<vtkWidgetRepresentation *>(r));
  }

  /**
   * Create the default widget representation if one is not set.
   */
  void CreateDefaultRepresentation() override;

  /**
   * Methods to change the whether the widget responds to interaction.
   * Overridden to pass the state to component widgets.
   */
  void SetProcessEvents(vtkTypeBool enabled) override;

protected:
  movableAxesWidget();
  ~movableAxesWidget() override;

private:
  enum WIDGETSTATE { Start = 0, Active };

  int state_;
  vtkCallbackCommand *keyEventCallbackCommand_;

  movableAxesWidget(const movableAxesWidget &) = delete;
  void operator=(const movableAxesWidget &) = delete;
};
