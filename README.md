## vtkAbstractWidget
The vtkAbstractWidget defines an API and implements methods common to all widgets using the interaction/representation design. In this design, the term interaction means that part of the widget that performs event handling, while the representation corresponds to a vtkProp (or the subclass vtkWidgetRepresentation) used to represent the widget. vtkAbstractWidget also implements some methods common to all subclasses.

Note that vtkAbstractWidget provides access to the vtkWidgetEventTranslator. This class is responsible for translating VTK events (defined in vtkCommand.h) into widget events (defined in vtkWidgetEvent.h). This class can be manipulated so that different VTK events can be mapped into widget events, thereby allowing the modification of event bindings. Each subclass of vtkAbstractWidget defines the events to which it responds.

// must override
CreateDefaultRepresentation()
SetProcessEvents()

https://vtk.org/doc/nightly/html/classvtkAbstractWidget.html

vtkCameraOrientationWidget 
https://vtk.org/doc/nightly/html/classvtkCameraOrientationWidget.html

vtkFollower
https://vtk.org/doc/nightly/html/classvtkFollower.html

vtkInteractorStyleTrackballActor
https://vtk.org/doc/nightly/html/classvtkInteractorStyleTrackballActor.html#details

vtkInteractorStyleTrackballActor example
https://kitware.github.io/vtk-examples/site/Python/Interaction/InteractorStyleTrackballActor/

VTK可移动三维坐标轴 vtkMovableAxesWidget
https://blog.csdn.net/weixin_38416696/article/details/124987253

[VTK] Draggable coordinate axis MovableAxesWidget
https://www.zditect.com/code/vtk-draggable-coordinate-axis-movableaxeswidget.html


VTKWidgets(Explain the concept and event translation table)
https://vtk.org/Wiki/VTKWidgets

VTK/Tutorials/New Pipeline
https://vtk.org/Wiki/VTK/Tutorials/New_Pipeline

Parallel Processing with VTK's SMP Framework
https://vtk.org/doc/nightly/html/md__builds_gitlab_kitware_sciviz_ci_Documentation_Doxygen_SMPTools.html