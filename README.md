## vtkAbstractWidget
The vtkAbstractWidget defines an API and implements methods common to all widgets using the interaction/representation design. In this design, the term interaction means that part of the widget that performs event handling, while the representation corresponds to a vtkProp (or the subclass vtkWidgetRepresentation) used to represent the widget. vtkAbstractWidget also implements some methods common to all subclasses.

Note that vtkAbstractWidget provides access to the vtkWidgetEventTranslator. This class is responsible for translating VTK events (defined in vtkCommand.h) into widget events (defined in vtkWidgetEvent.h). This class can be manipulated so that different VTK events can be mapped into widget events, thereby allowing the modification of event bindings. Each subclass of vtkAbstractWidget defines the events to which it responds.

### must override
* CreateDefaultRepresentation()
* SetProcessEvents()


## vtkWidgetRepresentation 

### must override
* SetRenderer() - the renderer in which the representations draws itself.
Typically the renderer is set by the associated widget.
Use the widget's SetCurrentRenderer() method in most cases;
otherwise there is a risk of inconsistent behavior as events
and drawing may be performed in different viewports.
* BuildRepresentation() - update the geometry of the widget based on its
current state.

### suggest override
* PlaceWidget() - given a bounding box (xmin,xmax,ymin,ymax,zmin,zmax), place
the widget inside of it. The current orientation of the widget
is preserved, only scaling and translation is performed.
StartWidgetInteraction() - generally corresponds to a initial event (e.g.,
mouse down) that starts the interaction process
with the widget.
* WidgetInteraction() - invoked when an event causes the widget to change
appearance.
* EndWidgetInteraction() - generally corresponds to a final event (e.g., mouse up)
and completes the interaction sequence.
* ComputeInteractionState() - given (X,Y) display coordinates in a renderer, with a
possible flag that modifies the computation,
what is the state of the widget?
* GetInteractionState() - return the current state of the widget. Note that the
value of "0" typically refers to "outside". The
interaction state is strictly a function of the
representation, and the widget/represent must agree
on what they mean.
* Highlight() - turn on or off any highlights associated with the widget.
Highlights are generally turned on when the widget is selected.


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

Virtual surgical planning using a DCIA free flap for mandibular reconstruction
https://www.youtube.com/watch?v=7RMermYOOYY&t=96s

VTK——可拖动的坐标轴MovableAxesWidget、商用
https://blog.51cto.com/u_15926338/5979963

VTKWidgets(Explain the concept and event translation table)
https://vtk.org/Wiki/VTKWidgets

VTK/Tutorials/New Pipeline
https://vtk.org/Wiki/VTK/Tutorials/New_Pipeline

Parallel Processing with VTK's SMP Framework
https://vtk.org/doc/nightly/html/md__builds_gitlab_kitware_sciviz_ci_Documentation_Doxygen_SMPTools.html