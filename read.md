## vtkAbstractWidget
The vtkAbstractWidget defines an API and implements methods common to all widgets using the interaction/representation design. In this design, the term interaction means that part of the widget that performs event handling, while the representation corresponds to a vtkProp (or the subclass vtkWidgetRepresentation) used to represent the widget. vtkAbstractWidget also implements some methods common to all subclasses.

Note that vtkAbstractWidget provides access to the vtkWidgetEventTranslator. This class is responsible for translating VTK events (defined in vtkCommand.h) into widget events (defined in vtkWidgetEvent.h). This class can be manipulated so that different VTK events can be mapped into widget events, thereby allowing the modification of event bindings. Each subclass of vtkAbstractWidget defines the events to which it responds.

// must override
CreateDefaultRepresentation()
SetProcessEvents()

https://vtk.org/doc/nightly/html/classvtkAbstractWidget.html
