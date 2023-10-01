#include <vtkActor.h>
#include <vtkAxesActor.h>
#include <vtkCamera.h>
#include <vtkCommand.h>
#include <vtkCubeAxesActor.h>
#include <vtkCubeSource.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSphereSource.h>

#include "movableAxesRepresentation.hpp"
#include "movableAxesWidget.hpp"

// This does the actual work.
// Callback for the interaction
class vtkLineCallback : public vtkCommand
{
   public:
    static vtkLineCallback *New()
    {
        return new vtkLineCallback;
    }

    virtual void Execute(vtkObject *caller, unsigned long, void *)
    {
        // movableAxesWidget *lineWidget =
        //     reinterpret_cast<movableAxesWidget *>(caller);

        // // Get the actual box coordinates of the line
        // vtkNew<vtkPolyData> polydata;
        // static_cast<movableAxesRepresentation
        // *>(lineWidget->GetRepresentation())
        //     ->GetPolyData(polydata);

        // // Display one of the points, just so we know it's working
        // double p[3];
        // polydata->GetPoint(0, p);
        // std::cout << "P: " << p[0] << " " << p[1] << " " << p[2] <<
        // std::endl;
    }
    vtkLineCallback()
    {
    }
};

int main(int, char *[])
{
    vtkNew<vtkNamedColors> colors;

    vtkNew<vtkSphereSource> sphereSource;
    sphereSource->Update();

    // Create a mapper and actor
    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputConnection(sphereSource->GetOutputPort());
    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(colors->GetColor3d("MistyRose").GetData());

    // A renderer and render window
    vtkNew<vtkRenderer> renderer;
    vtkNew<vtkRenderWindow> renderWindow;
    renderWindow->AddRenderer(renderer);
    renderWindow->SetWindowName("LineWidget2");

    // renderer->AddActor(actor);
    renderer->SetBackground(colors->GetColor3d("MidnightBlue").GetData());

    // An interactor
    vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
    renderWindowInteractor->SetRenderWindow(renderWindow);

    vtkNew<movableAxesWidget> myWidget;
    myWidget->SetInteractor(renderWindowInteractor);
    myWidget->CreateDefaultRepresentation();
    double bounds[6] = {-1.0, 1.0, -1.0, 1.0, -1.0, 1.0};
    myWidget->GetRepresentation()->PlaceWidget(bounds);

    // You could do this if you want to set properties at this point:
    // vtkNew<movableAxesRepresentation> lineRepresentation;
    // lineWidget->SetRepresentation(lineRepresentation);

    vtkNew<vtkLineCallback> lineCallback;

    // myWidget->AddObserver(vtkCommand::InteractionEvent, lineCallback);

    // vtkCubeAxesActor Actor
    auto cubeAxesActor = vtkSmartPointer<vtkCubeAxesActor>::New();
    cubeAxesActor->SetBounds(-5.0, 5.0, -5.0, 5.0, -5.0, 5.0);
    cubeAxesActor->SetCamera(renderer->GetActiveCamera());
    renderer->AddActor(cubeAxesActor);

    // Cube Actor
    vtkNew<vtkCubeSource> cube;
    cube->SetBounds(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
    cube->Update();

    vtkNew<vtkPolyDataMapper> cubeMapper;
    cubeMapper->SetInputData(cube->GetOutput());

    vtkNew<vtkActor> cubeActor;
    cubeActor->SetMapper(cubeMapper);
    cubeActor->GetProperty()->SetColor(colors->GetColor3d("Banana").GetData());
    cubeActor->SetVisibility(true);
    renderer->AddActor(cubeActor);

    // Axes Widget
    auto vtkAxes = vtkSmartPointer<vtkAxesActor>::New();
    auto axesWidget = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
    axesWidget->SetOrientationMarker(vtkAxes);
    axesWidget->SetDefaultRenderer(renderer);
    axesWidget->SetInteractor(renderWindowInteractor);
    axesWidget->SetViewport(0.0, 0.0, 0.20, 0.20);
    axesWidget->SetEnabled(1);
    axesWidget->InteractiveOff();

    // Reset camera
    double allBounds[6];
    renderer->ComputeVisiblePropBounds(allBounds);
    if (vtkMath::AreBoundsInitialized(allBounds))
    {
        renderer->ResetCamera(allBounds);
    }

    // Render
    renderWindow->Render();

    renderWindowInteractor->Initialize();
    renderWindow->Render();
    myWidget->On();

    // Begin mouse interaction
    renderWindowInteractor->Start();

    return EXIT_SUCCESS;
}
