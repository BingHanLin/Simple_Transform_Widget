#include <vtkActor.h>
#include <vtkAxesActor.h>
#include <vtkCamera.h>
#include <vtkCommand.h>
#include <vtkConeSource.h>
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
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkTransform.h>

#include "modernTransformRepresentation.hpp"
#include "transformWidget.hpp"

// This does the actual work.
// Callback for the interaction
class transformCallback : public vtkCommand
{
   public:
    static transformCallback *New()
    {
        return new transformCallback;
    }

    virtual void Execute(vtkObject *caller, unsigned long, void *)
    {
        vtkSmartPointer<transformWidget> widget =
            dynamic_cast<transformWidget *>(caller);

        if (widget != nullptr && actor_ != nullptr)
        {
            if (auto rep = dynamic_cast<transformRepresentation *>(
                    widget->GetRepresentation()))
            {
                vtkNew<vtkTransform> trans;
                rep->GetTransform(trans);

                vtkNew<vtkMatrix4x4> newMatrix;
                newMatrix->DeepCopy(trans->GetMatrix());

                actor_->SetUserMatrix(newMatrix);
            }
        }
    }

    void setActor(vtkSmartPointer<vtkActor> actor)
    {
        actor_ = actor;
    };

    transformCallback() : actor_(nullptr)
    {
    }

   private:
    vtkSmartPointer<vtkActor> actor_;
};

int main(int, char *[])
{
    vtkNew<vtkCamera> camera;
    camera->SetParallelProjection(true);

    vtkNew<vtkNamedColors> colors;

    // A mainRenderer and render window
    vtkNew<vtkRenderer> mainRenderer;
    mainRenderer->SetBackground(colors->GetColor3d("BurlyWood").GetData());
    mainRenderer->SetActiveCamera(camera);
    mainRenderer->SetLayer(0);

    vtkNew<vtkRenderer> frontRenderer;
    frontRenderer->SetActiveCamera(camera);
    frontRenderer->SetLayer(1);

    vtkNew<vtkRenderWindow> renderWindow;
    renderWindow->SetNumberOfLayers(2);
    renderWindow->AddRenderer(mainRenderer);
    renderWindow->AddRenderer(frontRenderer);

    // An interactor
    vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
    renderWindowInteractor->SetRenderWindow(renderWindow);

    vtkNew<vtkSphereSource> sphereSource;
    sphereSource->Update();

    // Create a mapper and actor
    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputConnection(sphereSource->GetOutputPort());
    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(colors->GetColor3d("MistyRose").GetData());

    // Cube Actor
    vtkNew<vtkConeSource> cone;
    cone->SetCenter(1.0, 1.0, 1.0);
    cone->SetHeight(5);
    cone->SetRadius(2.5);
    cone->Update();

    vtkNew<vtkPolyDataMapper> coneMapper;
    coneMapper->SetInputData(cone->GetOutput());

    vtkNew<vtkActor> coneActor;
    coneActor->SetMapper(coneMapper);
    coneActor->GetProperty()->SetColor(colors->GetColor3d("Banana").GetData());
    coneActor->SetVisibility(true);
    mainRenderer->AddActor(coneActor);

    // Add transformWidget
    vtkNew<modernTransformRepresentation> rep;
    vtkNew<transformWidget> myWidget;
    myWidget->SetInteractor(renderWindowInteractor);
    myWidget->SetRepresentation(rep);

    double bounds[6] = {-2.0, 4.0, -2.0, 4.0, -2.0, 4.0};
    myWidget->GetRepresentation()->PlaceWidget(bounds);

    vtkNew<transformCallback> callback;
    callback->setActor(coneActor);
    myWidget->AddObserver(vtkCommand::InteractionEvent, callback);
    myWidget->SetEnabled(1);

    // Add vtkCubeAxesActor
    auto cubeAxesActor = vtkSmartPointer<vtkCubeAxesActor>::New();
    cubeAxesActor->SetBounds(-5.0, 5.0, -5.0, 5.0, -5.0, 5.0);
    cubeAxesActor->SetCamera(mainRenderer->GetActiveCamera());
    mainRenderer->AddActor(cubeAxesActor);

    // Axes Widget
    auto vtkAxes = vtkSmartPointer<vtkAxesActor>::New();
    auto axesWidget = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
    axesWidget->SetOrientationMarker(vtkAxes);
    axesWidget->SetDefaultRenderer(mainRenderer);
    axesWidget->SetInteractor(renderWindowInteractor);
    axesWidget->SetViewport(0.0, 0.0, 0.20, 0.20);
    axesWidget->SetEnabled(1);
    axesWidget->InteractiveOff();

    // Reset camera
    double allBounds[6];
    mainRenderer->ComputeVisiblePropBounds(allBounds);
    if (vtkMath::AreBoundsInitialized(allBounds))
    {
        mainRenderer->ResetCamera(allBounds);
    }

    // Render
    renderWindow->Render();

    renderWindowInteractor->Initialize();
    renderWindow->Render();

    // Begin mouse interaction
    renderWindowInteractor->Start();

    return EXIT_SUCCESS;
}
