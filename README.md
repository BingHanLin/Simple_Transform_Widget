# Simple VTK Transform Widget

## Overview

This project implements a VTK (Visualization Toolkit) transform widget. This widget allows users to interactively manipulate 3D objects by applying translation and rotation.

![Transform widget](./asset/animation.gif)

## Features
* Interactive 3D transformation widget.
* Supports translation, rotation ans scale of 3D objects.
* Follow the widget/represent implementation to decouple Event Processing from Widget Geometry.
* MIT License for easy integration into your own projects.

## Usage
1. Clone the repository to your local machine.
    ```
    git clone https://github.com/yourusername/vtk-transformation-widget.git
    ```
2. Configure project through cmake command. You should replace `<YOUR_VTK_LIB_PATH>` with the actual path to your VTK library.
    ```
    cmake   -S . -B build \
            -D CMAKE_BUILD_TYPE=Release \
            -D CMAKE_INSTALL_PREFIX="./install" \
            -D VTK_DIR=<YOUR_VTK_LIB_PATH> \
            -G "Visual Studio 16 2019" -A x64
    ```
3. Build and install the project through cmake command.
    ```
    cmake --build build --config Release

    cmake --install build --config Release
    ```
4. Find and execute the main.exe in the `install` folder. To run it successfully,you may need to incdue VTK binary path in the enviroment variable `PATH `.

## License
This project is licensed under the MIT License - see the [LICENSE](https://opensource.org/license/mit/) file for details.

## Contributions
Contributions to this project are welcome. Please feel free to submit issues or pull requests to improve the functionality or fix any bugs.