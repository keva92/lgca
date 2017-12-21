# LGCA
LGCA is a C++ implementation of a Lattice Gas Cellular Automaton (LGCA) in 2D. Features comprise
* HPP and FHP-I lattice gas models
* Several apps with dedicated GUIs (including Kármán vortex street, pipe flow, and diffusion)
* On-line data visualization
* File I/O using vtkImageData (.vti) (supported by ParaView)

## Installation

Please make sure you have installed CMake version 3.8 or higher on your local machine. You will also need Git and some basic development tools (gcc, g++, make...).

Switch to the directory where you would like to install the code. Clone the Git repository of the project:
```
git clone https://github.com/keva92/lgca.git
```
Create a build directory and configure the build process:
```
cd lgca/
mkdir build
cd build/
ccmake ..
```
Build the project and and run one of the apps:
```
make
cd ../bin/
./lgca-pipe
```
