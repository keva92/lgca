# LGCA
LGCA is a C++ implementation of a Lattice Gas Cellular Automaton (LGCA) in 2D. Features comprise

* HPP and FHP-I, FHP-II, and FHP-III lattice gas models
* Several apps with dedicated GUIs (including Kármán vortex street, pipe flow, and diffusion)
* Shared memory parallelization
* On-line data visualization
* File I/O using vtkImageData (.vti) (supported by ParaView) as well as PNG images


## Build from Source

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
Build the project and run one of the apps:
```
make
cd ../bin/
./lgca-pipe
```

## Theory

Hopefully, some basic theory on lattice gas cellular automatons will follow soon. Until then, feel free to have a look at Dieter Wolf-Gladrow's "Lattice-Gas Cellular Automata and Lattice Boltzmann Models".


## Background

The initial code has been developed at Hamburg University of Technology (TUHH) in 2015 in the framework of the postgrad course "Application of Innovative CFD Methods in Research and Development" by the then students Niklas Kühl and Kerstin Vater. The project has been initiated and supervised by Christian F. Janßen, postdoc at the Institute of Fluid Dynamics and Ship Theory at that time. 
