# Viscous River

* Note: this project makes use of a lot material from the book ___Finite Element Methods for Engineers___ by ___Roger T. Fenner___.
* Fenner provides a `Fortran 77` solver for the problem we're looking at.
* My solution to the case study involves refactoring the software into __C++11__, utilising __OOD__ principles and the __Eigen Template Library__.

### TODO
* Add basic notes r.e. the underlying theory behind `FEM`.
* Data plotting. A `Python` script using `Matplotlib` or (https://www.paraview.org/)[ParavView] is a good choice for such applications; I'm not too concerned about this ATM.
* __Issue__: If you don't press 1 for `verbose mode`, no data is written to the files so press it.
* Use better Software Engineering principles and tools. I.e.
   - `CMake`. A `Makefile` is provided for the build process
   - Code organisation. I.e. classes into header files etc.

## Introduction

This software consists of a `C++ solver` which solves the `Navier-Stokes Equation(s)` for a `viscous fluid` flowing in a `channel/river`.
The domain is `discretised` via the `(Direct) Finite Element Method (FEM)`. The resulting `System of Algebraic Equations` is handled and solved via the `Eigen C++ template library`.
This software calculates the __Velocity__ and using this, the __Flow Rate__ through the channel.

## Model

* The model consists of a `channel` or rather, a `river`, flowing along the `z-axis (direction)`.
* Note: all physical quantities use `S. I. units`.
* The viscous fluid flows with `uniform velocity`, and has `no external forces` acting on it.
* With this is mind, the `Navier-Stokes Equation(s)` are reduced to:

	`grad(grad . v) = (1/mu)dp/dz`

* `grad`: the (vector) gradient operator. Note the left-hand side is collectively called the `Laplacian of the velocity field`.

* `mu`: `Dynamic (Shear) Viscosity`: Dynamic viscosity is a fluid's resistance to flow and shear when an external force is applied. Taken to be constant.

* `dp/dz`: This is the `Pressure Gradient`: the rate at which pressure changes over distance. Taken to be constant.

* __Dirichlet Boundary Conditions__ (BCs): At the bottom and sides of the channel/river: `v = 0 m/s`. At the top (surface): fixed `v = vz`.
* The shape of the `domain` is a rectangle of `height H` and `width W`. Values for these: `H = 1m` and `W = 2m`.
* `Dynamic Viscosity, mu = 1 N.s/m^2`.
* `Pressure Gradient, dp/dz = -6 N/m^3`.
* The ___Flow Rate___ through the channel/river is given by 

	`R = (1/2)*W*H*Vz*F_D - (1/2mu)*W*H^3*F_P`

* `F_D`: ___Drag Factor___; associated with the drag flow induced by the boundary velocity `Vz`.
* `F_P`: ___Pressure Flow Shape___; due to the pressure gradient `p'` or `dp/dz`. 
* The values of `F_D` & `F_P` provide an important ratio in modelling channel flow problems: the `aspect ratio`, `H/W`.
* Both the `F_D` & `F_P` factors are __Hyperbolic Functions__. 
* `H/W = 0.5` for this problem.
* Properties:

	- __Shallow and wider channels__: larger wetted perimeter leading to increased friction and slower flow velocities.
	- __Deeper and narrower channels__: smaller wetted perimeter resulting in less friction and faster flow velocities.
* The number of nodes (___vertices___) in the `x-direction` is given by `NodeX`. Those in the `y-direction`: `NodeY`.
* The __mesh__ of the problem domain consists of `6 x 6 nodal points`.
* The number of __triangular elements__ that make up the mesh is denoted by `elementNumber`. 50 will be used in total.

## Finite Element Method (FEM)

* This is a brief overview of what the method entails.

### Introduction

* TODO.

### The Variational Method: FEM

* TODO.

### Mesh Generation

* TODO.

## Requirements

* Developed and tested on `Linux (Ubuntu 20.04)`.
* Compiler: developed & tested with `g++ 13.1.0`.
* The `Eigen C++ Template Library` for the `Linear Algebra`; version used: `3.4.0`.
* `CMake` for building the software etc.
* Knowledge of `Applied Mathematics`; e.g. `FEM`, `Numerical Linear Algebra & ODEs/PDEs`.

## Eigen Library

* Eigen is a C++ library of template headers used extensively for `Numerical Linear Algebra`.

* Other uses:

	- `Matrix operations`,
	- `Vector operations`,
	- `Geometric transformations`.

### Eigen Installation

* This is a quick tutorial on how to install and run a program which uses `Eigen`. I've included this because I feel that a good tutorial for begginers is missing, and the ones that are available are lacking/missing key steps.

* [Download the latest version](http://eigen.tuxfamily.org/index.php?title=Main_Page#Download).
* Extract it and move it to `/home/user`.

`~/Downloads$ tar -xf eigen-3.4.0.tar.gz`
`~/Downloads$ mv eigen-3.4.0 /home/user`

* You'll need to make a `build directory` for `CMake`. 

`~$ cd eigen-3.4.0`
`~/eigen-3.4.0$ mkdir build && cd build`

* Run `CMake`

`~/eigen-3.4.0/build$ cmake ..`

* Now install `Eigen`

`~/eigen-3.4.0/build$ sudo make install`

* symlink or copy the Eigen folder into `/usr/local/include/`

`~/eigen-3.4.0/build$ sudo cp -r /usr/local/include/eigen3/Eigen /usr/local/include`

* And you're set! See the test program in `eigen_test` to ensure everything is working.

* Inside the `eigen_test` folder, execute the following at the CLI:

`$ g++ main.cpp -o main`


## Getting & Running the Software

* `$ git clone https://github.com/MRLintern/viscousRiver.git`
* `$ cd viscousRiver`
* `$ make main`
* `$ ./main`
* `Press 1 for verbose mode`.
* To delete the binary and files:
* `$ make clean`

* __TODO__; not using `CMake` yet.
* `$ mkdir build && cd build`
* `$ cmake ..`
