# Viscous River

Note: this is a work-in-progress; bare with me!

## Introduction

This software consists of a `C++ solver` which solves the `Navier-Stokes Equation(s)` for a `viscous fluid` flowing in a `channel/river`.
The domain is `discretised` via the `Finite Element Method (FEM)`. The resulting `System of Algebraic Equations` is handled and solved via the `Eigen C++ template library`.

## Model

* The model consists of a `channel` or rather, a `river`, flowing along the `z-axis (direction)`.
* Note: all physical quantities use `S. I. units`.
* The viscous fluid flows with `uniform velocity`, and has `no external forces` acting on it.
* With this is mind, the `Navier-Stokes Equation(s)` are reduced to:

	`grad(grad . v) = (1/mu)dp/dz`

* `grad`: the (vector) gradient operator. Note the left-hand side is collectively called the `Laplacian of the velocity field`.

* `mu`: `Dynamic (Shear) Viscosity`: Dynamic viscosity is a fluid's resistance to flow and shear when an external force is applied. Taken to be constant.

* `dp/dz`: This is the `Pressure Gradient`: the rate at which pressure changes over distance. Taken to be constant.

* `Boundary Conditions (BCs)`: At the bottom and sides of the channel/river: `v = 0 m/s`. At the top (surface): `v = vz`.

* The shape of the `domain` is a rectangle of `height H` and `width W`. Values for these: `H = 2m` and `W = 6m`.
* `Dynamic Viscosity, mu = 1 N.s/m^2`.
* `Pressure Gradient, dp/dz = -8 N/m^3`.
* The `mesh/discretisised domain` consists of `100 x 100` `nodes/points`.

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

`$ g++ main.cpp -o main`


## Getting & Running the Software

`$ git clone https://github.com/MRLintern/viscousRiver.git`

* Once the software has been developed, further instructions on building and running it will be available. 
