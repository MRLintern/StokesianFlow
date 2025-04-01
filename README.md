# Viscous River

* Note: this project makes use of a lot material from the book ___Finite Element Methods for Engineers___ by ___Roger T. Fenner___.
* Fenner provides a `Fortran 77` solver for the problem we're looking at.
* My solution to the case study involves refactoring the software into __C++14__, utilising __OOD__ principles, __Modern C++ Memory Management__ and the __Eigen Template Library__.

### TODO
* Data plotting. A `Python` script using `Matplotlib` or [ParaView](https://www.paraview.org/) is a good choice for such applications; I'm not too concerned about this ATM.
* Use better Software Engineering principles and tools. I.e.
   - Code organisation. E.g. classes into header files etc.
   - At present, I'm currently developing: `Vector.hpp`, `Vector.cpp`, `Solver.hpp`, `Solver.cpp` and `main.cpp`.

### TODO: Further Work
* __Software Considerations__: Increase the number of nodes significantly to increase the granularity of the mesh. Incorporating the __OpenMP__ or __MPI__ API might be helpful with this.
* __Software Considerations__: __Optimisation Flags__; experiment on optimising the compiler for different levels of optimisation.
* __Physics__: Increase Dynamic viscosity, `mu`.
* __Physics__: Vary `H/W` ratio.

## Introduction

This software, developed in Modern C++, models a Viscous Fluid flowing 
This software models 
The domain is `discretised` via the `Finite Element Method (FEM)`. The resulting `System of Algebraic Equations` is handled and solved via the `Eigen C++ template library`.
This software calculates the __Velocity Profile__ and using this, the __Volumetric Flow Rate__ through the channel.

## Model

* The model consists of a `channel` or rather, a `river`, flowing along the `z-axis (direction)`.
* Note: all physical quantities use `S. I. units`.
* The viscous fluid flows with `uniform velocity`, and has `no external forces` acting on it. The only force acting on the flow is the __Pressure Gradient__.
* If the rectangular channel is infinitely wide compared to its depth, `H << W`, the governing ODE is given by:

	`d^2 u/dy^2 = (1/mu)dp/dz`

* `d^2 u/d^2`: Second order derivative of the __Velcoity Field__.

* `mu`: `Dynamic (Shear) Viscosity`: Dynamic viscosity is a fluid's resistance to flow and shear when an external force is applied. Taken to be constant.

* `dp/dz`: This is the `Pressure Gradient`: the rate at which pressure changes over distance. Taken to be constant. In the software, this represented by `Pz`.

* Integrating the ODE with Boundary Conditions (BC, see below) to find the __Velocity Profile__ gives:

  	`u = yVz/H + Pz/2mu(y^2 - yH)`

* __Dirichlet Boundary Conditions__ (BCs): At the bottom and sides of the channel/river: `v = 0 m/s`. At the top (surface): fixed `v = Vz`.
* The shape of the `domain` is a rectangle of `height H` and `width W`. Values for these: `H = 1m` and `W = 2m`.
* `Dynamic Viscosity, mu = 1 N.s/m^2`.
* `Pressure Gradient, dp/dz = -6 N/m^3`.
* The ___Volumetric Flow Rate___ through the channel/river is given by 

	`R = (WHVz/)2F_D - (WPzH^3/12mu)F_P`

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
, 
### Introduction

* The Finite Element Method (FEM) is a general numerical method for solving partial differential equations in two- or three-space variables.
* The premise is very simple; continuous domains (geometries) are decomposed into discrete, connected regions (or __finite elements__).
* A typical approach for using FEM involves the following steps:
  	- ___Step 1___: Divide the domain of the problem into a collection of sub-domains (__finite elements__), forming a __mesh__, with each sub-domain represented by a set of __algebraic equations__ for the original differential equation(s).
  	- ___Step 2___: Systematically recombining all sets of the algebraic equations into a __global system of equations__ for the final calculation.
 
* The global system of equations takes the form:

  	__[K]{u} = {F}__,
  where __K__ is the ___stiffness matrix___, __u__ is the __nodal displacement vector__, the unknowns, and __F__ is the ___nodal forces vector___.

* For this project, __u__ represents the velocity field of the fluid and __F__ represents the force due to pressure; i.e. the pressure gradient.
* The stiffness matrix, __K__, in terms of fluid dynamics, represents the relationship between nodal displacements and applied forces, or the resistance of the fluid domain to deformation under external influences.

### Different Methods
* There are several types of Finite Element Method. The most common approaches include, for example:

#### The Direct Stiffness Method

* The __inverse__ of the  __Stiffness Matrix__, `K`, is performed directly; i.e. to find the __fluid flow vector__, __{u} = {F}[K]^-1__.

#### The Weighted Residuals Approach

* E.g. The ___Galerkin Method___

#### The Variational Approach

* e.g. The ___Rayleigh-Ritz Method___

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

* Run the executable:

`$ ./main`


## Getting & Running the Software

* `$ git clone https://github.com/MRLintern/viscousRiver.git`
* `$ cd viscousRiver`
* `$ mkdir build -p && cd build`
* `$ cmake ..`
* `$ cmake --build .`
* `$ ./main`
* `Press 1 for verbose mode`.

