## StokesianFlow

* Note: this project makes use of a lot material from the book ___Finite Element Methods for Engineers___ by ___Roger T. Fenner___.
* Fenner provides a `Fortran 77` solver for the problem we're looking at.
* My solution to the case study involves refactoring the software into __C++14__, utilising __OOD__ principles, __Modern C++ Memory Management__ and the __Eigen Template Library__.

## Introduction

* This software, developed in Modern C++, models a ___Stokes Flow___, a fluid dominated by __Viscosity__ rather than __Diffusion__.
* The domain is `discretised` via the `Finite Element Method`, `(FEM)`. The resulting `System of Algebraic Equations` is handled and solved via the `Eigen C++ template library`.
* This software calculates the __Velocity Profile__ and using this, the __Volumetric Flow Rate__ through the channel.

## Model

* The model consists of a `channel` or rather, a `river`, flowing along the `z-axis (direction)`.
* Note: all physical quantities use `S. I. units`.
* The fluid flows with `uniform velocity`, and has `no external forces` acting on it. The only force acting on the flow is the __Pressure Gradient__.
* If the rectangular channel is infinitely wide compared to its depth, `H << W`, the governing ODE is given by:

	`d^2 u/dy^2 = (1/mu)dp/dz`

* `d^2 u/d^2`: Second order derivative of the __Velcoity Field__.

* `mu`: `Dynamic (Shear) Viscosity`: Dynamic viscosity is a fluid's resistance to flow and shear when an external force is applied. Taken to be constant.

* `dp/dz`: This is the `Pressure Gradient`: the rate at which pressure changes over distance. Taken to be constant. In the software, this represented by `Pz`.

* Integrating the ODE with Boundary Conditions (BC, see below) to find the __Velocity Profile__ gives:

  	`u = yVz/H + Pz/2mu(y^2 - yH)`

* __Dirichlet Boundary Conditions (No-Slip)__ (BCs): At the bottom and sides of the channel/river: `v = 0 m/s`. At the top (surface): fixed `v = Vz`.
* The shape of the `domain` is a rectangle of `height H` and `width W`. Values for these: `H = 1m` and `W = 2m`.
* `Dynamic Viscosity, mu = 1 N.s/m^2`.
* `Pressure Gradient, dp/dz = -6 N/m^3`.
* The ___Volumetric Flow Rate___ through the channel/river is given by 

	`R = (WHVz/)2F_D - (WPzH^3/12mu)F_P`

* `F_D`: ___Drag Factor___; associated with the drag flow induced by the boundary velocity `Vz`.
* `F_P`: ___Pressure Flow Shape___; due to the pressure gradient `p'` or `dp/dz`. 
* The values of `F_D` & `F_P` provide an important ratio in modelling channel flow problems: the `aspect ratio`, `H/W`.
* Because the `H/W` ratio is provided, `F_D` & `F_P` aren't calculated and included in the solver; they're just mentioned for the purposes of interest
* Both the `F_D` & `F_P` factors are __Hyperbolic Functions__. These are complicated symbolically, and so are left out. See ___Finite Element Methods for Engineers, Fenner 1974___.
* `H/W = 0.5` for this problem.
* Properties:

	- __Shallow and wider channels__: larger wetted perimeter leading to increased friction and slower flow velocities.
	- __Deeper and narrower channels__: smaller wetted perimeter resulting in less friction and faster flow velocities.
* The number of nodes (___vertices___) in the `x-direction` is given by `NodeX`. Those in the `y-direction`: `NodeY`.
* The __mesh__ of the problem domain consists of `6 x 6 nodal points`.
* The number of __triangular elements__ that make up the mesh is denoted by `elementNumber`. 50 will be used in total.

## Finite Element Method (FEM)

* This is a brief overview of what the method entails.
* To understand the theory of the Finite Element Method applied to Stokes Flow, click here [View PDF](./Fem-Theory.pdf)
, 
### Introduction

* The Finite Element Method (FEM) is a general numerical method for solving partial differential equations in two- or three-space variables.
* The premise is very simple; continuous domains (geometries) are decomposed into discrete, connected regions (or __finite elements__).
* A typical approach for using FEM involves the following steps:
  	- ___Step 1___: Divide the domain of the problem into a collection of sub-domains (__finite elements__), forming a __mesh__, with each sub-domain represented by a system of __Linear Algebraic Equations__ for the original differential equation(s).
  	- ___Step 2___: Systematically recombining all sets of the algebraic equations into a __global system of equations__ for the final calculation.
 
* The global system of equations takes the form:

  	__[K]{u} = {F}__,
  where __K__ is the ___Stiffness Matrix___, __u__ is the __Nodal Displacement Vector__, the unknowns, and __F__ is the ___Nodal Forces Vector___.

* For this project, __u__ represents the velocity field of the fluid and __F__ represents the force due to pressure; i.e. the pressure gradient.
* The stiffness matrix, __K__, in terms of fluid dynamics, represents the relationship between nodal displacements and applied forces, or the resistance of the fluid domain to deformation under external influences.

### Different Methods
* There are several types of Finite Element Method. The most common approaches include, for example:

#### The Direct Stiffness Method

* The __inverse__ of the  __Stiffness Matrix__, `K`, is performed directly; i.e. to find the __Fluid flow Vector__, __{u} = {F}[K]^-1__. This is the approach used in ___this project___.
  	 #### 1. Element Stiffness Matrix
  	* Each __Finite Element__ (e.g., a bar, beam, or truss) has its own ___Local Stiffness Matrix___ `k`, which relates ___Nodal Displacements___ `u` to forces `f`: `f = ke`.
  	 #### 2. Assembly into Global System
  	* All ___Local Element Matrices___ are __assembled__ into a ___Global Stiffness Matrix___ `K` using __connectivity information__ (mapping __local DOFs__ to __global DOFs__). Note: `DOF = Degrees of Freedom`.
  	* The ___Global Force-Displacement System___: `Ku = F`, where `u` is the __vector__ of ___Global Displacements___, and `F` is the __vector__ of ___Global Applied Forces___.
  	* Note: in this project, `u` represents the `Fluid Velocity` and `F` represents the `Pressure Gradient`.
  	 #### 3. Boundary Conditions
  	* Essential (___Dirichlet___) ___Boundary Conditions___ are applied by modifying `K` and `F` appropriately.
  	 #### 4. Solve for Displacements
  	* After assembling the system and applying boundary conditions, solve the resulting linear system for `u`.
  	 #### 5. Compute Internal Forces
  	* Once ___Nodal Displacements___ are known, ___Internal Element Forces___ or ___Stresses___ can be computed using each ___Element’s Stiffness Relation___.   

#### The Weighted Residuals Approach

* E.g. The ___Galerkin Method___.
* __Foundation__: Derived from the method of ___Weighted Residuals___.
* __Key Idea__: The error (residual) between the ___exact___ and ___approximate solution___ is made __orthogonal__ to a set of ___Weight (test) Functions___.
* __Applications__: Widely used for both ___Linear___ and ___Nonlinear PDEs___ including __Heat Transfer__, __Elasticity__, and __Fluid Flow__.
* __Key Feature__: Works even when no variational (energy) principle exists; flexible and broadly applicable.

#### The Variational Approach

* e.g. The ___Rayleigh-Ritz Method___.
* __Foundation__: Based on the ___calculus of variations___ and ___Energy Minimization___.
* __Key Idea__: The solution is found by ___Minimizing the total potential energy___ __or__ a ___Related Functional___.
* __Requirements__: The __problem__ must be ___Self-Adjoint___ and have a ___Variational (energy) Formulation___.
* __Applications__: __Structural Mechanics__, __Vibration Analysis__, and __Problems with Conservative Forces__.
* __Key Feature__: Naturally provides ___Symmetric___ and ___Positive-Definite Stiffness Matrices___.

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

* `~/Downloads$ tar -xf eigen-3.4.0.tar.gz`
* `~/Downloads$ mv eigen-3.4.0 /home/user`

* You'll need to make a `build directory` for `CMake`. 

* `~$ cd eigen-3.4.0`
* `~/eigen-3.4.0$ mkdir build && cd build`

* Run `CMake`

* `~/eigen-3.4.0/build$ cmake ..`

* Now install `Eigen`

* `~/eigen-3.4.0/build$ sudo make install`

* symlink or copy the Eigen folder into `/usr/local/include/`

* `~/eigen-3.4.0/build$ sudo cp -r /usr/local/include/eigen3/Eigen /usr/local/include`

* And you're set! See the test program in `eigen_test` to ensure everything is working.

* Inside the `eigen_test` folder, execute the following at the CLI:

* `$ g++ main.cpp -o main`

* Run the executable:

* `$ ./main`


## Getting & Running the Software

* `$ git clone https://github.com/MRLintern/StokesianFlow.git`
* `$ cd StokesianFlow`
* `$ mkdir build -p && cd build`
* `$ cmake ..`
* `$ cmake --build .`
* `$ ./StokesianFlow`
* `Press 1 for verbose mode`.

## References
* ___Finite Element Methods for Engineers___ by ___Roger T. Fenner___.


## Results

* The __Volumetric Flow Rate__ came in at around `0.031 m^3/s`.
* The results will be generated and found in the `build` directory when you run the software. However, I have made a directory called `Results` where you can view the generated files along with a `Python script` for plotting and a `png` image file of the fluid flow.
* The __Velcoity Profile__ is provided in the file `uSolution.dat`. Instead of using the `Python script`, [ParaView](https://www.paraview.org/) can be used.
