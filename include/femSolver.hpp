// finite element solver variables and functions

#pragma once

#include <memory>
#include <array>
#include <Eigen/Dense>
#include "vector.hpp"

// -- type aliases

// a node is of type int
using node = int;

// an element of the mesh is of type int
using element = int;

// alias for element of the mesh; 3 nodes to 1 element; i.e. a triangular element
using matrixElement = Eigen::Matrix<int, Eigen::Dynamic, 3>;

// -- global variables

// node coordinates
extern int NodeX, NodeY;

// total number of nodes and elements that make up the mesh
extern int totalNodes, totalElements;

// Node pointer; holds coordinates of components
extern std::unique_ptr<Vector[]> Node;

// for tracking elements at nodes
extern matrixElement nodeElement;

// -- functions

// area of element; used for integration; finding volumetric flow rate etc
double area(element e);

// center of mass of elements of mesh
Vector com(element e);
