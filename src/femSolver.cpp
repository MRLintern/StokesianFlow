// implement femSolver functionality

#include <memory>
#include "femSolver.hpp"

// -- make a 6 x 6 mesh; 36 triangular elements

// number of nodes in x and y-directions respectively
int NodeX{6}, NodeY{6};

// total number of nodes
int totalNodes { (NodeX) * (NodeY) };

// total number of elements making up the mesh
int totalElements { (NodeX - 1) * (NodeY - 1) * 2 };

// allocate memory for coordinates of nodes
std::unique_ptr<Vector[]> Node { std::make_unique<Vector[]>(totalNodes) };

// track node associated with an element
matrixElement nodeElement;

// calculate the area of triangular element; uses determinant formula
// Note: numbers 0, 1, 2: represent vertices of triangular element; e is the element index
double area(element e) {

	// -- x-components

	double a1 { Node[nodeElement(e, 1)].getX() - Node[nodeElement(e, 0)].getX() }; // a1 = x1 - x0
	double a2 { Node[nodeElement(e, 0)].getX() - Node[nodeElement(e, 2)].getX() }; // a2 = x0 - x2

	// -- y-component

	double b1 { Node[nodeElement(e, 0)].getY() - Node[nodeElement(e, 1)].getY() }; // b1 = y0 - y1
	double b2 { Node[nodeElement(e, 2)].getY() - Node[nodeElement(e, 0)].getY() }; // b2 = y2 - y0

	return 0.5 * (a1 * b2 - a2 * b1);
}

// calculate center of mass (centroid) of a triangular element
// Note: each element is assumed to have uniform density
// nodeElement(e, i): get the global node index of the ith vertex of element e
// Node[index]: get coordinates of the node; represented as a Vector object

Vector com(element e) {

	return (Node[nodeElement(e, 0)] + Node[nodeElement(e, 1)] + Node[nodeElement(e, 2)]) / 3.0;
}
