// single source file; to be broken down into header, source and driver files

#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <vector>
#include <array>
#include <memory>
#include <Eigen/Dense>


class Vector2D {

public:

	// coordinates
	double x{}, y{};

	// initialise coordinates
	Vector2D() {x = 0.0; y = 0.0;};

	// initialise coordinates to those of initial coordinates
	Vector2D(double x0, double y0) {x = x0; y = y0;};

	// set up vector
	void setVector(double x0, double y0) {x = x0; y = y0;};

	// magnitude of vector; 2D Euclidean/Cartesian plane; p-norm, p = 2
	double vectorMag() {

		return sqrt(x*x + y*y);
	}

	// assignment and operator overloading
	double operator*(const Vector2D &otherVector);
	Vector2D operator+(const Vector2D &otherVector);
	Vector2D operator-(const Vector2D &otherVector);

	// for output; initial velocity at (x,y)
	friend std::ostream &operator<<(std::ostream &out, Vector2D v0) {

		out<<v0.x<<" "<<v0.y;

		return out;
	}

};

// overload binary - operator; difference between 2 vectors
Vector2D Vector2D::operator-(const Vector2D &otherVector) {

	return Vector2D(x - otherVector.x, y - otherVector.y);
}

// overload binary + operator; addition of 2 vectors
Vector2D Vector2D::operator+(const Vector2D &otherVector) {

	return Vector2D(x + otherVector.x, y + otherVector.y);
}

// overload * operator for dot product
double Vector2D::operator*(const Vector2D &otherVector) {

	return x*otherVector.x + y*otherVector.y;
}

// overload * operator; lhs = left-hand side; rhs = right-hand side
// lhsVector
Vector2D operator*(const double &lhsVector, const Vector2D &rhsVector) {

	// components of initial velocity vector
	Vector2D v0 = Vector2D(lhsVector*rhsVector.x, lhsVector*rhsVector.y);

	return v0;
}

// same as above but for rhsVector of initial velocity vector
Vector2D operator*(const Vector2D &lhsVector, const double rhsVector) {

	Vector2D v0 = Vector2D(rhsVector*lhsVector.x, rhsVector*lhsVector.y);

	return v0;
}

// overload / operator for initial velocity vector
Vector2D operator/(const Vector2D &lhsVector, const double &rhsVector) {

	Vector2D v0 = Vector2D(lhsVector.x/rhsVector, lhsVector.y/rhsVector);

	return v0;
}

// -- Geometry

// identifier for nodes of type int for the mesh
using node = int;

// identifier for elements of type int for the mesh
using element = int;

// a vector to hold coordinates of nodes
Vector2D *Node{nullptr};

// a custom matrix identifier; the # of rows not known at compile-time; # of columns (3) known at compile-time.
// the matrix identifier allows us to create a variable which will track the nodes associated with elements of the mesh
using matrixElement = Eigen::Matrix<int, Eigen::Dynamic, 3>;

// element which will track nodes associated with mesh elements
matrixElement nodeElement;

// function which calculates area of element(s)
double area(element e) {

	// x components
	double a1 = Node[nodeElement(e, 1)].x - Node[nodeElement(e, 0)].x;
	double a2 = Node[nodeElement(e, 0)].x - Node[nodeElement(e, 2)].x;
	// y components
	double b1 = Node[nodeElement(e, 0)].y - Node[nodeElement(e, 1)].y;
	double b2 = Node[nodeElement(e, 2)].y - Node[nodeElement(e, 0)].y;

	return 0.5*(a1*b2 - a2*b1);
}

// method to calculate the center of mass, com, of the element e
Vector2D com(element e) {

	return (Node[nodeElement(e, 0)] + Node[nodeElement(e, 1)] + Node[nodeElement(e, 2)])/3.0;
}

// all previous code to separated into .hpp files. Now enter main()

int main() {

	// should the user want more messages r.e. program and instructions etc
	int verbose{0}; 
	std::cout<<"Enter 1 for verbose mode:\n";
	std::cin>>verbose;

	// number of nodes in the x & y direction
	int NodeX{}, NodeY{};

	// H/W ratio
	double HW{0.5};

	// pressure gradient
	double Pz{-6.0};

	// intial velocity at the top of the channel
	double Vz{0.0};

	// file objects for output data; file extension: .dat
	std::ofstream meshPts, eleProp, uSol;
	// std::ofstream file3;
	// open files
	meshPts.open("meshPoints.dat"); // mesh coordinates of channel
	eleProp.open("elementProperties.dat"); // elements, center of mass of elements and area of elements
	//file3.open("c3.dat"); 
	uSol.open("uSolution.dat"); // for solution , u, at nodes

	// user to input some data r.e. nodes and elements
	std::cout<<"Enter the Number of Nodes in the x and y direction:\n";
	std::cin>>NodeX>>NodeY;

	// the total number of nodes for the mesh is given by NodeX x NodeY
	int totalNodes {NodeX*NodeY};

	// verbose mode
	std::cout<<"Enter 1 for verbose mode: "<<"\n";
	std::cin>>verbose;

	// total number of elements of the mesh
	int totalElements {(NodeX - 1)*(NodeY - 1)*2}; // a rectangular area

	// display number of nodes used
	std::cout<<"Number of Nodes used: "<<totalNodes<<"\n";

	// display the number of elements used
	std::cout<<"Number of Elements used: "<<totalElements<<"\n";

	// -- Vector & Matrix constructors for FEM; all of these are dynamic-size and of type double
	// ------------------------------------------------------------------------------------------

	Eigen::VectorXd u(totalNodes); // nodal vector; fluid flow rate
	Eigen::VectorXd f(totalNodes); // force acting on fluid; the only force we consider is the contribution of the pressure gradient

	// stiffness matrix; represents the system of linear equations that relates the nodal displacements of the fluid to the forces acting on them
	Eigen::MatrixXd K(totalNodes);

	// -------------------------------------------------------------------------------------------

	// initialise nodeElement with matrix identifier; keep track of nodes associated with elements
	nodeElement = matrixElement(totalElements, 3);

	// -- rectangle of nodes

	// allocate memory for a vector of nodes
	Node = new Vector2D[totalNodes];

	// starting node
	node nodeNumber{0};

	// iterate accross rectangle of nodes for (x, y)
	for (int iy{0}; iy < NodeY; iy++) {
		for (int ix{0}; ix < NodeX; ix++) {

			// make (x, y) values decimals; not going to be whole numbers
			double x {static_cast<double>(ix/NodeX - 1)};
			double y {static_cast<double>(iy/NodeY - 1)};

			// scaling y with H/W ratio
			y *= HW;

			// set points for vector
			Node[nodeNumber].setVector(x, y);
		}
	}

	// label elements and specify their nodes
	if (verbose) {std::cout<<"\n"<<"Element and Associated Nodes: "<<"\n\n";}

	// iterate through the rectangle so we can list the element and its node at (x, y)
	for (int iy{0}; iy < NodeY; iy++) {
		for (int ix{0}; ix < NodeX; ix++) {

			// nodes for (x,y)
			node i {iy*(NodeX - 1) + ix};
			node j {ix + NodeX*iy};

			// element number
			element elementNumber {2*i};

			// node/element positioning
			nodeElement(elementNumber, 0) = j;
			nodeElement(elementNumber, 1) = j + NodeX + 1;
			nodeElement(elementNumber, 2) = j + NodeX;

			// output elements and nodes
			if (verbose) {

				std::cout<<elementNumber<<" "<<nodeElement(elementNumber, 0)<<" "<<nodeElement(elementNumber, 1)<<" "<<nodeElement(elementNumber, 2)<<"\n";

			}

			// update element number
			elementNumber++;

			// now look at node/element position again
			nodeElement(elementNumber, 0) = j;
			nodeElement(elementNumber, 1) = j + 1;
			nodeElement(elementNumber, 2) = j + 1 + NodeX;

			// output elements and nodes
			if (verbose) {

				std::cout<<elementNumber<<" "<<nodeElement(elementNumber, 0)<<" "<<nodeElement(elementNumber, 1)<<" "<<nodeElement(elementNumber, 2)<<"\n";

			}
		}
	}

  if (verbose) {

		std::cout<<" ---------- "<<"\n";

		// print out mesh data to channel.dat file
		for (element e{0}; e < totalElements; e++) {

			meshPts<<Node[nodeElement(e, 0)]<<" "<<Node[nodeElement(e, 1)]<<"\n";
			meshPts<<Node[nodeElement(e, 1)]<<" "<<Node[nodeElement(e, 2)]<<"\n";
			meshPts<<Node[nodeElement(e, 2)]<<" "<<Node[nodeElement(e, 0)]<<"\n";
			// center of mass and area of elements
			eleProp<<e<<" "<<com(e)<<" "<<area(e)<<"\n";

			}
	}

	// see // -- Vector & Matrix constructors for FEM for more info
  // configure stiffness matrix for finite element method (FEM)
  // [K]{u} = {f}
  // K = stiffness matrix; contains the system of linear equations
  // u = nodal vector; fluid flow rate
  // f = force; pressure gradient acting on the fluid
  for (element e{0}; e < totalElements; e++) {

    // these parameters act as weights for calculating the approx. value of K
    // the values these parameters take on can be adjusted to balance accuracy and stability
    std::array<double, 3> beta{};
    std::array<double, 3> gamma{};


    // find beta and gamma
    for (int i{0}; i < 3; i++) {

    	int j {(i + 1) % 3};
    	int k {(i + 1) % 3};

    	beta[i] = Node[nodeElement(e, j)].y - Node[nodeElement(e, k)].y;
    	gamma[i] = Node[nodeElement(e, k)].x - Node[nodeElement(e, j)].x;

    	if (verbose) {std::cout<<e<<" "<<i<<" "<<" "<<beta[i]<<" "<<gamma[i]<<"\n";}

    	}

    	// find K for each 3-nodal point element; triangle in 2D mesh
    	for (int i{0}; i < 3; i++) {
    		for (int j{0}; j < 3; j++) {

    			node I {nodeElement(e, i)};
    			node J {nodeElement(e, j)};

    			// stiffness matrix
    			K(I, J) += (beta[i]*beta[j] + gamma[i]*gamma[j])/(4.0*area(e));
    		}
    	}

    }

    // find the determinant of the stiffness matrix K
    // if det[K] = 0; there are no unique solutions and u can't be found; we don't want this!

    if (verbose) {

    	std::cout<<" ------------ "<<"\n";
    	std::cout<<"K"<<"\n";
    	std::cout<<"det(K) = "<<K.determinant()<<"\n";

    	// value of stiffness matrix at (i, j)
    	for (node i{0}; i < totalNodes; i++) {
    		for (node j{0}; j < totalNodes; j++) {

    			std::cout<<i<<" "<<j<<" "<<K(i,j)<<"\n";
    		}
    	}

    	std::cout<<" ------------- "<<"\n";
    }

    // assemble force vector; i.e. the pressure gradient acting on the fluid; Pz = (1/mu)*dp/dz
    for (element e{0}; e < totalElements; e++) {

    	// holds total force contributions on elements
    	double F {-Pz*area(e)/3.0};

    	// total force for each element
    	f(nodeElement(e, 0)) += F;
    	f(nodeElement(e, 1)) += F;
    	f(nodeElement(e, 2)) += F;

    }

    // display force at nodes

    if (verbose) {

    	std::cout<<"f"<<"\n";

    	for (node i{0}; i < totalNodes; i++) {

    		std::cout<<i<<" "<<f(i)<<"\n";
    	}

    	std::cout<<" --------------- "<<"\n";

    }

    // -- Boundary Conditions (BCs) (Dirichlet); u = f = 0

    // nodes on the boundary of the domain: top & bottom of channel
    for (node n{0}; n < NodeX; n++) {

    	if (verbose) {std::cout<<"Boundary Nodes: "<<n<<" "<<n + NodeX*(NodeY - 1)<<"\n";}

    	// stiffness matrix for BCs
    	K(n, n) *= 1e12;
    	K(n + NodeX*(NodeY - 1), n + NodeX*(NodeY - 1)) *= 1e12;

    	// -- applied force/pressure gradient
    	// bottom of channel
    	f(n) = 0.0*K(n, n);

    	// top of channel
    	f(n + NodeX*(NodeY - 1)) = Vz*K(n + NodeX*(NodeY - 1), n + NodeX*(NodeY - 1));

    }

    // nodes on the boundary of the domain: left & right of channel
    for (node n{1}; n < NodeY - 1; n++) {

    	if (verbose) {std::cout<<"Boundary Nodes: "<<n*NodeX<<" "<<n*NodeX + NodeX - 1<<"\n";}

    	// stiffness matrix
    	K(n*NodeX, n*NodeX) *= 1e12;
    	K(n*NodeX + NodeX - 1, n*NodeX + NodeX - 1) *= 1e12;

    	// -- applied force/pressure gradient
    	// left hand side of channel
    	f(n*NodeX) = 0.0*K(n*NodeX, n*NodeX);

    	// right hand side of channel
    	f(n*NodeX + NodeX - 1) = 0.0;

    }

    // find non-zero entries of the stiffness matrix
    int num{0};

    // total K entries providing K isn't 0
    for (node n{0}; n < totalNodes; n++) {
    	for (node m{0}; m < totalNodes; m++) {

    		if (K(n, m) != 0.0) num++;
    	}
    }

    // print out non-zero entries of K
    std::cout<< --num<<"Non-Zero Elements in K"<<"\n";

    // find the fluid flow velocity by inverting the stiffness matrix
    // Note: {u} = {F}[K]^-1
    u = K.inverse()*f;

    // -- store fluid flow velocity solution u
    
    int nx{1};

    // solution u for each node stored in file
    for (node n{0}; n < totalNodes; n++) {

    	uSol<<Node[n]<<" "<<u(n)<<"\n";

    	nx++;

    	if (nx > NodeX) {

    		uSol<<" "<<"\n";

    		int nx{1}; // apparently the type has to be provided again? 
    	}
    }

    // -- now find the flow rate through the channel; integration

    // flow rate 
    double ut{0.0};

    // find flow rate accross all elements
    for (element e{0}; e < totalElements; e++) {

    	ut += area(e)*(u(nodeElement(e, 0)) + u(nodeElement(e, 1)) + u(nodeElement(e, 2)));
    }

    ut /= 3.0;

    // display result
    std::cout<<"Flow Rate (ut) = "<<ut<<"\n";

    // -- free up memory and close files

    delete[] Node;
    meshPts.close();
    eleProp.close();
    uSol.close();

    // messages to explain files
    std::cout<<"Mesh Coordinates in channel.dat"<<"\n";
    std::cout<<"totalElements, Center of Mass (of Elements) and Area of Elements in file c2.dat"<<"\n";
    std::cout<<"Nodes and Velocity at Nodes in file c4.dat"<<"\n";

}
