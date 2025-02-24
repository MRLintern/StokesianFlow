// single source file; to be broken down into header, source and driver files

#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <vector>
#include <array>
#include <Eigen/Dense>

// -- Vector.hpp

// vector class for vectors & associated objects
class Vector {

    public:

    // coordinates
    double x{}, y{};

    // initialise coordinates of vector object
    Vector () {x = 0.0; y = 0.0;}

    // initialize coordinates to initial coordinate values
    Vector (double x0, double y0) {x = x0; y = y0;}

    // setup vector object method
    void setVector (double x0, double y0) {x = x0; y = y0;}

    // metric: distance between points
    double distance() {

        return sqrt(x*x + y*y);
    }

    // -- assignment and operator overloading

    double operator*(const Vector &otherVector);
    Vector operator+(const Vector &otherVector);
    Vector operator-(const Vector &otherVector);

    // for output; initial velocity at (x,y)
	friend std::ostream &operator<<(std::ostream &out, Vector v0) {

		out<<v0.x<<" "<<v0.y;

		return out;
	}

}; // end of Vector class; to be placed in Vector.hpp

// -- Vector.cpp

// overload binary - operator; difference between 2 vectors
Vector Vector::operator-(const Vector &otherVector) {

    return Vector(x - otherVector.x, y - otherVector.y);
}

// overload binary + operator; addition of 2 vectors
Vector Vector::operator+(const Vector &otherVector) {

    return Vector(x + otherVector.x, y + otherVector.y);
}

// overload * operator for dot product
double Vector::operator*(const Vector &otherVector) {

    return x*otherVector.x + y*otherVector.y;
}

// overload * operator for v0; lhs = left-hand side; rhs = right-hand side
// lhsVector
Vector operator*(const double &lhsVector, const Vector &rhsVector) {

    Vector v0 = Vector(lhsVector*rhsVector.x, lhsVector*rhsVector.y);

    return v0; 
}

// overload * operator for v0; lhs = left-hand side; rhs = right-hand side
// rhsVector
Vector operator*(const Vector &lhsVector, const double &rhsVector) {

    Vector v0 = Vector(rhsVector*lhsVector.x, rhsVector*lhsVector.y);

    return v0;
}

// overload / operator for v0
Vector operator/(const Vector &lhsVector, const double &rhsVector) {

    Vector v0 = Vector(lhsVector.x/rhsVector, lhsVector.y/rhsVector);

    return v0;
}

// ------------------------------------------------------------------------------------------------

// -- Geometry.cpp

// identifier for nodes of type int for the mesh
using node = int;

// identifier for elements of type int for the mesh
using element = int;

// vector to hold coordinates of nodes
Vector *Node {nullptr};

// a custom matrix identifier; the # of rows not known at compile-time; # of columns (3) known at compile-time.
// the matrix identifier allows us to create a variable which will track the nodes associated with elements of the mesh
using matrixElement = Eigen::Matrix<int, Eigen::Dynamic, 3>;

// element which will track nodes associated with mesh elements
matrixElement nodeElement;

// function which calculates area of element(s)
double area(element e) {

    // x component
    double a1 {Node[nodeElement(e, 1)].x - Node[nodeElement(e, 0)].x};
    double a2 {Node[nodeElement(e, 0)].x - Node[nodeElement(e, 2)].x};

    // y component
    double b1 {Node[nodeElement(e, 0)].y - Node[nodeElement(e, 1)].y};
    double b2 {Node[nodeElement(e, 2)].y - Node[nodeElement(e, 0)].y};

    return 0.5*(a1*b2 - a2*b1);
}

// function calculates center of mass of elements
Vector com(element e) {

    return (Node[nodeElement(e, 0)] + Node[nodeElement(e, 1)] + Node[nodeElement(e, 2)])/3.0;

}

// main.cpp

int main() {

    // verbose: this setting will provide the option do display more data
    int verbose {0};

    std::cout<<"Do want to use verbose mode? Enter 1 for yes, 0 for no\n";
    std::cin>>verbose;

    // nodes in x and y directions
    int NodeX {6}, NodeY {6}; // 6 x 6 mesh 

    // H/W ratio
    double HW {0.5};

    // pressure gradient Pz = (1/mu) * dp/dz
    double Pz {-6};

    // initial velocity at the top of the channel
    double Vz {0.0};

    // file objects
    std::ofstream meshPts, eleProps, uSol;

    // files to open for calculated data
    meshPts.open("meshPoints.dat"); // coordinates in the mesh
    eleProps.open("elementProps.dat"); // element properties; com, area etc
    uSol.open("uSolution.dat"); // solution, u, at nodes in mesh

    // find total number of nodes used
    int totalNodes {NodeX*NodeY};

    // find total number of elements used
    int totalElements {(NodeX - 1)*(NodeY - 1)*2};

    // output to screen total number of nodes & elements used
    std::cout<<"Number of nodes used: "<<totalNodes<<"\n";
    std::cout<<"Number of elements used: "<<totalElements<<"\n";

    // -- Vector & Matrix constructors for FEM; all of these are dynamic-size and of type double
	// ------------------------------------------------------------------------------------------

    // f: force acting on fluid; the only force we consider is the contribution of the pressure gradient
    // u: fluid flow velocity
    Eigen::VectorXd u(totalNodes), f(totalNodes);

    // stiffness matrix; represents the system of linear equations that relates the nodal displacements of the fluid to the forces acting on them
    Eigen::MatrixXd K(totalNodes, totalNodes);

    // ------------------------------------------------------------------------------------------------------------------------

    // initialise nodeElement with matrix identifier; keep track of nodes associated with elements
    nodeElement = matrixElement(totalElements, 3);

    // -- rectangle of nodes

    // allocate memory for a vector of nodes
    Node = new Vector[totalNodes];

    // node number; starting position
    node nodeNumber{0};

    // iterate across rectangle of nodes for (x,y)
    for (int iy{0}; iy < NodeY; iy++) {
        for (int ix{0}; ix < NodeX; ix++) {

            // make (x,y) values decimals; discretization leaves numbers as decimals
            double x {static_cast<double>(ix/NodeX - 1)};
            double y {static_cast<double>(iy/NodeY - 1)};

            // scaling y with H/W ratio
            y *= HW;

            // set node numbers with vector points x & y
            Node[nodeNumber].setVector(x, y);

            // keep adding node numbers until rectangle is full
            nodeNumber++;
        }
    }

    // label elements and specify their nodes
    if (verbose) {std::cout<<"\n"<<"Element and Associated Nodes: "<<"\n\n";}

    // iterate through the rectangle so we can list the elements and their nodes at (x,y)
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

            // update element numbers and nodes
            elementNumber++;
            nodeElement(elementNumber, 0) = j;
            nodeElement(elementNumber, 1) = j + 1;
            nodeElement(elementNumber, 2) = j + 1 + NodeX;

            // output updated element numbers and associated nodes
            if (verbose) {

                std::cout<<elementNumber<<" "<<nodeElement(elementNumber, 0)<<" "<<nodeElement(elementNumber, 1)<<" "<<nodeElement(elementNumber, 2)<<"\n";

            }
        }
    }

    // data to files
    if (verbose) {

        std::cout<<" ------------------ "<<"\n";

        // print out mesh points (nodes) to meshPts
        for (element e{0}; e < totalElements; e++) {

            meshPts<<Node[nodeElement(e, 0)]<<" "<<Node[nodeElement(e, 1)]<<"\n";
            meshPts<<Node[nodeElement(e, 1)]<<" "<<Node[nodeElement(e, 2)]<<"\n";
            meshPts<<Node[nodeElement(e, 2)]<<" "<<Node[nodeElement(e, 0)]<<"\n";
            // center of mass, com, and area of elements
            eleProps<<e<<" "<<com(e)<<" "<<area(e)<<"\n";

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

        // find beta and gamma values
        for (int i{0}; i < 3; i++) {

            int j {(i + 1) % 3};
            int k {(i + 1) % 3};

            beta[i] = Node[nodeElement(e, j)].y - Node[nodeElement(e, k)].y;
            gamma[i] = Node[nodeElement(e, k)].x - Node[nodeElement(e, j)].x;

            // print out values for beta, gamma associated element
            if (verbose) {std::cout<<e<<" "<<i<<" "<<beta[i]<<" "<<gamma[i]<<"\n";}
        }

        // find K for each 3-nodal point element; triangle in 2D mesh
        for (int i{0}; i < 3; i++) {
            for (int j{0}; j < 3; j++) {

                node I {nodeElement(e, i)};
                node J {nodeElement(e, j)};

                // stiffness matrix K
                K(I, J) += (beta[i]*beta[j] + gamma[i]*gamma[j])/(4.0*area(e));
            }
        }
    }

    // find the determinant of the stiffness matrix K
    // if det[K] = 0, => there are no unique solutions and u can't be found

    if (verbose) {

        std::cout<<" -------------- "<<"\n";
        std::cout<<"K"<<"\n";
        std::cout<<"det(k) = "<<K.determinant()<<"\n";

        // value of K at (i,j)
        for (node i{0}; i < totalNodes; i++) {
            for (node j{0}; j < totalNodes; j++) {

                std::cout<<i<<" "<<j<<" "<<K(i, j)<<"\n";
            }
        }

        std::cout<<" --------------- "<<"\n";
    }

    // assemble force vector: i.e. the pressure gradient acting on the fluid; Pz = (1/mu)*dp/dp
    for (element e{0}; e < totalElements; e++) {

        // total force contribution on elements
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

        std::cout<<" ------------------- "<<"\n";
    }

    // -- Boundary conditions, BCs: Dirichlet; u = f = 0

    // nodes on the boundary of the domain: top & bottom of channel
    for (node n{0}; n < NodeX; n++) {

        if (verbose) {std::cout<<"Boundary Nodes: "<<n<<" "<<n + NodeX*(NodeY - 1)<<"\n";}

        // stiffness matrix for BCs
        K(n, n) *= 1e12;
        K(n + NodeX*(NodeY - 1), n + NodeX*(NodeY - 1)) *= 1e12;

        // -- applied force/pressure gradient
        // left hand side of channel
        f(n*NodeX) = 0.0*K(n*NodeX, n*NodeX);

        // right hand side of channel
        f(n*NodeX + NodeX - 1) = 0.0;
    }

    // -- find non-zero entries of the stiffness matrix
    // holds entries of K matrix
    int num{0};

    for (node n{0}; n < totalNodes; n++) {
        for (node m{0}; m < totalNodes; m++) {

            // add up non-zero entries for K
            if (K(n, m) != 0.0) num++;
        }
    }

    // print out non-zero entries for K
    std::cout<< --num<<"Non-Zero Elements in K"<<"\n";

    // find the fluid flow velocity u by inverting K
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

            nx = 1;
        }
    }

    // -- now find the flow rate through the channel; integration

    // flow rate
    double ut{0.0};

    // find flow rate across all elements
    for (element e{0}; e < totalElements; e++) {

        ut += area(e)*(u(nodeElement(e, 0)) + u(nodeElement(e, 1)) + u(nodeElement(e, 2)));

    }

    ut /= 3.0;

    // display flow rate to screen
    std::cout<<"Flow Rate = "<<ut<<"\n";

    // -- free up memory and close files

    delete[] Node;
    meshPts.close();
    eleProps.close();
    uSol.close();

    // output messages to explain files
    std::cout<<"Mesh Coordinates stored in meshPoints.dat"<<"\n";
    std::cout<<"Total Elements, Center of Mass (of Elements) and Area of Elements stored inside elementProps.dat"<<"\n";
    std::cout<<"Nodes and Velocity at Nodes stored inside uSolution.dat"<<"\n";

    // end of program

}








