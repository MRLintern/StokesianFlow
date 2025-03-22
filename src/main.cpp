#include <iostream>
#include <cmath>
#include <fstream>
#include <random>
#include <vector>
#include <array>
#include <memory>
#include <Eigen/Dense>

// Vector.hpp

class Vector {

    public:

        // coordinates
        double x{}, y{};

        // initialise coordinates
        Vector() { x = 0.0; y = 0.0; };

        // initialise coordinates to initial coordinate values
        Vector (double x0, double y0) { x = x0; y = y0; };

        // set up vector
        void setVector(double x0, double y0) { x = x0; y = y0; }

        // size of the vector; p-norm; p = 2 for this problem; Euclidean norm
        double norm() {

            return sqrt(x*x + y*y);
        }

        // -- assignment and operator overloading
        double operator*(const Vector& otherVector);
        Vector operator+(const Vector& otherVector);
        Vector operator-(const Vector& otherVector);

        // for output; initial velocity at (x, y)
        friend std::ostream& operator<<(std::ostream& output, Vector v0) {

            output<<v0.x<<" "<<v0.y;

            return output;
        }
        
};

// Vector.cpp

// overload binary - operator; difference between 2 vectors
Vector Vector::operator-(const Vector& otherVector) {

    return Vector(x - otherVector.x, y - otherVector.y);
}

// overload binary + operator; addition of 2 vectors
Vector Vector::operator+(const Vector& otherVector) {

    return Vector(x + otherVector.x, y + otherVector.y);
}

// overload * operator for dot product
double Vector::operator*(const Vector& otherVector) {

    // dot product
    return x*otherVector.x + y*otherVector.y;
}

// overload * operator; left hand side vector 
Vector operator*(const double &lhsVector, const Vector& rhsVector) {

    // components of initial velocity vector
    Vector v0 = Vector(lhsVector*rhsVector.x, lhsVector*rhsVector.y);

    return v0;
}

// overload * operator; right hand side vector
Vector operator*(const Vector& lhsVector, const double &rhsVector) {

    // components of initial velocity vector
    Vector v0 = Vector(rhsVector*lhsVector.x, rhsVector*lhsVector.y);

    return v0;
}

// overload / operator for initial velocity vector
Vector operator/(const Vector& lhsVector, const double &rhsVector) {

    // components of initial velocity vector
    Vector v0 = Vector(lhsVector.x/rhsVector, lhsVector.y/rhsVector);

    return v0;
}

// end of Vector.cpp

// ------------------------------------------------------------------------------------------------

// Geometry.cpp

// identifier for nodes of type int for the mesh
using node = int;

// identifier for elements of type int for the mesh
using element = int;

// total number of elements that make up the mesh
int totalElements {(NodeX - 1)*(NodeY - 1)*2};

// number of nodes in the x & y direction
int NodeX {6}, NodeY {6};

// total number of nodes used for mesh
int totalNodes {NodeX*NodeY};

// vector to hold coordinates of nodes
std::unique_ptr<Vector[]> Node = std::make_unique<Vector[]>(totalNodes);

// a custom matrix identifier; # of rows at compile-time; # of columns (3) at compile-time
// the matrix identifier allows us to create a variable which will track the nodes associated with elements of the mesh
using matrixElement = Eigen::Matrix<int, Eigen::Dynamic, 3>;

// element which will track nodes associated with mesh elements
matrixElement nodeElement;

// function calculates the area (dx*dy) of elements making up the mesh
// this is used for integrating the equation governing fluid flow rate
// through the (elements that make up the) channel/river
double area(element e) {

    // x-components
    double a1 {Node[nodeElement(e, 1)].x - Node[nodeElement(e, 0)].x};
    double a2 {Node[nodeElement(e, 0)].x - Node[nodeElement(e, 2)].x};

    // y-components
    double b1 {Node[nodeElement(e, 0)].y - Node[nodeElement(e, 1)].y};
    double b2 {Node[nodeElement(e, 2)].y - Node[nodeElement(e, 0)].y};

    return 0.5*(a1*b2 - a2*b1);
}

// method to calculate the center of mass, com, of the elements e
Vector com(element e) {

    return (Node[nodeElement(e, 0)] + Node[nodeElement(e, 1)] + Node[nodeElement(e, 2)])/3.0;
}

// end of Geometry.cpp

// ----------------------------------------------------------------

// main

int main() {

    // allows data to be printed to the CLI; always choose 1!
    int verbose {0};

    
    // H/W ratio
    double HW {0.5};

    // Pressure gradient; the only force that acts on the fluid
    double Pz {-6};

    // initial velocity at the top of the channel
    double Vz {0.0};

    // file objects for output data; file extension: .dat
    std::ofstream meshPts, eProps, uSol;

    // open files
    meshPts.open("meshPoints.dat");
    eProps.open("elementProperties.dat");
    uSol.open("uSolution.dat");

    std::cout<<"Enter 1 for verbose mode, 0 if not\n";
    std::cin>>verbose;

    // output total number of nodes and elements to the CLI
    std::cout<<"Total Number of Nodes: "<<totalNodes<<"\n";
    std::cout<<"Total Number of Elements: "<<totalElements<<"\n";

    // -- Vector & Matrix Constructors for FEM: all of these are dynamic-size and of type double
    // --------------------------------------------------------------------------------------------------------------------

    // nodal vector; force acting on fluid; pressure gradient
    Eigen::VectorXd f(totalNodes);

    // nodal vector; fluid velocity
    Eigen::VectorXd u(totalNodes);

    // stiffness matrix K; represents the system of linear equations that relates the nodal displacements of the fluid to the forces acting on them
    Eigen::MatrixXd K(totalNodes, totalNodes);

    // element properties for nodes
    nodeElement = matrixElement(totalElements, 3);

    // -- rectangle of nodes

    // allocate memory for a vector of nodes
    //Node = new Vector[totalNodes];

    // starting node
    node nodeNumber {0};

    // iterate across rectangle of nodes for (x, y)
    for (int iy {0}; iy < NodeY; iy++) {
        for (int ix {0}; ix < NodeX; ix++) {

            // coordinates in mesh; convert ints to doubles
            double x {static_cast<double>(ix)/(NodeX - 1)};
            double y {static_cast<double>(iy)/(NodeY - 1)};
            
            // scaling y with H/W ratio
            y *= HW;

            // set points for vector
            Node[nodeNumber].setVector(x, y);

            // update node number
            nodeNumber++;
        }
    }

    // label elements and specify their nodes
    if (verbose) {std::cout<<"Element and Associated Nodes: "<<"\n";}

    // iterate through the rectangle so we can list the elements their nodes at (x, y)
    for (int iy {0}; iy < NodeY - 1; iy++) {
        for (int ix {0}; ix < NodeX - 1; ix++) {

            // node x direction
            node i {iy*(NodeX - 1) + ix};

            // element number 
            element eNumber = 2*i;

            // node y direction
            node j {ix + NodeX*iy};

            // node/element position
            nodeElement(eNumber, 0) = j;
            nodeElement(eNumber, 1) = j + NodeX + 1;
            nodeElement(eNumber, 2) = j + NodeX;

            // output elements and nodes
            if (verbose) {std::cout<<eNumber<<" "<<nodeElement(eNumber, 0)<<" "<<nodeElement(eNumber, 1)<<" "<<nodeElement(eNumber, 2)<<"\n";}

            // update element number
            eNumber++;

            // now consider node/element position after updating
            nodeElement(eNumber, 0) = j;
            nodeElement(eNumber, 1) = j + 1;
            nodeElement(eNumber, 2) = j + 1 + NodeX;

            // output updated elements and nodes
            if (verbose) {std::cout<<eNumber<<" "<<nodeElement(eNumber, 0)<<" "<<nodeElement(eNumber, 1)<<" "<<nodeElement(eNumber, 2)<<"\n";}

        }
    }

    if (verbose) {

        std::cout<<" ----------------- "<<"\n";

        // print out mesh data to meshPoints.dat
        for (element e {0}; e < totalElements; e++) {

            meshPts<<Node[nodeElement(e, 0)]<<" "<<Node[nodeElement(e, 1)]<<"\n";
            meshPts<<Node[nodeElement(e, 1)]<<" "<<Node[nodeElement(e, 2)]<<"\n";
            meshPts<<Node[nodeElement(e, 2)]<<" "<<Node[nodeElement(e, 0)]<<"\n";

            // print out element, center of mass and area of elements to elementProperties.dat
            eProps<<e<<" "<<com(e)<<" "<<area(e)<<"\n";
        }
    }

    // see Vector & Matrix Constructors for FEM for my info
    // stiffness matrix K
    // [K]{u} = {f}
    // [K] = stiffness matrix; contains the system of linear equations
    // {u} = fluid velocity
    // {f} = force; pressure gradient acting on fluid
    for (element e {0}; e < totalElements; e++) {
        
        // beta & gamma; shape functions
        // - functions which interpolate the solution between the discrete values obtained at the mesh nodes.
        std::array<double, 3> beta{};
        std::array<double, 3> gamma{};

        // find beta and gamma for K
        for (int i = 0; i < 3; i++) {

            int j {(i + 1) % 3};
            int k {(i + 2) % 3};

            beta[i] = Node[nodeElement(e, j)].y - Node[nodeElement(e, k)].y;
            gamma[i] = Node[nodeElement(e, k)].x - Node[nodeElement(e, j)].x;

            // display element, beta and gamma values
            if (verbose) {std::cout<<e<<" "<<i<<" "<<beta[i]<<" "<<gamma[i]<<"\n";}

        }

        for (int i {0}; i < 3; i++) {
            for (int j {0}; j < 3; j++) {

                node I {nodeElement(e, i)};
                node J {nodeElement(e, j)};

                // find K at (I, J)
                K(I, J) += (beta[i]*beta[j] + gamma[i]*gamma[j])/(4.0*area(e));
            }
        }
    }

    // find the determinant of the stiffness matrix K
    // Note: if | K | = 0; there are no unique solutions and u can't be found
    if (verbose) {

        std::cout<<" --------------- "<<"\n";
        std::cout<<"K "<<"\n";

        // method (from Eigen) finds the determinant
        std::cout<<"det(K) = "<<K.determinant()<<"\n";

        // value of stiffness matrix K at (i, j)
        for (node i {0}; i < totalNodes; i++) {
            for (node j {0}; j < totalNodes; j++) {

                std::cout<<i<<" "<<j<<" "<<K(i, j)<<"\n";
            }
        }

        std::cout<<" --------------- "<<"\n";
    }

    // construct the force vector; pressure gradient
    for (element e {0}; e < totalElements; e++) {

        // holds total force contributions on elements
        double F {-Pz*area(e)/3.0};

        // forces for triangular elements/vertices
        f(nodeElement(e, 0)) += F;
        f(nodeElement(e, 1)) += F;
        f(nodeElement(e, 2)) += F;
    }

    // display force at nodes
    if (verbose) {

        std::cout<<"f\n";

        for (node i {0}; i < totalNodes; i++) {std::cout<<i<<" "<<f(i)<<"\n";}

        std::cout<<" --------------- "<<"\n";
    }

    // -- Dirichlet Boundary Conditions (BCs); u = f = 0

    // nodes on the boundary of the domain: top and bottom of the channel
    for (node n {0}; n < NodeX; n++) {

        if (verbose) {

            std::cout<<"Boundary Nodes: "<<n<<" "<<n + NodeX*(NodeY - 1)<<"\n";

        }

        // -- stiffness matrix K for BCs

        // bottom of channel
        K(n, n) *= 1e12;

        // top of channel
        K(n + NodeX*(NodeY - 1), n + NodeX*(NodeY - 1)) *= 1e12;

        // -- applied force/pressure gradient

        // bottom of channel
        f(n) = 0.0*K(n,n);

        // top of channel
        f(n + NodeX*(NodeY - 1)) = Vz*K(n + NodeX*(NodeY - 1), n + NodeX*(NodeY - 1));   
    }

    // nodes on the boundary of the domain: left and right of channel
    for (node n = 1; n < NodeY - 1; n++) {

        if (verbose) {
            
            std::cout<<"Boundary Nodes: "<<n*NodeX<<" "<<n*NodeX + NodeX - 1<<"\n";
        }

        // -- stiffness matrix K

        // left hand side of the channel
        K(n*NodeX, n*NodeX) *= 1e12;

        // right hand side of the channel
        K(n*NodeX + NodeX - 1, n*NodeX + NodeX - 1) *= 1e12;

        // -- applied force/pressure gradient

        // left hand side of the channel
        f(n*NodeX) = 0.0*K(n*NodeX, n*NodeX);

        // right hand side of the channel
        f(n*NodeX + NodeX - 1) = 0.0;
    }

    // -- find non-zero entries of the stiffness matrix K

    int num{0};

    // iterate through the rectangle of nodes to find non-zero entries
    for (node n {0}; n < totalNodes; n++) {
        for (node m {0}; m < totalNodes; m++) {

            // find the sum of non-zero entries until the whole rectangle has been covered
            if (K(n, m) != 0.0) num++;
        }
    }

    // print out non-zero entries of K
    std::cout<<"Non-Zero Elements in K: "<<--num<<"\n";

    // find the fluid velocity solution by inverting the stiffness matrix
    // Note: {u} = {f}[K]^-1
    // Note: inverse() method from Eigen
    // Note: for such a basic model, solving for {u} by direct inversion of [K] is fine.
    // However, in reality the order of [K] can be very large and it becomes impractical to use .inverse().
    // For such a situation, an iterative technique will be needed; e.g. Jacobi Method, Gauss-Seidal Method or Successive Over-Relaxation (SOR) Method.
    u = K.inverse()*f;

    // fluid velocity solution in uSolution.dat
    int nx {1};

    for (node n {0}; n < totalNodes; n++) {

        // fluid velocity at nodes
        uSol<<Node[n]<<" "<<u(n)<<"\n";

        nx++;

        if (nx > NodeX) {

            uSol<<" "<<"\n";
            nx = 1;
        }
    }

    // calculate the fluid flow rate; rate at which fluid moves through elements of mesh
    double ut {0.0};

    for (element e {0}; e < totalElements; e++) {

        ut += area(e)*(u(nodeElement(e, 0)) + u(nodeElement(e, 1)) + u(nodeElement(e, 2)));

    }

    ut /= 3.0;

    std::cout<<"Flow-Rate = "<<ut<<"\n";

    // -- make sure files are closed
   
    meshPts.close();
    eProps.close();
    uSol.close();

    // simple messages to the CLI to tell user where calculated data has been placed
    std::cout<<"Mesh Coordinates written to meshPoints.dat\n";
    std::cout<<"Element e, Center of Mass & Area of Elements written to file elementProperties.dat\n";
    std::cout<<"Node Number and Fluid Velocity at Nodes written to file uSolution.dat\n";

    // end of main
}
