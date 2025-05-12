// driver for viscousRiver

#include <iostream>
#include <fstream>
#include <array>

#include "vector.hpp"
#include "femSolver.hpp"


int main() {

    // allows data to be printed to the CLI; always choose 1!
    int verbose {0};

    
    // height : width of channel ratio: H/W ratio
    double HW {0.5};

    // Pressure gradient; the only force that acts on the fluid
    double Pz {-6};

    // initial velocity at the top of the channel
    double Vz {0.0};

    // -- file objects for output data; file extension: .csv
    // meshPts: mesh coordinates
    // eProps: element properties; node number of element etc
    // uSol: fluid velocity at nodes

    std::ofstream meshPts, eProps, uSol;

    // open files
    meshPts.open("meshPoints.csv");
    eProps.open("elementProperties.csv");
    uSol.open("uSolution.csv");

    std::cout<<"Enter 1 for verbose mode, 0 for not\n";
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

    // -- rectangle of nodes; i.e. rectangular mesh

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

            // set points for vector components
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
        
        // beta and gamma arrays
        // ---------------------

        // shape functions
        // ---------------
        // these interpolate the solution (e.g. fluid velocity) across the element
        // beta and gamma are components of the shape functions derivatives; i.e. how they change across the element
        // beta: difference in the y-coordinates of the nodes of the triangular element; the gradient of the shape function in the y-direction
        // gamma: difference in the x-coordinates of the nodes of the triangular element; the gradient of the shape function in the x-direction

        std::array<double, 3> beta{}; // i.e. beta[i] = y[j] - y[k]
        std::array<double, 3> gamma{}; // i.e. gamma[i] = x[k] - x[j]

        // find beta and gamma for K
        for (int i = 0; i < 3; i++) {

            int j {(i + 1) % 3};
            int k {(i + 2) % 3};

            // gradient of shape function in the y-direction
            beta[i] = Node[nodeElement(e, j)].getY() - Node[nodeElement(e, k)].getY();

            // gradient of shape function in the x-direction
            gamma[i] = Node[nodeElement(e, k)].getX() - Node[nodeElement(e, j)].getX();

            // display element, beta and gamma values at CLI
            if (verbose) {std::cout<<e<<" "<<i<<" "<<beta[i]<<" "<<gamma[i]<<"\n";}

        }

        for (int i {0}; i < 3; i++) {
            for (int j {0}; j < 3; j++) {

                // node positions for K
                node I {nodeElement(e, i)};
                node J {nodeElement(e, j)};

                // find K at (I, J)
                K(I, J) += (beta[i]*beta[j] + gamma[i]*gamma[j])/(4.0*area(e));
            }
        }
    }

    // find the determinant of the stiffness matrix K
    // Note: if | K | = 0; there are no unique solutions and u can't be found; i.e. we can't find the fluid velocity
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

    // -- Dirichlet Boundary Conditions (DBCs); u = f = 0

    // nodes on the boundary of the domain: top and bottom of the channel
    for (node n {0}; n < NodeX; n++) {

        if (verbose) {

            std::cout<<"Boundary Nodes: "<<n<<" "<<n + NodeX*(NodeY - 1)<<"\n";

        }

        // -- stiffness matrix K for DBCs

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

        // LHS of the channel
        K(n*NodeX, n*NodeX) *= 1e12;

        // RHS of the channel
        K(n*NodeX + NodeX - 1, n*NodeX + NodeX - 1) *= 1e12;

        // -- applied force/pressure gradient

        // LHS of the channel
        f(n*NodeX) = 0.0*K(n*NodeX, n*NodeX);

        // RHS of the channel
        f(n*NodeX + NodeX - 1) = 0.0;
    }

    // -- find non-zero entries of the stiffness matrix K

    int num{0};

    // iterate through the rectangle of nodes to find non-zero entries
    for (node n {0}; n < totalNodes; n++) {
        for (node m {0}; m < totalNodes; m++) {

            // find the sum of non-zero entries until the whole rectangle has been covered; i.e. until n = m
            if (K(n, m) != 0.0) { num++; }
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

    // find u
    u = K.inverse()*f;

    // -- fluid velocity solution in uSolution.csv

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

    // -- calculate the fluid flow rate; rate at which fluid moves through elements of mesh
    double vfr {0.0};

    for (element e {0}; e < totalElements; e++) {

        vfr += area(e)*(u(nodeElement(e, 0)) + u(nodeElement(e, 1)) + u(nodeElement(e, 2)));

    }

    // fluid-flow rate for vertices of triangular elements    
    vfr /= 3.0;

    std::cout<<"Flow-Rate = "<<vfr<<"\n";

    // -- make sure files are closed
   
    meshPts.close();
    eProps.close();
    uSol.close();

    // -- simple messages to the CLI to tell user where calculated data has been placed

    std::cout<<"Mesh Coordinates written to meshPoints.csv\n";
    std::cout<<"Element e, Center of Mass & Area of Elements written to file elementProperties.csv\n";
    std::cout<<"Node Number and Fluid Velocity at Nodes written to file uSolution.csv\n";

    
}
