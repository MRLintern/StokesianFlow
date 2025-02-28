#include <iostream>
#include <cmath>
#include <fstream>
#include <random>
#include <vector>
#include <array>
#include <Eigen/Dense>

// Vector.hpp

class Vector {

    public:

        double x{};
        double y{};

        Vector() {x = 0.0; y = 0.0;};

        Vector (double x0, double y0) {

            x = x0; y = y0;
        };

        void setVector(double x0, double y0) {

            x = x0; y = y0;

        }

        double metric() {

            return sqrt(x*x + y*y);
        }

    double operator*(const Vector &otherVector);
    Vector operator+(const Vector &otherVector);
    Vector operator-(const Vector &otherVector);
    friend std::ostream &operator<<(std::ostream &output, Vector v0) {

        output<<v0.x<<" "<<v0.y;
        return output;
    }
};

// end of Vector.hpp

// ------------------------------------------------------------------------------------------------

// Vector.cpp

Vector Vector::operator-(const Vector &otherVector) {

    return Vector(x - otherVector.x, y - otherVector.y);
}

Vector Vector::operator+(const Vector &otherVector) {

    return Vector(x + otherVector.x, y + otherVector.y);
}

double Vector::operator*(const Vector &otherVector) {

    return x*otherVector.x + y*otherVector.y;
}

Vector operator*(const double &lhs, const Vector &rhs) {

    Vector v0 = Vector(lhs*rhs.x, lhs*rhs.y);

    return v0;
}

Vector operator*(const Vector &lhs, const double &rhs) {

    Vector v0 = Vector(rhs*lhs.x, rhs*lhs.y);

    return v0;
}

Vector operator/(const Vector &lhs, const double &rhs) {

    Vector v0 = Vector(lhs.x/rhs, lhs.y/rhs);

    return v0;
}

// end of Vector.cpp

// ------------------------------------------------------------------------------------------------

// Geometry.cpp

using node = int;
using element = int;
Vector *Node {nullptr};
using matrixElement = Eigen::Matrix<int, Eigen::Dynamic, 3>;
matrixElement nodeElement;

double area(element e) {

    double a1 {Node[nodeElement(e, 1)].x - Node[nodeElement(e, 0)].x};
    double a2 {Node[nodeElement(e, 0)].x - Node[nodeElement(e, 2)].x};
    double b1 {Node[nodeElement(e, 0)].y - Node[nodeElement(e, 1)].y};
    double b2 {Node[nodeElement(e, 2)].y - Node[nodeElement(e, 0)].y};

    return 0.5*(a1*b2 - a2*b1);
}

Vector com(element e) {

    return (Node[nodeElement(e, 0)] + Node[nodeElement(e, 1)] + Node[nodeElement(e, 2)])/3.0;
}

// end of Geometry.cpp

// ----------------------------------------------------------------

// main

int main() {

    int verbose {0};
    int NodeX {6}, NodeY {6};
    double HW {0.5};
    double Pz {-6};
    double Vz {0.0};
    std::ofstream meshPts, eProps, uSol;
    meshPts.open("meshPoints.dat");
    eProps.open("elementProperties.dat");
    uSol.open("uSolution.dat");

    std::cout<<"Enter 1 for verbose mode, 0 if not\n";
    std::cin>>verbose;

    int totalNodes {NodeX*NodeY};
    int totalElements {(NodeX - 1)*(NodeY - 1)*2};

    std::cout<<"Total Number of Nodes: "<<totalNodes<<"\n";
    std::cout<<"Total Number of Elements: "<<totalElements<<"\n";

    Eigen::VectorXd f(totalNodes), u(totalNodes);
    Eigen::MatrixXd K(totalNodes, totalNodes);
    nodeElement = matrixElement(totalElements, 3);

    Node = new Vector[totalNodes];

    node nodeNumber {0};

    for (int iy {0}; iy < NodeY; iy++) {
        for (int ix {0}; ix < NodeX; ix++) {

            
            double x = (double) ix/(NodeX - 1);
           
            //# double x = static_cast<double>(ix/(NodeX - 1));
            double y = (double) iy/(NodeY - 1);
            
            //# double y = static_cast<double>(iy/(NodeY - 1));
            
            y *= HW;
            Node[nodeNumber].setVector(x, y);
            nodeNumber++;
        }
    }

    if (verbose) {std::cout<<"Element and Associated Nodes: "<<"\n";}

    for (int iy = 0; iy < NodeY - 1; iy++) {
        for (int ix = 0; ix < NodeX - 1; ix++) {

            node i = iy*(NodeX - 1) + ix;
            element eNumber = 2*i;
            node j = ix + NodeX*iy;

            nodeElement(eNumber, 0) = j;
            nodeElement(eNumber, 1) = j + NodeX + 1;
            nodeElement(eNumber, 2) = j + NodeX;

            if (verbose) {std::cout<<eNumber<<" "<<nodeElement(eNumber, 0)<<" "<<nodeElement(eNumber, 1)<<" "<<nodeElement(eNumber, 2)<<"\n";}

            eNumber++;

            nodeElement(eNumber, 0) = j;
            nodeElement(eNumber, 1) = j + 1;
            nodeElement(eNumber, 2) = j + 1 + NodeX;

            if (verbose) {std::cout<<eNumber<<" "<<nodeElement(eNumber, 0)<<" "<<nodeElement(eNumber, 1)<<" "<<nodeElement(eNumber, 2)<<"\n";}

        }
    }

    if (verbose) {
        std::cout<<" ----------------- "<<"\n";

        for (element e = 0; e < totalElements; e++) {

            meshPts<<Node[nodeElement(e, 0)]<<" "<<Node[nodeElement(e, 1)]<<"\n";
            meshPts<<Node[nodeElement(e, 1)]<<" "<<Node[nodeElement(e, 2)]<<"\n";
            meshPts<<Node[nodeElement(e, 2)]<<" "<<Node[nodeElement(e, 0)]<<"\n";
            eProps<<e<<" "<<com(e)<<" "<<area(e)<<"\n";
        }
    }

    // stiffness matrix K
    for (element e = 0; e < totalElements; e++) {
        double beta[3], gamma[3];
        // std::array<double, 3> beta{}
        // std::array<double, 3> gamma{};
        for (int i = 0; i < 3; i++) {

            int j = (i + 1) % 3;
            int k = (i + 2) % 3;

            beta[i] = Node[nodeElement(e, j)].y - Node[nodeElement(e, k)].y;
            gamma[i] = Node[nodeElement(e, k)].x - Node[nodeElement(e, j)].x;

            if (verbose) {std::cout<<e<<" "<<i<<" "<<beta[i]<<" "<<gamma[i]<<"\n";}

        }

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {

                node I = nodeElement(e, i);
                node J = nodeElement(e, j);
                K(I, J) += (beta[i]*beta[j] + gamma[i]*gamma[j])/(4.0*area(e));
            }
        }
    }

    if (verbose) {

        std::cout<<" --------------- "<<"\n";
        std::cout<<"K "<<"\n";
        std::cout<<"det(K) = "<<K.determinant()<<"\n";

        for (node i = 0; i < totalNodes; i++) {
            for (node j = 0; j < totalNodes; j++) {

                std::cout<<i<<" "<<j<<" "<<K(i, j)<<"\n";
            }
        }

        std::cout<<" --------------- "<<"\n";
    }

    // force vector
    for (element e = 0; e < totalElements; e++) {

        double F = -Pz*area(e)/3.0;

        f(nodeElement(e, 0)) += F;
        f(nodeElement(e, 1)) += F;
        f(nodeElement(e, 2)) += F;
    }

    if (verbose) {

        std::cout<<"f\n";

        for (node i = 0; i < totalNodes; i++) {std::cout<<i<<" "<<f(i)<<"\n";}

        std::cout<<" --------------- "<<"\n";
    }

    // boundary conditions
    for (node n = 0; n < NodeX; n++) {

        if (verbose) {

            std::cout<<"Boundary Nodes: "<<n<<" "<<n + NodeX*(NodeY - 1)<<"\n";

        }

        K(n, n) *= 1e12;
        K(n + NodeX*(NodeY - 1), n + NodeX*(NodeY - 1)) *= 1e12;
        f(n) = 0.0*K(n,n);
        f(n + NodeX*(NodeY - 1)) = Vz*K(n + NodeX*(NodeY - 1), n + NodeX*(NodeY - 1));   
    }

    for (node n = 1; n < NodeY - 1; n++) {

        if (verbose) {
            
            std::cout<<"Boundary Nodes: "<<n*NodeX<<" "<<n*NodeX + NodeX - 1<<"\n";
        }

        K(n*NodeX, n*NodeX) *= 1e12;
        K(n*NodeX + NodeX - 1, n*NodeX + NodeX - 1) *= 1e12;
        f(n*NodeX) = 0.0*K(n*NodeX, n*NodeX);
        f(n*NodeX + NodeX - 1) = 0.0;
    }

    int num{0};

    for (node n = 0; n < totalNodes; n++) {
        for (node m = 0; m < totalNodes; m++) {

            if (K(n, m) != 0.0) num++;
        }
    }

    std::cout<<"Non-Zero Elements in K: "<<--num<<"\n";

    // solution
    u = K.inverse()*f;

    // store solution
    int nx = 1;

    for (node n = 0; n < totalNodes; n++) {

        uSol<<Node[n]<<" "<<u(n)<<"\n";

        nx++;

        if (nx > NodeX) {

            uSol<<" "<<"\n";
            nx = 1;
        }
    }

    // obtain flow rate
    double ut = 0.0;

    for (element e = 0; e < totalElements; e++) {

        ut += area(e)*(u(nodeElement(e, 0)) + u(nodeElement(e, 1)) + u(nodeElement(e, 2)));

    }

    ut /= 3.0;

    std::cout<<"Flow-Rate = "<<ut<<"\n";

    delete[] Node;

    meshPts.close();
    eProps.close();
    uSol.close();

    std::cout<<"Mesh Coordinates written to meshPoints.dat\n";
    std::cout<<"Element e, Center of Mass & Area of Elements written to file elementProperties.dat\n";
    std::cout<<"Node Number and Fluid Velocity at Nodes written to file uSolution.dat\n";

    return 0;

}
