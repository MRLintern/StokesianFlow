#include <iostream>
#include <cmath>

class Vector {

private:

	// vector components, x & y
	double x{0.0}, y{0.0};

public:

	// -- constructors

	Vector() = default;
	Vector(double x0, double y0) : x(x0), y(y0) {}

	// -- accessor functions

	// x component of vector
	double getX() const { return x; }

	// y component of vector
	double getY() const { return y; }

	// -- setter functions

	// x component
	void setX(double x0) { x = x0; }

	// y component
	void setY(double y0) { y = y0; }

	// set components for vector
	void setVector(double x0, double y0) {

		x = x0;
		y = y0;
	}

	// Euclidean norm; size of vector
	double EuclidNorm() const {

		return sqrt(x*x + y*y);
	}

	// -- operator overloading

	// multiplication
	double operator*(const Vector& otherVector) const;

	// addition
	Vector operator+(const Vector& otherVector) const;

	// subtraction
	Vector operator-(const Vector& otherVector) const;

	// overloading insertion operator; for printing out Vector objects
	friend std::ostream& operator<<(std::ostream& output, const Vector& v0);

	// multiplication
	friend Vector operator*(const double& lhsVec, const Vector& rhsVec);
	friend Vector operator*(const Vector& lhsVec, const double& rhsVec);

	// division
	friend Vector operator/(const Vector& lhsVec, const double& rhsVec);

};

