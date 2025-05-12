// implements the Vector class

#include <cmath>
#include "vector.hpp"

// dot product
double Vector::operator*(const Vector& otherVector) const {

	return x * otherVector.x + y * otherVector.y;
}

// addition
Vector Vector::operator+(const Vector& otherVector) const {

	return Vector(x + otherVector.x, y + otherVector.y);
}

// subtraction
Vector Vector::operator-(const Vector& otherVector) const {

	return Vector(x - otherVector.x, y - otherVector.y);
}

// output/insertion operator
std::ostream& operator<<(std::ostream& output, const Vector& v0) {

	output<<v0.getX()<<" "<<v0.getY()<<"\n";

	return output;
}

// scalar * vector multiplication
Vector operator*(const double& lhsVec, const Vector& rhsVec) {
    return Vector(lhsVec * rhsVec.getX(), lhsVec * rhsVec.getY());
}

// vector * scalar multiplication
Vector operator*(const Vector& lhsVec, const double& rhsVec) {
    return Vector(lhsVec.getX() * rhsVec, lhsVec.getY() * rhsVec);
}

// vector division by a scalar
Vector operator/(const Vector& lhsVec, const double& rhsVec) {

	return Vector(lhsVec.getX() / rhsVec, lhsVec.getY() / rhsVec);
}
