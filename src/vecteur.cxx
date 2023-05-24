/**
 * @file vecteur.cxx
 * @brief Defines the Vecteur class for representing a 3D vector.
 * @date May 2023
 * @authors Achraf Kerzazi & Mohamed Errazki
 */

#include "vecteur.hxx"
#include <iostream>
#include <cmath>

Vecteur::Vecteur(double x, double y, double z) {
    this->x = x; 
    this->y = y; 
    this->z = z;
}

Vecteur & Vecteur::operator+=(const Vecteur &v) {
    x += v.getX(); y += v.getY(); z += v.getZ();
    return *this;
}

Vecteur & Vecteur::operator-=(const Vecteur &v) {
    x -= v.getX(); y -= v.getY(); z -= v.getZ();
    return *this;
}

Vecteur & Vecteur::operator*=(double p) {
    x *= p; y *= p; z *= p;
    return *this;
}

Vecteur Vecteur::operator+(const Vecteur &v) {
    Vecteur v1(x + v.getX(), y + v.getY(), z + v.getZ());
    return v1;
}

Vecteur Vecteur::operator-(const Vecteur &v) {
    Vecteur v1(x - v.getX(), y - v.getY(), z - v.getZ());
    return v1;
}

Vecteur Vecteur::operator*(double p) {
    Vecteur v1(p * x, p * y, p * z);
    return v1;
}

Vecteur & Vecteur::operator=(const Vecteur &v) {
    x = v.getX(); y = v.getY(); z = v.getZ();
    return *this;
}

double &Vecteur::operator()(int i) {
    if (i == 0) return x;
    if (i == 1) return y;
    if (i == 2) return z;
    throw std::out_of_range("Vecteur::operator() : index out of range");
}

double Vecteur::operator()(int i) const {
    if (i == 0) return x;
    if (i == 1) return y;
    if (i == 2) return z;
    throw std::out_of_range("Vecteur::operator() : index out of range");
}

Vecteur operator*(double p, const Vecteur &v) {
    Vecteur v1(p * v.getX(), p * v.getY(), p * v.getZ());
    return v1;
}

Vecteur operator+(const Vecteur &v1, const Vecteur &v2) {
    Vecteur v3(v1.getX() + v2.getX(), v1.getY() + v2.getY(), v1.getZ() + v2.getZ());
    return v3;
}

Vecteur operator-(const Vecteur &v1, const Vecteur &v2) {
    Vecteur v3(v1.getX() - v2.getX(), v1.getY() - v2.getY(), v1.getZ() - v2.getZ());
    return v3;
}

std::ostream & operator<<(std::ostream &os, const Vecteur &v) {
    os << "(" << v.getX() << ", " << v.getY() << ", " << v.getZ() << ")";
    return os;
}

double Vecteur::distance(const Vecteur &v) {
    return sqrt(pow(x - v.getX(), 2) + pow(y - v.getY(), 2) + pow(z - v.getZ(), 2));
}

double Vecteur::norme() {
    return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
}

double Vecteur::getX() const {return x;}
double Vecteur::getY() const {return y;}
double Vecteur::getZ() const {return z;}

void Vecteur::setX(double x) {this->x = x;}
void Vecteur::setY(double y) {this->y = y;}
void Vecteur::setZ(double z) {this->z = z;}