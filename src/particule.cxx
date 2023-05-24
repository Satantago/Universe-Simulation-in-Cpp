/**
 * @file particule.cxx
 * @brief Defines the Particule class for representing a particle.
 * @date May 2023
 * @authors Achraf Kerzazi & Mohamed Errazki
 */

#include "particule.hxx"
#include <math.h>
#include <vector>


double Particule::distance(const Particule &p) {
    return sqrt(pow(x.getX() - p.get_position().getX(), 2) + pow(x.getY() - p.get_position().getY(), 2) + pow(x.getZ() - p.get_position().getZ(), 2));
}

int Particule::get_id() {return id;}

Vecteur Particule::get_speed() {return v;}

void Particule::set_speed(Vecteur &vi) {v = vi;}

Vecteur Particule::get_position() const {return x;}

void Particule::set_position(Vecteur &xi) {x = xi;}

Particule::Particule(int _id, double _mass) {
    id = _id;
    mass = _mass;
}

double Particule::get_mass() {
    return mass;
}


/*
-------------------------------------------------
-------------------------------------------------
LAB 4
-------------------------------------------------
-------------------------------------------------
*/

double Particule::lennardJonesPotential(const Particule &p, double epsilon, double sigma) {
    double r = distance(p);
    return 4 * epsilon * (pow(sigma/r, 12) - pow(sigma/r, 6));
}

Vecteur Particule::computeLennardForce(const Particule &p, double epsilon, double sigma) {
    double r = distance(p);
    double f = 24 * epsilon * pow(1/r, 2) * pow(sigma/r, 6) * (1 - 2*pow(sigma/r, 6));
    return (p.x - x) * f;
}

Vecteur Particule::computeLennardAll(std::vector<Particule> particules, double epsilon, double sigma) {
    Vecteur force = Vecteur(0, 0, 0);
    for (Particule &p : particules){
        if(p.get_id() != id){
            force = force + computeLennardForce(p, epsilon, sigma);
        }
    }
    return force;
}

Particule::Particule(int _id, Vecteur &_x, Vecteur &_v, double _mass) {
    id = _id;
    x = _x;
    v = _v;
    mass = _mass;
}