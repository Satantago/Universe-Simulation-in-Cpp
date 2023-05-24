#include "vecteur.hxx"
#include "particule.hxx"
#include "univers.hxx"
#include <iostream>

int main(){
    Vecteur v1 = Vecteur(12);
    Vecteur v2 = Vecteur(45, 69);
    Vecteur v3 = Vecteur(48, 78, 32);

    Vecteur b = Vecteur();

    b = v1 + b + v3 * 2;
    b = v2 * 4 + b;

    std::cout << b.getX() << ", " << b.getY() << ", " << b.getZ() << std::endl;

    return EXIT_SUCCESS;
}