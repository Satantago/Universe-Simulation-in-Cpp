#include "particule.hxx"
#include "vecteur.hxx"
#include <unordered_map>
#include "univers.hxx"

int main() {
    int dim = 2;
    double t = 468.5;
    double eta_t = 0.015;
    double epsilon = 5;
    double sigma = 1;
    double rcut = 2.5*sigma;
    double g = 12;
    Vecteur ld = Vecteur{250, 40, rcut};
    std::unordered_map<int, Particule> particules;
    double dist = 1.12246204831*sigma;

    Vecteur pos_soleil = Vecteur{0, 0, 0};
    Vecteur vit_soleil = Vecteur{0, 0, 0};
    Particule soleil = Particule(1, pos_soleil, vit_soleil, 1);
    particules.insert({1, soleil});

    Vecteur pos_terre = Vecteur{0, 1, 0};
    Vecteur vit_terre = Vecteur{-1, 0, 0};
    Particule terre = Particule(2, pos_terre, vit_terre, 3.0e-6);
    particules.insert({2, terre});

    Vecteur pos_jupiter = Vecteur{0, 5.36, 0};
    Vecteur vit_jupiter = Vecteur{-0.425, 0, 0};
    Particule jupiter = Particule(3, pos_jupiter, vit_jupiter, 9.55e-4);
    particules.insert({3, jupiter});

    Vecteur pos_haley = Vecteur{34.75, 0, 0};
    Vecteur vit_haley = Vecteur{0, 0.0296, 0};
    Particule haley = Particule(4, pos_haley, vit_haley, 1.0e-14);
    particules.insert({4, haley});

    Univers u = Univers(dim, t, eta_t, sigma, epsilon, rcut, ld, particules, g);
    double tmax = t;
    u.set_t(0);
    while (u.get_t() < tmax)
    {   
        u.next();
        u.show_state();
    }

}