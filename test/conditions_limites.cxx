#include <iostream>
#include "univers.hxx"
#include "particule.hxx"
#include "vecteur.hxx"
#include "grid.hxx"
#include <unordered_map>

using namespace std;

int main() {
    int dim = 2;
    double t = 19.5;
    double eta_t = 0.5;
    double epsilon = 5;
    double sigma = 1;
    double mass = 1;
    double rcut = 2.5*sigma;
    double g = 12;
    Vecteur ld = Vecteur{250, 40, rcut};
    unordered_map<int, Particule> particules;

    double dist = 1.12246204831*sigma;

    Vecteur x = Vecteur{15, 0, 0};
    Vecteur v = Vecteur{0, -50, 0};
    Particule p1 = Particule(1, x, v, mass);

    x = Vecteur{100, 15, 0};
    v = Vecteur{50, 0, 0};
    Particule p2 = Particule(2, x, v, mass);

    x = Vecteur{0, 100, 0};
    v = Vecteur{0, 50, 0};
    Particule p3 = Particule(3, x, v, mass);

    x = Vecteur{85, 100, 0};
    v = Vecteur{0, 50, 0};
    Particule p4 = Particule(4, x, v, mass);
        
    particules.insert({1, p1});
    particules.insert({2, p2});
    particules.insert({3, p3});
    particules.insert({4, p4});

    Univers u = Univers(dim, t, eta_t, sigma, epsilon, rcut, ld, particules, g);

    double tmax = t;
    u.set_t(0);
    while (u.get_t() < tmax)
    {   // reflection
        // u.next_cond_limites(1); 

        // periodic
        // u.next_cond_limites(2);

        // absorption
        // u.next_cond_limites(3);

        // potentiel reflexion
        u.next_potentiel_reflexion();
    }
}