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
    double eta_t = 0.00005;
    double epsilon = 5;
    double sigma = 1;
    double mass = 1;
    double rcut = 2.5*sigma;
    double g = 12;
    Vecteur ld = Vecteur{250, 40, rcut};
    unordered_map<int, Particule> particules;

    double dist = 1.12246204831*sigma;

    Vecteur v_rectangle {0, 0, 0};

    int id = 1;

    double x_base = 0;
    double y_base = 0;

    for (int i = 0; i < 160; i++) {
        for (int j = 0; j < 40; j++) {
            Vecteur pos = Vecteur(x_base + dist * i, y_base + dist * j, 0);
            Particule p = Particule(id, pos, v_rectangle, mass);
            particules.insert({id, p});
            id++;
        }
    }


    Vecteur v_carre {0, 10, 0};
    x_base = 60*dist;
    y_base = 45*dist;
    for (int i = 0; i < 40; i++) {
        for (int j = 0; j < 40; j++) {
            Vecteur pos = Vecteur{x_base + dist * i, y_base + dist * j, 0};
            Particule p = Particule(id, pos, v_carre, mass);
            particules.insert({id, p});
            id++;
        }
    }
    
    Univers u = Univers(dim, t, eta_t, sigma, epsilon, rcut, ld, particules, g);

    double tmax = t;
    u.set_t(0);
    while (u.get_t() < tmax)
        u.next_grid();
}