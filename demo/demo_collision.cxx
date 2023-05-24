#include <iostream>
#include "univers.hxx"
#include "particule.hxx"
#include "vecteur.hxx"
#include "grid.hxx"
#include <unordered_map>
#include <random>

using namespace std;

int main() {
    int dim = 2;
    double t = 29.5;
    double eta_t = 0.00005;
    double epsilon = 1;
    double sigma = 1;
    double mass = 1;
    double rcut = 2.5*sigma;
    double g = -12;
    Vecteur ld = Vecteur{250, 180, rcut};
    double Ecd = 0.005;

    unordered_map<int, Particule> particules;

    double dist = 1.12246204831*sigma;

    Vecteur v_eau {0, 0, 0};

    int id = 1;

    double x_base = 0;
    double y_base = 0;

    for (int i = 0; i < 749; i++) {
        for (int j = 0; j < 23; j++) {
            Vecteur pos = Vecteur(x_base + dist * i, y_base + dist * j, 0);
            Particule p = Particule(id, pos, v_eau, mass);
            particules.insert({id, p});
            id++;
        }
    }


    Vecteur v_cercle {0, 10, 0};
    x_base = 375*dist;
    y_base = 40*dist;
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1.0, 1.0);

    double x;
    double y;
    for (int i = 0; i < 395; i++) {
        x = 1; y = 1;
        do {
            x = dis(gen);
            y = dis(gen);
        } while (x*x + y*y > 1);
        Vecteur pos = Vecteur{x_base + x*50*dist, y_base + y*10*dist, 0};
        Particule p = Particule(id, pos, v_cercle, mass);
        particules.insert({id, p});
        id++;
    }
    
    Univers u = Univers(dim, t, eta_t, sigma, epsilon, rcut, ld, particules, g);

    double tmax = t;
    int n = 0;
    u.set_t(0);
    while (u.get_t() < tmax) {
        
        n++;
        if (n % 1000 == 0) {
            double Ec = 0;
            for (auto &e : u.get_particules()) {
                Particule &p = e.second;
                Ec += 0.5*p.get_mass()*(p.get_speed().norme());        
            }
            double beta = sqrt(Ecd/Ec);
            for (auto &e : u.get_particules()) {
                Particule &p = e.second;
                Vecteur new_speed = p.get_speed()*beta;
                p.set_speed(new_speed);
            }
        }

        u.next_potentiel_reflexion();
    }
}