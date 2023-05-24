#include <set>
#include <list>
#include <deque>
#include <vector>
#include <unordered_map>
#include "particule.hxx"
#include "vecteur.hxx"

using namespace std;

#include <random>
#include <iostream>
#include <chrono>


int main() {
    
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> pos(0.0, 100.0);
    uniform_real_distribution<double> vel(0.0, 10.0);
    uniform_real_distribution<double> mass(0.0, 5.0);
    uniform_real_distribution<double> force(0.0, 2.5);

    vector<Particule> particleList;

    auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < 100000; i++) {
        Vecteur position = Vecteur(pos(mt), pos(mt), pos(mt));
        Vecteur vitesse = Vecteur(vel(mt), vel(mt), vel(mt));
        Particule p = Particule(i+1, position, vitesse, mass(mt));
        particleList.push_back(p);
    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
}