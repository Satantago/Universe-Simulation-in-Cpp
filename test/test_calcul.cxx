#include <iostream>
#include "univers.hxx"
#include <random>
#include <chrono>

//test de perfomance en insertion
int main(){
    Univers u(3);

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    for(int k = 3; k < 8; k++) {
        for (int i = 0; i < 1 << k; i++) {
            double x = distribution(generator);
            double y = distribution(generator);
            double z = distribution(generator);
            Vecteur v(x, y, z);
            Particule p(i);
            p.set_position(v);

            u.add_particule(p);
        }
        auto start = std::chrono::steady_clock::now();
        u.next();
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_sec = end - start;
        std::cout << "k = " << k << " elapsed time: " << elapsed_sec.count() << "s\n";
    }
    return EXIT_SUCCESS;
}