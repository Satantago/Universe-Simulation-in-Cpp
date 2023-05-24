#include <gtest/gtest.h>
#include <particule.hxx>
#include <univers.hxx>
#include <random>

TEST(UniversTest, UniversCreation){
    Univers u(3);

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    for(int i = 0; i < 1<<5; i++) {
        double x = distribution(generator);
        double y = distribution(generator);
        double z = distribution(generator);
        Vecteur v(x, y, z);
        Particule p(i);
        p.set_position(v);

        u.add_particule(p);
    }
}
