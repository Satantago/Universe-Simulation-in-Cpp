#include <gtest/gtest.h>
#include <particule.hxx>

TEST(ParticuleTest, ParticuleCreation){
    Particule p = Particule(1);
}

TEST(PariculeTest, ParticuleSetters){
    Vecteur x = Vecteur(2,1,3);
    Vecteur v = Vecteur(1, 2,3);
    Particule p = Particule(1, x, v);
    Particule p1 = Particule(2, v, x);

    Vecteur pos = Vecteur(1, 1, 1);

    p.set_position(pos);

    p1.set_speed(pos);

    EXPECT_EQ(p.get_position().x, 1);
    EXPECT_EQ(p.get_position().y, 1);
    EXPECT_EQ(p.get_position().z, 1);

    EXPECT_EQ(p1.get_speed().x, 1);
    EXPECT_EQ(p1.get_speed().y, 1);
    EXPECT_EQ(p1.get_speed().z, 1);
}
