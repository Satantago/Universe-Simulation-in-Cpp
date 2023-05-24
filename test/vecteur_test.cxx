#include <gtest/gtest.h>
#include "vecteur.hxx"

TEST(VecteurTest, VecteurCreation){
    Vecteur a = Vecteur();
    Vecteur b = Vecteur(1);
    Vecteur c = Vecteur(-6, 8);
    Vecteur d = Vecteur(7, -96, 78);

    EXPECT_EQ(a.x, 0);
    EXPECT_EQ(a.y, 0);
    EXPECT_EQ(a.z, 0);

    EXPECT_EQ(b.x, 1);
    EXPECT_EQ(b.y, 0);
    EXPECT_EQ(b.z, 0);

    EXPECT_EQ(c.x, -6);
    EXPECT_EQ(c.y, 8);
    EXPECT_EQ(c.z, 0);

    EXPECT_EQ(d.x, 7);
    EXPECT_EQ(d.y, -96);
    EXPECT_EQ(d.z, 78);
}

TEST(VecteurTest, BinaryOperation){
    Vecteur a = Vecteur(4, 8, 7);
    Vecteur b = Vecteur( -1, 5, -6);
    Vecteur c = a + b;
    EXPECT_EQ(c.x, 3);
    EXPECT_EQ(c.y, 13);
    EXPECT_EQ(c.z, 1);

    c = a - b;
    EXPECT_EQ(c.x, 5);
    EXPECT_EQ(c.y, 3);
    EXPECT_EQ(c.z, 13);

    c = 3 * a;
    EXPECT_EQ(c.x, 12);
    EXPECT_EQ(c.y, 24);
    EXPECT_EQ(c.z, 21);
}

TEST(VecteurTest, SingleOperation){
    Vecteur a = Vecteur(-3, 15, -23);
    Vecteur b = Vecteur( 6, -9, 45);

    a += b;
    EXPECT_EQ(a.x, 3);
    EXPECT_EQ(a.y, 6);
    EXPECT_EQ(a.z, 22);

    a -= 2 * b;
    EXPECT_EQ(a.x, -9);
    EXPECT_EQ(a.y, 24);
    EXPECT_EQ(a.z, -68);

    a *= -1;
    EXPECT_EQ(a.x, 9);
    EXPECT_EQ(a.y, -24);
    EXPECT_EQ(a.z, 68);
}

TEST(VecteurTest, Distance){
    Vecteur a = Vecteur();
    Vecteur b = Vecteur(1);

    EXPECT_EQ(a.distance(b), 1);

    a = Vecteur(1, 1, 0);

    EXPECT_EQ(a.distance(b), 1);

    a = Vecteur(0, 3, 0);
    b = Vecteur(0, 0, -4);

    EXPECT_EQ(a.distance(b), 5);
}