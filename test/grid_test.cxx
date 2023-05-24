#include <gtest/gtest.h>
#include <grid.hxx>
#include <particule.hxx>

TEST(CellTest, CellCreation){
    Cell cl = Cell(0, 0, 0, Vecteur(), Vecteur());
}

TEST(CellTest, CellOps){
    Cell cl = Cell();
    Particule p = Particule(1);
    Particule b = Particule(2);

    cl.insert(p);
    EXPECT_EQ(cl.contains(p), true);
    EXPECT_EQ(cl.contains(b), false);

    cl.insert(b);
    cl.remove(p);
    EXPECT_EQ(cl.contains(p), false);
}

TEST(GridTest, GridCreation){
    Grid gl = Grid(Vecteur(), 1, std::vector<Particule>{});
}
