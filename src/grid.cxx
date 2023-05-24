/**
 * @file grid.cxx
 * @brief Defines the Grid class for representing the grid used for spatial partitioning
 *        and the Cell class for representing a cell in the grid.
 * 
 * @date May 2023
 * @authors Achraf Kerzazi & Mohamed Errazki
 */


#include <unordered_map>
#include <vector>
#include "grid.hxx"
#include "particule.hxx"
#include "vecteur.hxx"
#include <iostream>
#include <cmath>

Vecteur Cell::getBottomLeft() {return bottom_left;}
Vecteur Cell::getCellSize() {return cell_size;}
Vecteur Cell::getIndices() {return Vecteur(i, j, k);}

Cell::Cell(int i, int j, int k, Vecteur cell_size, Vecteur bottom_left) {
    this->i = i;
    this->j = j;
    this->k = k;
    this->cell_size = cell_size;
    this->bottom_left = bottom_left;
    this->particles = std::unordered_map<int, Particule>();
}


bool Cell::contains(Particule & p) {
    if (this->particles.find(p.get_id()) != this->particles.end())
        return true;
    return false;
}

void Cell::insert(Particule &p) {
    particles.insert({p.get_id(), p});
}
    
void Cell::remove(Particule &p) {
    particles.erase(p.get_id());
}

Vecteur Cell::getCenter() {
    return bottom_left + 0.5*cell_size;
}

Cell::Cell() {
    this->i = 0;
    this->j = 0;
    this->k = 0;
    this->cell_size = Vecteur(0,0,0);
    this->bottom_left = Vecteur(0,0,0);
    this->particles = std::unordered_map<int, Particule>();
}

Grid::Grid(Vecteur ld, double rcut, std::unordered_map<int, Particule> particleList) {
    double xmin = std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::min();
    double ymin = std::numeric_limits<double>::max();
    double ymax = std::numeric_limits<double>::min();
    double zmin = std::numeric_limits<double>::max();
    double zmax = std::numeric_limits<double>::min();


    // Find xmin, xmax, ymin, ymax, zmin and zmax
    for (auto& pair : particleList) {
        Particule &particle = pair.second;
        if (particle.get_position().getX() < xmin) xmin = particle.get_position().getX();
        if (particle.get_position().getX() > xmax) xmax = particle.get_position().getX();
        if (particle.get_position().getY() < ymin) ymin = particle.get_position().getY();
        if (particle.get_position().getY() > ymax) ymax = particle.get_position().getY();
        if (particle.get_position().getZ() < zmin) zmin = particle.get_position().getZ();
        if (particle.get_position().getZ() > zmax) zmax = particle.get_position().getZ();
    }
    
    xmax += rcut;
    xmin -= rcut;
    ymax += rcut;
    ymin -= rcut;
    zmax += rcut;
    zmin -= rcut;

    this->ld = ld;
    this->rcut = rcut;

    int nx = floor(ld.getX() / rcut);
    int ny = floor(ld.getY() / rcut);
    int nz = floor(ld.getZ() / rcut);

    double cell_size_x = (xmax - xmin) / nx;
    double cell_size_y = (ymax - ymin) / ny;
    double cell_size_z = (zmax - zmin ) / nz;
    Vecteur cell_size = Vecteur(cell_size_x, cell_size_y, cell_size_z);

    std::vector<std::vector<std::vector<Cell>>> grid(nx, std::vector<std::vector<Cell>>(ny, std::vector<Cell>(nz)));
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                Vecteur bottom_left = Vecteur(xmin + i * cell_size_x, ymin + j * cell_size_y, zmin + k * cell_size_z);
                grid[i][j][k] = Cell(i, j, k, cell_size, bottom_left);
            }
        }
    }

    for (auto &pair : particleList) {
        Particule &p = pair.second;
        int i = floor((p.get_position().getX() - xmin) / cell_size_x);
        int j = floor((p.get_position().getY() - ymin) / cell_size_y);
        int k = floor((p.get_position().getZ() - zmin) / cell_size_z);
        
        int dec_x = 0;
        int dec_y = 0;
        int dec_z = 0;
        if (i == nx) dec_x = -1;
        if (j == ny) dec_y = -1;
        if (k == nz) dec_z = -1;

        grid[i+dec_x][j+dec_y][k+dec_z].insert(p);
    }
    this->grid = grid;
}

Grid::Grid() {
    this->ld = Vecteur(10, 10, 10);
    this->rcut = 1;
    this->grid = std::vector<std::vector<std::vector<Cell>>>();
}


Vecteur Grid::getCellSize() {
    return this->grid[0][0][0].getCellSize();
}

Cell& Grid::containingCell(Particule &p) {
    Vecteur mini = this->grid[0][0][0].getBottomLeft();
    double xmin = mini.getX();
    double ymin = mini.getY();
    double zmin = mini.getZ();

    Vecteur cell_size = this->grid[0][0][0].getCellSize();
    double cell_size_x = cell_size.getX();
    double cell_size_y = cell_size.getY();
    double cell_size_z = cell_size.getZ();

    Vecteur pos = p.get_position();
    int i = floor((pos.getX() - xmin) / cell_size_x);
    int j = floor((pos.getY() - ymin) / cell_size_y);
    int k = floor((pos.getZ() - zmin) / cell_size_z);

    int nx = this->grid.size();
    int ny = this->grid[0].size();
    int nz = this->grid[0][0].size();

    // if it evers goe out of grid, we store in in the last cell
    if (i >= nx) i = nx-1;
    if (j >= ny) j = ny-1;
    if (k >= nz) k = nz-1;

    if (i < 0) i = 0;
    if (j < 0) j = 0;
    if (k < 0) k = 0;

    return this->grid[i][j][k];
}

void Grid::updateCell(Particule &p1, Particule &p2) {
    Cell &c1 = this->containingCell(p1);
    Cell &c2 = this->containingCell(p2);

    if (c1.getIndices().getX() == c2.getIndices().getX() 
        && c1.getIndices().getY() == c2.getIndices().getY()
        && c1.getIndices().getZ() == c2.getIndices().getZ()) {
        return;
    }

    c1.remove(p1);
    c2.insert(p2);
}

std::vector<Cell> Grid::getNeighbors(Cell &c) {
    std::vector<Cell> neighbors;
    int i = c.getIndices().getX();
    int j = c.getIndices().getY();
    int k = c.getIndices().getZ();

    int nx = this->grid.size();
    int ny = this->grid[0].size();
    int nz = this->grid[0][0].size();

    for (int x = -1; x <= 1; x++) {
        if (i + x < 0 || i + x >= nx) continue;
        for (int y = -1; y <= 1; y++) {
            if (j + y < 0 || j + y >= ny) continue;
            for (int z = -1; z <= 1; z++) {
                if (k + z < 0 || k + z >= nz) continue;
                neighbors.push_back(this->grid[i + x][j + y][k + z]);
            }
        }
    }
    return neighbors;
}