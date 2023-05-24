/**
 * @file univers.cxx
 * @brief Implements the Univers class for representing the universe in which particles interact.
 *        Universe creation, grid configuration, advancing system in time, outputting system state
 *        to file, handling boundary conditions, and visualizing system state.
 * @date May 2023
 * @authors Achraf Kerzazi & Mohamed Errazki
 */


#include "univers.hxx"
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include "particule.hxx"
#include "vecteur.hxx"
#include "grid.hxx"
#include <array>


void Univers::set_t(double ti) {t = ti;}
void Univers::set_eta(double eta) {eta_t = eta;}
void Univers::set_speed(int index, Vecteur &v) {particules.at(index).set_speed(v);}
std::unordered_map<int, Particule> Univers::get_particules() {return particules;}
double Univers::get_t() {return t;}
double Univers::get_eta() {return eta_t;}
Vecteur Univers::get_position(int index) {return particules.at(index).get_position();}
Vecteur Univers::get_speed(int index) {return particules.at(index).get_speed();}
double Univers::get_mass(int index) {return particules.at(index).get_mass();}
int Univers::get_nb() {return particules.size();}

std::unordered_map<int, Vecteur> Univers::force_field() {
    std::unordered_map<int, Vecteur> forces;

    for(auto &pair : particules) {
        Particule &e = pair.second;
        Vecteur force = Vecteur();
        for (auto &pair2 : particules){
            Particule &p = pair2.second;
            if(p.get_id() != e.get_id()){
                double rij = e.distance(p);
                double f = (p.get_mass() * e.get_mass()) / pow(rij, 3);
                force = force + (p.get_position() + e.get_position() * (-1)) * f;
            }
        }
        forces.insert({e.get_id(), force});
    }
    return forces;
}

void Univers::next(){

    t += eta_t;
    std::unordered_map<int, Vecteur> forces_old = force_field();
    for(auto &pair : particules){
        Particule &e = pair.second;
        Vecteur pos = e.get_position() + (e.get_speed() + forces_old.at(e.get_id()) * (0.5 / e.get_mass()) * eta_t) * eta_t;
        e.set_position(pos);
    }

    std::unordered_map<int, Vecteur> forces = force_field();

    for(auto &pair : particules){
        Particule &e = pair.second;
        Vecteur spe = e.get_speed() + (forces.at(e.get_id()) + forces_old.at(e.get_id())) * eta_t *  (0.5 / e.get_mass());
        e.set_speed(spe);
    }
}

Univers::Univers(int d, double time, double eta) {
    dim = d;
    t = time;
    eta_t = eta;
}

void Univers::add_particule(Particule &p) {
    particules.insert({p.get_id() ,p});
}

void Univers::show_state() {
    std::cout << "t = " << t << std::endl;

    for(auto &pair : particules){
        Particule &e = pair.second;
        std::cout << "Particule id " << e.get_id() << std::endl;
        std::cout << " Position:" << e.get_position().getX() << ", " <<  e.get_position().getY() << ", " << e.get_position().getZ() << std::endl;
        std::cout << " Vitesse:"  << e.get_speed().getX()    << ", " <<  e.get_speed().getY()    << ", " << e.get_speed().getZ()    << std::endl;
        std::cout << "\n";
    }
}

/*
----------------------------------------------------
----------------------------------------------------
----------------------------------------------------
LAB 4
----------------------------------------------------
----------------------------------------------------
----------------------------------------------------
*/

Univers::Univers(int dim, double t, double eta_t, double sigma, double epsilon, double rcut, Vecteur ld, std::unordered_map<int, Particule> particleList, double g) {
    this->sigma = sigma;
    this->epsilon = epsilon;
    this->rcut = rcut;
    this->ld = ld;
    this->particules = particleList;
    this->grid = Grid(ld, rcut, particleList);
    this->dim = dim;
    this->t = t;
    this->eta_t = eta_t;
    this->g = g;
}


std::unordered_map<int, Vecteur> Univers::force_field_grid() {
    std::unordered_map<int, Vecteur> forces;

    for(auto &pair : particules) {
        Particule &e = pair.second;
        Vecteur force = Vecteur();
        Cell cell = grid.containingCell(e);
        std::vector<Cell> neighbors = grid.getNeighbors(cell);
        for (Cell &c: neighbors) {
            Vecteur cell_center = c.getCenter();
            if (cell_center.distance(e.get_position()) < rcut) {
                for (auto &p : c.particles) {
                    if (p.first != e.get_id()) {
                        force += e.computeLennardForce(p.second, sigma, epsilon);
                    }
                }
            }
        }
        
        forces.insert({e.get_id(), force});
    }
    return forces;
}

void Univers::next_grid() {
    t += eta_t;
    std::unordered_map<int, Vecteur> forces_old_leonnard = force_field_grid();

    for(auto &pair : particules){
        Particule &e = pair.second;
        Vecteur pos = e.get_position() + (e.get_speed() + forces_old_leonnard.at(e.get_id()) * (0.5 / e.get_mass()) * eta_t) * eta_t;
        Particule p = e;
        p.set_position(pos);
        this->grid.updateCell(e,p);
        e = p;
    }

    std::unordered_map<int, Vecteur> forces_new_leonnard = force_field_grid();
    for (auto &pair: particules) {
        Particule &e = pair.second;
        Vecteur spe = e.get_speed() + (forces_old_leonnard.at(e.get_id()) + forces_new_leonnard.at(e.get_id())) * eta_t *  (0.5 / e.get_mass());
        e.set_speed(spe);
    }
}


/*
----------------------------------------------------
----------------------------------------------------
----------------------------------------------------
LAB 5
----------------------------------------------------
----------------------------------------------------
----------------------------------------------------
*/

void Univers::state_vis(int ds) {
    std::ofstream file("state" + std::to_string(ds) + ".vtu");

    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n";
    file << "<UnstructuredGrid>\n";
    file << "<Piece NumberOfPoints=\"" << this->get_nb() << "\" NumberOfCells=\"0\">\n";

    file << "<Points>\n";
    file << "<DataArray name=\"Position\" type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for(int i = 0; i < this->get_nb(); i++){
        file << this->get_position(i).getX() << " " << this->get_position(i).getY() << " " << this->get_position(i).getZ() << " ";
    }
    file << "\n" << "</DataArray>\n" << "</Points>\n";

    file << "<PointData Vectors=\"vector\">\n";
    file << "<DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for(int i = 0; i < this->get_nb(); i++){
        file << this->get_speed(i).getX() << " " << this->get_speed(i).getY() << " " << this->get_speed(i).getZ() << " ";
    }
    file << "\n" << "</DataArray>\n";

    file << "<DataArray type=\"Float32\" Name=\"Masse\" format=\"ascii\">\n";
    for(int i = 0; i < this->get_nb(); i++){
        file << this->get_mass(i) << " ";
    }
    file << "\n" << "</DataArray>\n" << "</PointData>";

    file << "<Cells>\n"
            "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n"
            "        </DataArray>\n"
            "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n"
            "        </DataArray>\n"
            "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n"
            "        </DataArray>\n"
            "      </Cells>";
    file << "</Piece>\n" << "</UnstructuredGrid>\n" << "</VTKFile>";
}



/*
----------------------------------------------------
----------------------------------------------------
----------------------------------------------------
LAB 6
----------------------------------------------------
----------------------------------------------------
----------------------------------------------------
*/

std::array<double, 6> Univers::getLimitCoordinates() {
    Vecteur bottom_left = this->grid.grid[0][0][0].getBottomLeft();
    double xmin = bottom_left.getX();
    double ymin = bottom_left.getY();
    double zmin = bottom_left.getZ();

    int nx = this->grid.grid.size();
    int ny = this->grid.grid[0].size();
    int nz = this->grid.grid[0][0].size();

    Vecteur cell_size = this->grid.grid[0][0][0].getCellSize();
    Vecteur top_right = this->grid.grid[nx-1][ny-1][nz-1].getBottomLeft() + cell_size;
    double xmax = top_right.getX();
    double ymax = top_right.getY();
    double zmax = top_right.getZ();

    std::array<double, 6> limits = {xmin, xmax, ymin, ymax, zmin, zmax};
    return limits;
}


void Univers::next_cond_limites(int k) {
    t += eta_t;
    std::unordered_map<int, Vecteur> forces_old_leonnard = force_field_grid();

    for(auto &pair : particules){
        Particule &e = pair.second;
        Particule original_e = e;
        Particule p = e;
        Vecteur pos = e.get_position() + (e.get_speed() + forces_old_leonnard.at(e.get_id()) * (0.5 / e.get_mass()) * eta_t) * eta_t;

        p.set_position(pos);
        Vecteur new_pos = Vecteur();
        if (k == 1)       new_pos = handleReflexion(e, p); // this modifies e, that's why need original
        else if (k == 2)  new_pos = handleTransportation(e,p);
        else              new_pos = handleAbsorption(e, p);
        p.set_position(new_pos);
        
        if (k == 1 || k == 2) this->grid.updateCell(original_e,p);
        e = p;
    }

    std::unordered_map<int, Vecteur> forces_new_leonnard = force_field_grid();
    for (auto &pair: particules) {
        Particule &e = pair.second;
        Vecteur spe = e.get_speed() + (forces_old_leonnard.at(e.get_id()) + forces_new_leonnard.at(e.get_id())) * eta_t *  (0.5 / e.get_mass());
        e.set_speed(spe);
    }
}

Vecteur Univers::handleReflexion(Particule &old, Particule &new_pos) {
    double x0 = old.get_position().getX();
    double y0 = old.get_position().getY();
    double z0 = old.get_position().getZ();

    double x1 = new_pos.get_position().getX();
    double y1 = new_pos.get_position().getY();
    double z1 = new_pos.get_position().getZ();

    std::array<double, 6> limits = this->getLimitCoordinates();
    double xmin = limits[0];         double xmax = limits[1];
    double ymin = limits[2];         double ymax = limits[3];
    double zmin = limits[4];         double zmax = limits[5];

    Vecteur intermediaire_potentiel = Vecteur();

    bool collision_x = false;
    bool collision_y = false;
    bool collision_z = false;

    double t_x = std::numeric_limits<double>::max();
    double t_y = std::numeric_limits<double>::max();
    double t_z = std::numeric_limits<double>::max();

    if (x1 < xmin) {
        intermediaire_potentiel.setX(xmin);
        collision_x = true;
        /*
        xmin = x0 + tc(x1-x0)
        tc = (xmin - x0) / (x1 - x0)
        */
        t_x = (xmin - x0) / (x1 - x0);
        
    } else if (x1 > xmax) {
        intermediaire_potentiel.setX(xmax);
        collision_x = true;
        t_x = (xmax - x0) / (x1 - x0);
    } else {
        intermediaire_potentiel.setX(x1);
    }


    if (y1 < ymin) {
        intermediaire_potentiel.setY(ymin);
        collision_y = true;
        t_y = (ymin - y0) / (y1 - y0);
    } else if (y1 > ymax) {
        intermediaire_potentiel.setY(ymax);
        collision_y = true;
        t_y = (ymax - y0) / (y1 - y0);
    } else {
        intermediaire_potentiel.setY(y1);
    }


    if (z1 < zmin) {
        intermediaire_potentiel.setZ(zmin);
        collision_z = true;
        t_z = (zmin - z0) / (z1 - z0);
    } else if (z1 > zmax) {
        intermediaire_potentiel.setZ(zmax);
        collision_z = true;
        t_z = (zmax - z0) / (z1 - z0);
    } else {
        intermediaire_potentiel.setZ(z1);
    }

    Vecteur intermediaire_final = Vecteur();
    bool x = false;
    bool y = false;
    bool z = false;
    if (collision_x || collision_y || collision_z) {
        if (t_x <= t_y && t_x <= t_z) {
            x = true;
            intermediaire_final = Vecteur(intermediaire_potentiel.getX(), y0 + (y1 - y0) * t_x, z0 + (z1 - z0) * t_x);
        }

        else if (t_y <= t_x && t_y <= t_z) {
            y = true;
            intermediaire_final = Vecteur(x0 + (x1 - x0) * t_y, intermediaire_potentiel.getY(), z0 + (z1 - z0) * t_y);
        }

        else if (t_z <= t_x && t_z <= t_y) {
            z = true;
            intermediaire_final = Vecteur(x0 + (x1 - x0) * t_z, y0 + (y1 - y0) * t_z, intermediaire_potentiel.getZ());
        }

        Vecteur new_diff = new_pos.get_position() - intermediaire_final;
        if (x) new_diff.setX(-new_diff.getX());
        if (y) new_diff.setY(-new_diff.getY());
        if (z) new_diff.setZ(-new_diff.getZ());

        Vecteur new_new = intermediaire_final + new_diff;
        new_pos.set_position(new_new);
        old.set_position(intermediaire_final);

        return handleReflexion(old, new_pos);
    }

    else {
        intermediaire_final = intermediaire_potentiel;
    }

    return intermediaire_final;
}

Vecteur Univers::handleTransportation(Particule &old, Particule &new_pos) {
    double x1 = new_pos.get_position().getX();
    double y1 = new_pos.get_position().getY();
    double z1 = new_pos.get_position().getZ();

    std::array<double, 6> limits = this->getLimitCoordinates();
    double xmin = limits[0];         double xmax = limits[1];
    double ymin = limits[2];         double ymax = limits[3];
    double zmin = limits[4];         double zmax = limits[5];

    Vecteur other_side = Vecteur();

  
    if (x1 < xmin)          other_side.setX(xmax - (xmin - x1));
    else if (x1 > xmax)     other_side.setX(xmin + (x1 - xmax));
    else                    other_side.setX(x1);

    if (y1 < ymin)          other_side.setY(ymax - (ymin - y1));
    else if (y1 > ymax)     other_side.setY(ymin + (y1 - ymax));
    else                    other_side.setY(y1);

    if (z1 < zmin)          other_side.setZ(zmax - (zmin - z1));
    else if (z1 > zmax)     other_side.setZ(zmin + (z1 - zmax));
    else                    other_side.setZ(z1);

    return other_side;
}

Vecteur Univers::handleAbsorption(Particule &old, Particule &new_pos) {
    double x1 = new_pos.get_position().getX();
    double y1 = new_pos.get_position().getY();
    double z1 = new_pos.get_position().getZ();

    std::array<double, 6> limits = this->getLimitCoordinates();
    double xmin = limits[0];         double xmax = limits[1];
    double ymin = limits[2];         double ymax = limits[3];
    double zmin = limits[4];         double zmax = limits[5];

    Cell& containing_cell = this->grid.containingCell(old);
    
    if (x1 < xmin || x1 > xmax || y1 < ymin || y1 > ymax || z1 < zmin || z1 > zmax) {
        containing_cell.remove(old);
        particules.erase(old.get_id());
        return Vecteur(0,0,0);
    }

    else {
        return new_pos.get_position();
    }
}

double Univers::computeForceBordure(Particule &p, double rcut, double min, double max, int i) {
    double coordonne = p.get_position()(i);

    double distance_min = coordonne - min + 1e-4;
    double distance_max = max - coordonne + 1e-4;

    double res = 0;
    if (distance_min < rcut) {
        // if distance is very small, it needs to rebound in direction ex. res is positive
        res += -24 * epsilon * (1/ (2*distance_min)) * pow((sigma/(2*distance_min)), 6) * (1 - 2 * pow(sigma/(2*distance_min), 6));
    }
    if (distance_max < rcut) {
        // if distance is very small, it needs to rebound in direction -ex. res is positive
        res += 24 * epsilon * (1/ (2*distance_max)) * pow((sigma/(2*distance_max)), 6) * (1 - 2 * pow(sigma/(2*distance_max), 6));
    }
    return res;
}

Vecteur Univers::forceBordures(Particule &p) {

    std::array<double, 6> limits = this->getLimitCoordinates();
    double xmin = limits[0];         double xmax = limits[1];
    double ymin = limits[2];         double ymax = limits[3];
    double zmin = limits[4];         double zmax = limits[5];

    double rcut_bordure = 1.12246204831*sigma;

    Vecteur x = Vecteur(1,0,0);
    Vecteur y = Vecteur(0,1,0);
    Vecteur z = Vecteur(0,0,1);

    double dir_x = computeForceBordure(p, rcut_bordure, xmin, xmax, 0);
    Vecteur forceX = dir_x * x;

    Vecteur forceY = Vecteur(0,0,0);
    if (this->dim >= 2) {
        double dir_y = computeForceBordure(p, rcut_bordure, ymin, ymax, 1);
        forceY = dir_y * y;
    }

    Vecteur forceZ = Vecteur(0,0,0);
    if (this->dim == 3) {
        double dir_z = computeForceBordure(p, rcut_bordure, zmin, zmax, 2);
        forceZ = dir_z * z;
    }

    return forceX + forceY + forceZ;
}

std::unordered_map<int, Vecteur> Univers::force_field_potentiel_reflexion() {
    std::unordered_map<int, Vecteur> forces;

    for(auto &pair : particules) {
        Particule &e = pair.second;
        Vecteur force = forceBordures(e) - e.get_mass()*this->g*Vecteur(0,1,0);
        Cell cell = grid.containingCell(e);
        std::vector<Cell> neighbors = grid.getNeighbors(cell);
        for (Cell &c: neighbors) {
            Vecteur cell_center = c.getCenter();
            if (cell_center.distance(e.get_position()) < rcut) {
                for (auto &p : c.particles) {
                    if (p.first != e.get_id()) {
                        force += e.computeLennardForce(p.second, sigma, epsilon);
                    }
                }
            }
        }
        forces.insert({e.get_id(), force});
    }
    return forces;
}




void Univers::next_potentiel_reflexion() {
    t += eta_t;
    std::unordered_map<int, Vecteur> forces_old_leonnard = force_field_potentiel_reflexion();

    for(auto &pair : particules){
        Particule &e = pair.second;
        Particule original_e = e;
        Vecteur pos = e.get_position() + (e.get_speed() + forces_old_leonnard.at(e.get_id()) * (0.5 / e.get_mass()) * eta_t) * eta_t;
        Particule p = e;
        p.set_position(pos);

        handleReflexion(e, p);
        this->grid.updateCell(e,p);
        e = p;
    }

    std::unordered_map<int, Vecteur> forces_new_leonnard = force_field_grid();
    for (auto &pair: particules) {
        Particule &e = pair.second;
        Vecteur spe = e.get_speed() + (forces_old_leonnard.at(e.get_id()) + forces_new_leonnard.at(e.get_id())) * eta_t *  (0.5 / e.get_mass());
        e.set_speed(spe);
    }
}