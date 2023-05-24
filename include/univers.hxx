/**
 * @file univers.hxx
 * @brief Defines the Univers class for representing the universe in which particles interact.
 * @date May 2023
 * @authors Achraf Kerzazi & Mohamed Errazki
 */

#ifndef LAB3_UNIVERS_HXX
#define LAB3_UNIVERS_HXX

#include <iostream>
#include <vector>
#include <unordered_map>
#include "particule.hxx"
#include "vecteur.hxx"
#include "grid.hxx"
#include <array>

/**
 * @brief Represents the universe in which particles interact.
 */
class Univers {
private:
    std::unordered_map<int, Particule> particules; /**< The collection of particles in the universe. */
    int dim; /**< The dimension of the universe. */
    double t; /**< The current time in the universe. */
    double eta_t; /**< The time step size. */
    double sigma; /**< The sigma value for force calculations. */
    double epsilon; /**< The epsilon value for force calculations. */
    Vecteur ld; /**< Parameters that help define the grid. */
    double rcut; /**< The cutoff distance for force calculations. */
    Grid grid; /**< The grid used for spatial partitioning. */
    double g; /**< The acceleration due to gravity. */

public:

    /**
     * @brief Constructs an instance of the Univers class.
     * @param dim The dimension of the universe.
     * @param time The initial time (default: 0).
     * @param eta The time step size (default: 0.0015).
     */
    Univers(int dim, double time = 0, double eta = 0.0015);

    /**
     * @brief Constructs an instance of the Univers class.
     * @param dim The dimension of the universe.
     * @param t The initial time.
     * @param eta_t The time step size.
     * @param sigma The sigma value for force calculations.
     * @param epsilon The epsilon value for force calculations.
     * @param rcut The cutoff distance for force calculations.
     * @param ld The length and width of the universe.
     * @param particules The collection of particles.
     * @param g The acceleration due to gravity.
     */
    Univers(int dim, double t, double eta_t, double sigma, double epsilon, double rcut, Vecteur ld, std::unordered_map<int, Particule> particules, double g);

    /**
     * @brief Adds a particle to the universe.
     * @param p The particle to be added.
     */
    void add_particule(Particule &p);

    /**
     * @brief Returns the current time in the universe.
     * @return The current time.
     */
    double get_t();

    /**
     * @brief Sets the current time in the universe.
     * @param t The current time.
     */
    void set_t(double t);

    /**
     * @brief Returns the collection of particles in the universe.
     * @return The collection of particles.
     */
    std::unordered_map<int, Particule> get_particules();

    /**
     * @brief Returns the time step size.
     * @return The time step size.
     */
    double get_eta();

    /**
     * @brief Sets the time step size.
     * @param eta The time step size.
     */
    void set_eta(double eta);

    /**
     * @brief Sets the speed of a particle.
     * @param index The index of the particle.
     * @param v The new speed vector.
     */
    void set_speed(int index, Vecteur &v);

    /**
     * @brief Returns the position of a particle.
     * @param index The index of the particle.
     * @return The position vector.
     */
    Vecteur get_position(int index);


    /**
     * @brief Returns the speed of a particle.
     * @param index The index of the particle.
     * @return The speed vector.
     */
    Vecteur get_speed(int index);

    /**
     * @brief Returns the mass of a particle.
     * @param index The index of the particle.
     * @return The mass.
     */
    double get_mass(int index);


    /**
     * @brief Get the number of particles in the universe.
     * 
     * @return int 
     */
    int get_nb();

    void state_vis(int ds);

    /**
     * @brief Displays the current state of the universe.
     */
    void show_state();

    /**
     * @brief Returns the limit coordinates of the universe.
     * @return An array containing the minimum and maximum coordinates in each dimension.
     */
    std::array<double, 6> getLimitCoordinates();

    /**
     * @brief Calculates the force field acting on the particles using a naive approach.
     * @return A map of particle IDs to force vectors.
     */
    std::unordered_map<int, Vecteur> force_field();

    /**
     * @brief Calculates the force field acting on the particles using the grid method.
     * @return A map of particle IDs to force vectors.
     */
    std::unordered_map<int, Vecteur> force_field_grid();

    /**
     * @brief Calculates the force field acting on the particles using the potential reflection too.
     * @return A map of particle IDs to force vectors.
     */
    std::unordered_map<int, Vecteur> force_field_potentiel_reflexion();

    /**
     * @brief Advances the simulation to the next time step using the naive approach.
     */
    void next();

    /**
     * @brief Advances the simulation to the next time step using the grid method.
     */
    void next_grid();

    /**
     * @brief Advances the simulation to the next time step with boundary conditions.
     * @param k handles conditions limites using reflection if k = 1,
     *          transportation if 2, absorption if else.
     */
    void next_cond_limites(int k);

    /**
     * @brief Advances the simulation to the next time step using the potential reflection method.
     */
    void next_potentiel_reflexion();

    /**
     * @brief Handles the reflection of particles at the boundaries.
     * @param p The particle to handle.
     * @param e The new particle, same one but with different position.
     * @return The force vector acting on the particle due to the reflection.
     */
    Vecteur handleReflexion(Particule &p, Particule &e);

    /**
     * @brief Handles the transportation of particles across boundaries.
     * @param p The particle to handle.
     * @param border_p The particle representing the boundary.
     * @return The force vector acting on the particle due to the transportation.
     */
    Vecteur handleTransportation(Particule &p, Particule &border_p);

    /**
     * @brief Handles the absorption of particles by boundaries.
     * @param p The particle to handle.
     * @param border_p The particle representing the boundary.
     * @return The force vector acting on the particle due to the absorption.
     */
    Vecteur handleAbsorption(Particule &p, Particule &border_p);

    /**
     * @brief Computes the force acting on a particle near the boundaries.
     * @param p The particle.
     * @param rcut The cutoff distance.
     * @param min The minimum boundary coordinate.
     * @param max The maximum boundary coordinate.
     * @param i The dimension index.
     * @return The norm force acting on the particle in dimension i.
     */
    double computeForceBordure(Particule &p, double rcut, double min, double max, int i);

    /**
     * @brief Computes the total force acting on a particle due to the boundaries.
     * @param p The particle.
     * @return The force acting on the particle in all dimensions.
     */
    Vecteur forceBordures(Particule &p);
};

#endif // LAB3_UNIVERS_HXX
