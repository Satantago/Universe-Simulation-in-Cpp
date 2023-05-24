/**
 * @file particule.hxx
 * @brief Defines the Particule class for representing a particle.
 * @date May 2023
 * @authors Achraf Kerzazi & Mohamed Errazki
 */

#ifndef LAB3_PARTICULE_HXX
#define LAB3_PARTICULE_HXX

#include "vecteur.hxx"
#include <iostream>
#include <vector>

/**
 * @brief Represents a particle.
 */
class Particule {
private:
    int id; /**< The particle ID. */
    Vecteur x; /**< The position vector of the particle. */
    Vecteur v; /**< The velocity vector of the particle. */
    double mass; /**< The mass of the particle. */
    std::string type; /**< The type of the particle. */

public:
    /**
     * @brief Returns the ID of the particle.
     * @return The ID of the particle.
     */
    int get_id();

    /**
     * @brief Returns the position vector of the particle.
     * @return The position vector.
     */
    Vecteur get_position() const;

    /**
     * @brief Sets the position vector of the particle.
     * @param v The new position vector.
     */
    void set_position(Vecteur &v);

    /**
     * @brief Returns the velocity vector of the particle.
     * @return The velocity vector.
     */
    Vecteur get_speed();

    /**
     * @brief Sets the velocity vector of the particle.
     * @param v The new velocity vector.
     */
    void set_speed(Vecteur &v);

    /**
     * @brief Returns the mass of the particle.
     * @return The mass of the particle.
     */
    double get_mass();

    /**
     * @brief Constructs a particle with the given ID and mass.
     * @param id The ID of the particle.
     * @param mass The mass of the particle (default: 1).
     */
    explicit Particule(int id, double mass = 1);

    /**
     * @brief Constructs a particle with the given ID, position, velocity, and mass.
     * @param id The ID of the particle.
     * @param x The position vector of the particle.
     * @param v The velocity vector of the particle.
     * @param mass The mass of the particle (default: 1).
     */
    explicit Particule(int id, Vecteur &x, Vecteur &v, double mass = 1);

    /**
     * @brief Calculates the distance between this particle and another particle.
     * @param p The other particle.
     * @return The distance between the particles.
     */
    double distance(const Particule &p);

    /**
     * @brief Calculates the Lennard-Jones potential between this particle and another particle.
     * @param p The other particle.
     * @param epsilon The epsilon value for the Lennard-Jones potential.
     * @param sigma The sigma value for the Lennard-Jones potential.
     * @return The Lennard-Jones potential energy.
     */
    double lennardJonesPotential(const Particule &p, double epsilon, double sigma);

    /**
     * @brief Computes the Lennard-Jones force vector between this particle and another particle.
     * @param p The other particle.
     * @param epsilon The epsilon value for the Lennard-Jones potential.
     * @param sigma The sigma value for the Lennard-Jones potential.
     * @return The force vector between the particles.
     */
    Vecteur computeLennardForce(const Particule &p, double epsilon, double sigma);

    /**
     * @brief Computes the total Lennard-Jones force vector between this particle and a vector of particles.
     * @param particles The vector of particles.
     * @param epsilon The epsilon value for the Lennard-Jones potential.
     * @param sigma The sigma value for the Lennard-Jones potential.
     * @return The total force vector between this particle and the vector of particles.
     */
    Vecteur computeLennardAll(std::vector<Particule> particles, double epsilon, double sigma);
};

#endif // LAB3_PARTICULE_HXX
