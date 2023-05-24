/**
 * @file vecteur.hxx
 * @brief Defines the Vecteur class for representing a 3D vector.
 * @date May 2023
 * @authors Achraf Kerzazi & Mohamed Errazki
 */

#ifndef LAB3_INFRA_VECTEUR_HXX
#define LAB3_INFRA_VECTEUR_HXX

#include <iostream>

/**
 * @brief Represents a 3D vector.
 */
class Vecteur {
private:
    double x; /**< The x-component of the vector. */
    double y; /**< The y-component of the vector. */
    double z; /**< The z-component of the vector. */

public:
    /**
     * @brief Default constructor for the Vecteur class.
     * @param x The x-component of the vector (default: 0).
     * @param y The y-component of the vector (default: 0).
     * @param z The z-component of the vector (default: 0).
     */
    explicit Vecteur(double x = 0, double y = 0, double z = 0);

    /**
     * @brief Gets the x-component of the vector.
     * @return The x-component of the vector.
     */
    double getX() const;

    /**
     * @brief Gets the y-component of the vector.
     * @return The y-component of the vector.
     */
    double getY() const;

    /**
     * @brief Gets the z-component of the vector.
     * @return The z-component of the vector.
     */
    double getZ() const;


    void setX(double x);
    void setY(double y);
    void setZ(double z);

    /**
     * @brief Calculates the distance between this vector and another vector.
     * @param other The other vector.
     * @return The distance between the two vectors.
     */
    double distance(const Vecteur &other);

    /**
     * @brief Calculates the norm (magnitude) of the vector.
     * @return The norm of the vector.
     */
    double norme();

    // Overloaded operators
    Vecteur operator+(const Vecteur &);
    Vecteur operator-(const Vecteur &);
    Vecteur &operator+=(const Vecteur &);
    Vecteur &operator-=(const Vecteur &);
    Vecteur operator*(double);
    Vecteur &operator*=(double);
    Vecteur &operator=(const Vecteur &);
    double &operator()(int);
    double operator()(int) const;

    // Friend functions
    friend std::ostream &operator<<(std::ostream &, const Vecteur &);
    friend Vecteur operator*(double, const Vecteur &);
    friend Vecteur operator+(const Vecteur &, const Vecteur &);
    friend Vecteur operator-(const Vecteur &, const Vecteur &);
};

#endif // LAB3_INFRA_VECTEUR_HXX
