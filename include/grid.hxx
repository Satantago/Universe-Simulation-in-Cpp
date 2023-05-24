/**
 * @file   grid.hxx
 * @brief Defines the Grid class and the Cell class for managing a grid of cells in a simulation.
 * @authors Achraf Kerzazi & Mohamed Errazki 
 * @date   May 2023
 */ 



#ifndef GRID_HPP
#define GRID_HPP

#include <unordered_map>
#include <vector>
#include "particule.hxx"
#include "vecteur.hxx"

/**
 * @brief Represents a cell in the grid.
*/
class Cell {
    private:
        Vecteur bottom_left; /**< The bottom-left position of the cell.*/
        int i,j,k;           /**< The indices of the cell in the grid.*/
        Vecteur cell_size;  /**< The size of the cell.*/

    public:
        /**
         * @brief Constructs a Cell object with specified indices, cell size, and bottom-left position.
         * @param i The index in the x-axis.
         * @param j The index in the y-axis.
         * @param k The index in the z-axis.
         * @param cell_size The size of the cell.
         * @param bottom_left The bottom-left position of the cell.
        */
        Cell(int i, int j, int k, Vecteur cell_size, Vecteur bottom_left);

        /**
         * @brief Constructs a Cell object with default values.
        */
        Cell();

        std::unordered_map<int, Particule> particles; /**< The particles contained within the cell. */
        
        /**
         * @brief Checks if the cell contains a given particle.
         * @param particule The particle to check.
         * @return True if the cell contains the particle, false otherwise.
        */
        bool contains(Particule &);

        /**
         * @brief Inserts a particle into the cell.
         * @param particule The particle to insert.
        */
        void insert(Particule &);

        /**
         * @brief Removes a particle from the cell.
         * @param particule The particle to remove.
        */
        void remove(Particule &);

        /**
         * @brief Gets the center of the cell.
         * @return The center of the cell.
        */
        Vecteur getCenter();

        /**
         * @brief Gets the bottom-left position of the cell.
         * @return The bottom-left position of the cell.
        */
        Vecteur getBottomLeft();

        /**
         * @brief Gets the size of the cell.
         * @return The size of the cell.
        */
        Vecteur getCellSize();

        /**
         * @brief Gets the indices of the cell in the grid.
         * @return The indices of the cell in the grid.
        */
        Vecteur getIndices();


};

/**
 * @brief Represents a grid of cells.
*/
class Grid {
public:
    Vecteur ld; /**< The length and width of the grid. */
    double rcut; /**< The cutoff distance. */
    std::unordered_map<int, Particule> particles; /**< The particles in the grid. */
    std::vector<std::vector<std::vector<Cell>>> grid; /**< The cells of the grid. */

    /**
     * @brief Constructs a Grid object with specified length, width, cutoff distance, and particles.
     * @param ld Parameters of the universe.
     * @param rcut The cutoff distance.
     * @param particles The particles in the grid.
     */
    Grid(Vecteur ld, double rcut, std::unordered_map<int, Particule> particles);

    /**
     * @brief Default constructor for the Grid class.
     */
    Grid();

    /**
     * @brief Returns the neighboring cells of a given cell.
     * @param c The cell for which to retrieve the neighbors.
     * @return The neighboring cells of the given cell.
     */
    std::vector<Cell> getNeighbors(Cell &c);

    /**
     * @brief Returns the size of each cell in the grid.
     * @return The size of each cell.
     */
    Vecteur getCellSize();

    /**
     * @brief Returns the cell that contains a given particle.
     * @param p The particle to find the containing cell.
     * @return The cell that contains the given particle.
     */
    Cell& containingCell(Particule &p);

    /**
     * @brief Updates the cell that contains a given particle after its position has changed.
     * @param old_particle The particle before the update.
     * @param new_particle The particle after the update.
     */
    void updateCell(Particule &old_particle, Particule &new_particle);
};

#endif