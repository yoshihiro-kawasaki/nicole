/**
 * @file element_manager.hpp
 * @brief Manages a collection of chemical elements
 * @date 2025-02-12
 * @author Y. Kawasaki
 */

#ifndef ELEMENT_MANAGER_HPP
#define ELEMENT_MANAGER_HPP

#include "../nicole_defs.hpp"
#include "element.hpp"

namespace nicole
{
    /**
     * @class ElementManager
     * @brief A class to manage multiple elements and provide access functions.
     */
    class ElementManager
    {
    public:

        /**
         * @brief Constructor that initializes elements from an input file.
         * @param input Reference to InputConfig, which provides configuration details.
         */
        explicit ElementManager(InputConfig& input);

        /**
         * @brief Get the name of an element by its index.
         * @param index The index of the element.
         * @return The name of the element. Returns an empty string if not found.
         */
        std::string GetElementName(std::size_t index) const;

        /**
         * @brief Get the atomic mass of an element by its index.
         * @param index The index of the element.
         * @return The atomic mass in atomic mass units (amu). Returns 0.0 if not found.
         */
        double GetElementMass(std::size_t index) const;

        /**
         * @brief Get the total number of elements managed.
         * @return The number of elements.
         */
        std::size_t GetNumberOfElements() const { return number_of_elements_; }

        /**
         * @brief Display all elements currently stored.
         */
        void DisplayElements() const;

        /**
         * @brief Find the index of an element given its name.
         * @param element_name The name of the element to search for.
         * @return The index of the element, or 9999 (kNotFoundSpecies) if not found.
         */
        int FindIndexElement(const std::string& element_name) const;

        /**
         * @brief Add a new element to the manager.
         * @param element A shared pointer to the element to be added.
         */
        void AddElement(std::shared_ptr<Element> element);

    private:

        /**
         * @brief Read element data from a file.
         * @param filename The name of the file to read.
         */
        void ReadElementFile(const std::string& filename);

        // List of elements
        std::vector<std::shared_ptr<Element>> element_list_;

        // Number of elements managed
        std::size_t number_of_elements_;

    };
}

#endif /* ELEMENT_MANAGER_HPP */