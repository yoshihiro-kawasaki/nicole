/**
 * @file element.hpp
 * @brief Defines the Element class for managing elemental data.
 * 
 * This header file defines the `Element` class, which is responsible for holding data 
 * related to a chemical element, such as its name, mass, and a unique identifier. 
 * 
 * @date 2025-02-11
 * @author Y. Kawasaki
 */
#ifndef ELEMENT_HPP
#define ELEMENT_HPP


#include "../nicole_defs.hpp"


namespace nicole
{

    /**
     * @class Element
     * @brief A class that holds data for a chemical element.
     * 
     * The `Element` class stores the name, mass, and a unique identifier for each element. 
     * It provides getter methods to retrieve the element's data and ensures that each 
     * element has a unique ID.
     * 
     * The class also manages the total number of elements created through the static member 
     * `number_of_total_elements_`.
     */
    class Element
    {
    public:

        /**
         * @brief Constructs an Element object with the specified name and mass.
         * 
         * This constructor initializes the element's name, mass, and assigns a unique ID.
         * 
         * @param name The name of the element
         * @param mass The mass of the element in amu
         */
        Element(const std::string& name, const double mass);

        /**
         * @brief Gets the unique identifier of the element.
         * 
         * Each element is assigned a unique ID when it is created. This method retrieves 
         * that ID.
         * 
         * @return The unique ID of the element.
         */
        std::size_t GetID() const { return id_; }

        /**
         * @brief Gets the name of the element.
         * 
         * This method returns the name of the element.
         * 
         * @return The name of the element.
         */
        const std::string& GetName() const { return name_; }

        /**
         * @brief Gets the mass of the element.
         * 
         * This method returns the mass of the element in amu.
         * 
         * @return The mass of the element.
         */
        double GetMass() const { return mass_; }

    protected:

        // Unique identifier for the element.
        std::size_t id_;

        // The name of the element.
        std::string name_;
        
        // The mass of the element in amu.
        double mass_;

    private:

        // Total number of elements created.
        static std::size_t number_of_total_elements_;

    };

}

#endif /* ELEMENT_HPP */