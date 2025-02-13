/**
 * @file element.cpp
 * @brief Implementation of the 'Element' class.
 * 
 * This source file contains the implementation of the `Element` class, which manages
 * the data related to a chemical element. The class holds the element's name, mass, 
 * and a unique identifier. This file defines the constructor that initializes the 
 * element's data and assigns it a unique ID.
 * 
 * @date 2025-02-11
 * @author Y. Kawasaki
 */

#include "element.hpp"

namespace nicole
{
    // 静的メンバーの初期化
    std::size_t Element::number_of_total_elements_ = 0;

    /**
     * @brief Constructs an Element object with a specified name and mass.
     * 
     * The constructor initializes the element with the provided name and mass, and 
     * assigns a unique ID based on the total number of elements created so far. The 
     * `number_of_total_elements_` static variable is incremented with each new element 
     * creation to ensure each element has a unique identifier.
     * 
     * @param name The name of the element.
     * @param mass The mass of the element in amu.
     */
    Element::Element(const std::string& name, const double mass)
            : name_(name), mass_(mass)
        { 
            number_of_total_elements_++;
            id_ = number_of_total_elements_;
        }
}