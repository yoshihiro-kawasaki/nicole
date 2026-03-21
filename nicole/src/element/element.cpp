#include "nicole/element/element.hpp"

namespace nicole {
    std::size_t Element::number_of_total_elements_ = 0;


    Element::Element(const std::string& name, const Real mass)
        : name_(name), mass_(mass)
    { 
        number_of_total_elements_++;
        id_ = number_of_total_elements_;
    }
}
