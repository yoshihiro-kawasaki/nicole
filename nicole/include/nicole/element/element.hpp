#ifndef ELEMENT_HPP_
#define ELEMENT_HPP_

#include "nicole/nicole_defs.hpp"

namespace nicole {
    class Element {
    public:
        Element(const std::string& name, const Real mass);
        ElementID GetID() const { return id_; }
        const std::string& GetName() const { return name_; }
        Real GetMass() const { return mass_; }
    private:
        ElementID id_;
        std::string name_;
        Real mass_;
        static std::size_t number_of_total_elements_;
    };
}

#endif /* ELEMENT_HPP_ */
