#ifndef ELEMENT_MANAGER_HPP_
#define ELEMENT_MANAGER_HPP_

#include <cstddef>
#include <memory>
#include <string>

#include "nicole/config/input_config.hpp"
#include "nicole/element/element.hpp"
#include "nicole/nicole_defs.hpp"

namespace nicole {
    class ElementManager {
    public:
        ElementManager(InputConfig& input);

        // Getter
        std::string GetElementName(ElementID id) const;
        double GetElementMass(ElementID id) const;
        std::size_t GetNumberOfElements() const { return number_of_elements_; }

        ElementID FindIdElement(const std::string& element_name) const;
        void AddElement(std::shared_ptr<Element> element);

        void DisplayElements() const;

    private:
        void ReadElementFile(const std::string& filename);

        std::vector<std::shared_ptr<Element>> element_list_;
        std::size_t number_of_elements_;
    };
}

#endif /* ELEMENT_MANAGER_HPP */
