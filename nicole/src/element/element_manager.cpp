#include <fstream>
#include <iomanip>
#include <iostream>

#include "nicole/element/element_manager.hpp"
#include "nicole/utils/string_utils.hpp"

namespace nicole {
    ElementManager::ElementManager(InputConfig& input) : number_of_elements_(0) {
        ReadElementFile(input.GetString("element_file"));
    }


    std::string ElementManager::GetElementName(ElementID id) const {
        if (id < element_list_.size()) {
            return element_list_[id]->GetName();
        }
        return nicole::kEmptyString;
    }


    Real ElementManager::GetElementMass(ElementID id) const {
        if (id < element_list_.size()) {
            return element_list_[id]->GetMass();
        }
        return 0.0;
    }


    ElementID ElementManager::FindIdElement(const std::string& element_name) const {
        for (std::size_t i = 0; i < element_list_.size(); ++i) {
            if (element_list_[i]->GetName() == element_name) {
                return static_cast<ElementID>(i);
            }
        }
        return kNotFoundElement;
    }


    void ElementManager::AddElement(std::shared_ptr<Element> element) {
        element_list_.push_back(element);
        number_of_elements_++;  // 元素を追加した際にカウントを増やす
    }


    void ElementManager::DisplayElements() const {
        for (const auto& element : element_list_) {
            std::cout << "ID: "   << std::setw(4) << element->GetID()   << "  "
                      << "Name: " << std::setw(4) << element->GetName() << "  "
                      << "Mass: " << std::setw(4) << element->GetMass() 
                      << std::endl;
        }
    }


    void ElementManager::ReadElementFile(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return;
        }

        std::string line;
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#' || line[0] == '!') continue;

            std::string line;
            while (std::getline(file, line)) {
                if (line.empty() || line[0] == '#' || line[0] == '!') continue;

                std::vector<std::string> split_line = string_utils::Split(line, ' ');
                std::string element_name = split_line[0];
                Real element_mass = std::stoi(split_line[1]);

                element_list_.push_back(std::make_shared<Element>(element_name, element_mass));
            }
        }

        number_of_elements_ = element_list_.size();
        return;
    }
} // namespace nicole
