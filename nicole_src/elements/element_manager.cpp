/**
 * @file element_manager.cpp
 * @brief Implementation of class ElementManager
 * @date 2025-02-12
 * @author Y. Kawasaki
 */

#include "element_manager.hpp"

namespace nicole
{
    /**
     * @brief Constructor that initializes elements from an input file.
     * @param input Reference to InputConfig, which provides the element file path.
     */
    ElementManager::ElementManager(InputConfig& input)
        : number_of_elements_(0)
    {
        std::string element_filename = input.GetString("element_file");
        ReadElementFile(element_filename);
    }

    /**
     * @brief Get the name of an element by its index.
     * @param index The index of the element.
     * @return The element name, or an empty string if index is out of range.
     */
    std::string ElementManager::GetElementName(std::size_t index) const 
    {
        if (index < element_list_.size()) {
            return element_list_[index]->GetName();
        }
        return nicole::kEmptyString;
    }

    /**
     * @brief Get the atomic mass of an element by its index.
     * @param index The index of the element.
     * @return The atomic mass in atomic mass units (amu), or 0.0 if index is out of range.
     */
    double ElementManager::GetElementMass(std::size_t index) const 
    {
        if (index < element_list_.size()) {
            return element_list_[index]->GetMass();
        }
        return 0.0;
    }

    /**
     * @brief Display all elements currently stored.
     */
    void ElementManager::DisplayElements() const {
        for (const auto& element : element_list_) {
            std::cout << "ID: "   << std::setw(4) << element->GetID()   << "  "
                      << "Name: " << std::setw(4) << element->GetName() << "  "
                      << "Mass: " << std::setw(4) << element->GetMass() 
                      << std::endl;
        }
    }

    /**
     * @brief Find the index of an element by its name.
     * @param element_name The name of the element to search for.
     * @return The index of the element, or 9999 (kNotFoundSpecies) if not found.
     */
    int ElementManager::FindIndexElement(const std::string& element_name) const 
    {
        for (std::size_t i = 0; i < element_list_.size(); ++i) {
            if (element_list_[i]->GetName() == element_name) {
                return static_cast<int>(i);
            }
        }
        return kNotFoundElement;  // 見つからなかった場合
    }

    /**
     * @brief Add a new element to the manager.
     * @param element A shared pointer to the element to be added.
     */
    void ElementManager::AddElement(std::shared_ptr<Element> element) 
    {
        element_list_.push_back(element);
        number_of_elements_++;  // 元素を追加した際にカウントを増やす
    }

    /**
     * @brief Read element data from a file.
     * @param filename The name of the file to read.
     */
    void ElementManager::ReadElementFile(const std::string& filename)
    {
        // open file
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return;
        }

        std::string line;
        while (std::getline(file, line)) { // read file

            // Skip empty lines or comment lines
            if (line.empty() || line[0] == '#' || line[0] == '!') continue;

            // 空白で分割
            std::istringstream iss(line);
            std::string element_name;
            double element_mass;

            // Read element name and mass
            if (!(iss >> element_name >> element_mass)) {
                std::cerr << "Warning: Invalid line format: " << line << std::endl;
                continue; // 不正な行はスキップ
            }

            // Store element
            element_list_.emplace_back(std::make_shared<Element>(
                element_name, 
                element_mass
                )
            );
        }

        // Update element count
        number_of_elements_ = element_list_.size();
    }
    
} // namespace nicole