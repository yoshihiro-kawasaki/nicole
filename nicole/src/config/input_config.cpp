#include <fstream>
#include <iostream>

#include "nicole/config/input_config.hpp"
#include "nicole/utils/string_utils.hpp"

namespace nicole {
    InputConfig::InputConfig(const std::string& filename) : input_filename_(filename) {
        ReadFile(filename);
    }


    void InputConfig::ReadFile(const std::string& filename) {
        std::ifstream file(filename, std::ios::in);
        if (!file) {
            std::cerr << "Error opening file: " << filename << std::endl;
            return;
        }

        std::string line;
        while (std::getline(file, line)) {
            if (line[0] == '#' || line[0] == '!' || line.empty()) continue;

            std::vector<std::string> split_line = string_utils::Split(line, ' ', true);
            if (line.size() >= 3) {
                dict_[split_line[0]] = split_line[2];
            } else {
                std::cerr << "Invalid line in file: " << line << std::endl;
            }
        } 

        file.close();
    }


    bool InputConfig::Contains(const std::string& key) const {
        return dict_.find(key) != dict_.end();
    }


    void InputConfig::Insert(const std::string& key, const std::string& item) {
        dict_[key] = item;
    }


    const std::string& InputConfig::GetString(const std::string& key) const {
        auto it = dict_.find(key);
        if (it != dict_.end()) {
            return it->second;
        } else {
            throw std::runtime_error("Key not found: " + key);
        }
    }


    Real InputConfig::GetReal(const std::string& key) const {
        auto it = dict_.find(key);
        if (it != dict_.end()) {
            try {
                return std::stod(it->second); // Convert string to double
            } catch (const std::invalid_argument& e) {
                throw std::runtime_error("Invalid double value for key: " + key + ", value: " + it->second);
            } catch (const std::out_of_range& e) {
                throw std::runtime_error("Out of range double value for key: " + key + ", value: " + it->second);
            }
        } else {
            throw std::runtime_error("Key not found: " + key); // Key not found
        }
    }


    int InputConfig::GetInt(const std::string& key) const {
        auto it = dict_.find(key);
        if (it != dict_.end()) {
            try {
                return std::stoi(it->second); // Convert string to int
            } catch (const std::invalid_argument& e) {
                throw std::runtime_error("Invalid int value for key: " + key + ", value: " + it->second);
            } catch (const std::out_of_range& e) {
                throw std::runtime_error("Out of range int value for key: " + key + ", value: " + it->second);
            }
        } else {
            throw std::runtime_error("Key not found: " + key); // Key not found
        }
    }


    std::size_t InputConfig::GetSizeT(const std::string& key) const {
        auto it = dict_.find(key);
        if (it != dict_.end()) {
            try {
                return std::stoull(it->second); // Convert string to size_t
            } catch (const std::invalid_argument& e) {
                throw std::runtime_error("Invalid size_t value for key: " + key + ", value: " + it->second);
            } catch (const std::out_of_range& e) {
                throw std::runtime_error("Out of range size_t value for key: " + key + ", value: " + it->second);
            }
        } else {
            throw std::runtime_error("Key not found: " + key); // Key not found
        }
    }


    bool InputConfig::GetBool(const std::string& key) const {
        auto it = dict_.find(key);
        if (it != dict_.end()) {
            return string_utils::ConvertStringToBool(it->second);
        } else {
            throw std::runtime_error("Key not found: " + key); // Key not found
        }
    }
}
