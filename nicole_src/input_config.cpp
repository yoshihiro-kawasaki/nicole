/**
 * @file input_config.cpp
 * @brief implimentation of class InputConfig
 * @date 2025-02-11
 * @author Y. Kawasaki
 */


#include "input_config.hpp"


/**
 * @brief Constructs an InputConfig object and reads the configuration file.
 * 
 * This constructor initializes the input_filename_ member variable with the given filename
 * and then calls the ReadFile function to parse the contents of the file.
 * 
 * @param filename The name of the configuration file to be read.
 */
InputConfig::InputConfig(const std::string& filename)
    : input_filename_(filename) 
{
    ReadFile(filename);
}


/**
 * @brief Reads the configuration file and stores key-value pairs.
 * 
 * This function opens the specified configuration file, processes each line, and stores the
 * key-value pairs found in the file into the internal dictionary (dict_).
 * 
 * Lines that start with '#' or '!' are treated as comments and are ignored. Empty lines are also skipped.
 * The file is expected to contain key-value pairs separated by spaces.
 * 
 * @param filename The name of the configuration file to read from.
 * 
 * @throws std::runtime_error If the file cannot be opened or if the file format is incorrect.
 */
void InputConfig::ReadFile(const std::string& filename) 
{
    // file open
    std::ifstream file(filename, std::ios::in);
    if (!file) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // read file
    std::string str;
    while (std::getline(file, str)) {
        // Lines that start with '#' or '!' are treated as comments and are ignored.
        // Empty lines are also skipped.
        if (str[0] == '#' || str[0] == '!' || str.empty()) continue;

        std::vector<std::string> line = string_utils::Split(str, ' ', true);
        if (line.size() >= 3) {
            dict_[line[0]] = line[2];
        } else {
            std::cerr << "Invalid line in file: " << str << std::endl;
        }
    } 

    file.close();
}


/**
 * @brief Checks if the key exists in the configuration dictionary.
 * 
 * This function checks whether a given key is present in the configuration dictionary.
 * It returns true if the key exists, false otherwise.
 * 
 * @param key The key to check for existence.
 * @return True if the key exists, false otherwise.
 */
bool InputConfig::Contains(const std::string& key) const 
{
    return dict_.find(key) != dict_.end();
}


/**
 * @brief Inserts a key-value pair into the configuration dictionary.
 * 
 * This function inserts or updates a key-value pair in the configuration dictionary.
 * If the key already exists, its value will be updated.
 * 
 * @param key The key to insert or update.
 * @param item The value to associate with the key.
 */
void InputConfig::Insert(const std::string& key, const std::string& item) 
{
    dict_[key] = item;
}


/**
 * @brief Retrieves a string value associated with a given key.
 * 
 * This function retrieves the value corresponding to the specified key from the
 * configuration dictionary. If the key does not exist, an exception is thrown.
 * 
 * @param key The key whose value is to be retrieved.
 * @return The string value associated with the key.
 * 
 * @throws std::runtime_error If the key is not found in the dictionary.
 */
const std::string& InputConfig::GetString(const std::string& key) 
{
    auto it = dict_.find(key);
    if (it != dict_.end()) {
        return dict_[key];
    } else {
        throw std::runtime_error("Key not found: " + key);
    }
}


/**
 * @brief Retrieves a double value associated with a given key.
 * 
 * This function retrieves the value corresponding to the specified key from the
 * configuration dictionary and converts it to a double. If the key does not exist
 * or if the value cannot be converted to a double, an exception is thrown.
 * 
 * @param key The key whose value is to be retrieved and converted to double.
 * @return The double value associated with the key.
 * 
 * @throws std::runtime_error If the key is not found or the value cannot be converted to double.
 */
double InputConfig::GetDouble(const std::string& key) 
{
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


/**
 * @brief Retrieves an integer value associated with a given key.
 * 
 * This function retrieves the value corresponding to the specified key from the
 * configuration dictionary and converts it to an integer. If the key does not exist
 * or if the value cannot be converted to an integer, an exception is thrown.
 * 
 * @param key The key whose value is to be retrieved and converted to int.
 * @return The integer value associated with the key.
 * 
 * @throws std::runtime_error If the key is not found or the value cannot be converted to int.
 */
int InputConfig::GetInt(const std::string& key) 
{
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


/**
 * @brief Retrieves a size_t value associated with a given key.
 * 
 * This function retrieves the value corresponding to the specified key from the
 * configuration dictionary and converts it to a `size_t`. If the key does not exist
 * or if the value cannot be converted to a `size_t`, an exception is thrown.
 * 
 * @param key The key whose value is to be retrieved and converted to `size_t`.
 * @return The `size_t` value associated with the key.
 * 
 * @throws std::runtime_error If the key is not found or the value cannot be converted to `size_t`.
 */
std::size_t InputConfig::GetSizeT(const std::string& key)
{
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


/**
 * @brief Retrieves a boolean value associated with a given key.
 * 
 * This function retrieves the value corresponding to the specified key from the
 * configuration dictionary and converts it to a boolean value. If the key does not exist,
 * an exception is thrown.
 * 
 * @param key The key whose value is to be retrieved and converted to bool.
 * @return The boolean value associated with the key.
 * 
 * @throws std::runtime_error If the key is not found in the dictionary.
 */
bool InputConfig::GetBool(const std::string& key) 
{
    
    if (Contains(key)) {
        return string_utils::ConvertStringToBool(dict_[key]);
    } else {
        throw std::runtime_error("Key not found: " + key); // Key not found
    }
}