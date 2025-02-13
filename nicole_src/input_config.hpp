/**
 * @file input_config.hpp
 * @brief Manages input configuration by reading key-value pairs from a file.
 * 
 * This class provides functionality to read configuration files and retrieve values
 * as strings, integers, doubles, or booleans.
 * 
 * @date 2025-02-11
 * @author Y. Kawasaki
 */

#ifndef INPUT_CONFIG_HPP
#define INPUT_CONFIG_HPP

#include "utils/string_utils.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <algorithm>
#include <stdexcept>

/**
 * @class InputConfig
 * @brief A class for handling input configuration files.
 */
class InputConfig
{
public:

    /**
     * @brief Constructor that reads the configuration file.
     * @param filename Path to the configuration file.
     */
    InputConfig(const std::string& filename);

    /**
     * @brief Get the name of the input file.
     * @return The input file name as a string.
     */
    const std::string& GetFileName() { return input_filename_; }

    /**
     * @brief Reads a configuration file and stores key-value pairs.
     * @param filename Path to the configuration file.
     */
    void ReadFile(const std::string& filename);

    /**
     * @brief Checks if the given key exists in the configuration.
     * @param key The key to search for.
     * @return True if the key exists, false otherwise.
     */
    bool Contains(const std::string& key) const;

    /**
     * @brief Inserts a new key-value pair into the configuration.
     * @param key The key string.
     * @param item The corresponding value string.
     */
    void Insert(const std::string& key, const std::string& item);

    /**
     * @brief Retrieves the value associated with a key as a string.
     * @param key The key to look up.
     * @return The value as a string.
     * @throws std::runtime_error if the key does not exist.
     */
    const std::string& GetString(const std::string& key);

    /**
     * @brief Retrieves the value associated with a key as a double.
     * @param key The key to look up.
     * @return The value as a double.
     * @throws std::runtime_error if the key does not exist or is not a valid number.
     */
    double GetDouble(const std::string& key);

    /**
     * @brief Retrieves the value associated with a key as an integer.
     * @param key The key to look up.
     * @return The value as an integer.
     * @throws std::runtime_error if the key does not exist or is not a valid integer.
     */
    int GetInt(const std::string& key);

    /**
     * @brief Retrieves the value associated with a key as an std::size_t.
     * @param key The key to look up.
     * @return The value as an std::size_t.
     * @throws std::runtime_error if the key does not exist or is not a valid std::size_t
     */
    std::size_t GetSizeT(const std::string& key);

    /**
     * @brief Retrieves the value associated with a key as a boolean.
     * @param key The key to look up.
     * @return True if the value is "true" or "1", false otherwise.
     * @throws std::runtime_error if the key does not exist or is not a valid boolean.
     */
    bool GetBool(const std::string& key);

private:

    // Stores key-value pairs read from the configuration file.
    std::map<std::string, std::string> dict_;

    // The name of the input configuration file.
    std::string input_filename_;

};

#endif /* INPUT_CONFIG_HPP */