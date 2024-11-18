#include "utils.hpp"

#include <iostream>
#include <fstream>

namespace Utility {

std::vector<std::string> Split(const std::string &s, char delim)
{
    std::vector<std::string> elems;
    std::string item;
    for (char ch: s) {
        if (ch == delim) {
            if (!item.empty())
                elems.push_back(item);
            item.clear();
        }
        else {
            item += ch;
        }
    }
    if (!item.empty())
        elems.push_back(item);
    return elems;
}

// class InputConfigure

InputConfigure::InputConfigure()
{

}

InputConfigure::InputConfigure(std::string file_name)
{
    std::ifstream file(file_name, std::ios::in);
    if (!file) {
        std::cout << "error open file : " << file_name << std::endl;
        exit(1);
    }

    std::string str;
    std::vector<std::string> line;

    while (std::getline(file, str)) {
        std::cout << str[0] << std::endl;
        if (str[0] == '#') continue;
        
        line = Split(str, ' ');
        dict[line[0]] = line[2];
    } 

    file.close(); 
}

InputConfigure::~InputConfigure()
{

}

bool InputConfigure::Contains(const std::string &key)
{
    return dict.find(key) != dict.end();
}

void InputConfigure::Insert(std::string key, std::string item)
{
    dict[key] = item;
}

void InputConfigure::ReadFile(std::string file_name)
{
    std::ifstream file(file_name, std::ios::in);
    if (!file) {
        std::cout << "error open file : " << file_name << std::endl;
        exit(1);
    }

    std::string str;
    std::vector<std::string> line;

    while (std::getline(file, str)) {

        if (str[0] == '#') continue;
        
        line = Split(str, ' ');
        dict[line[0]] = line[2];
    } 

    file.close();

    return;
}

size_t InputConfigure::Size()
{
    return dict.size();
}

std::string InputConfigure::Get(std::string key)
{
    if (Contains(key)) {
        return dict[key];
    } else {
        std::cout << "## No item" << std::endl;
        return "";
    }
}

// template <>
// bool dictionary::get<bool>(std::string key, bool value) {
//     if (!contains(key)) {
//         return value;
//     }
//     if (("yes" == dict[key]) || ("Yes" == dict[key])) {
//         return true;
//     } else {
//         return false;
//     }
// }

// template <>
// bool dictionary::get(std::string key) {
//     return get<bool>(key, false);
// }

// template <>
// int dictionary::get(std::string key, int value) {
//     if (!contains(key)) {
//         return value;
//     }
//     return stoi(dict[key]);
// }

// template <>
// int dictionary::get(std::string key) {
//     return get<int>(key, 0);
// }

// template <>
// double dictionary::get(std::string key, double value) {
//     if (!contains(key)) {
//         return value;
//     }
//     return stod(dict[key]);
// }

// template <>
// double dictionary::get(std::string key) {
//     return get<double>(key, 0.0);
// }

}