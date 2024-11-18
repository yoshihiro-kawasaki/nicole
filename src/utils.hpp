#ifndef UTILITY_HPP_
#define UTILITY_HPP_

#include <iostream>
#include <string>
#include <vector>
#include <map>

namespace Utility 
{

std::vector<std::string> Split(const std::string &s, char delim);

// template <class T>
class Dictionary
{
public:
    // explicit operator bool() const { return valid; }
    bool Contains(const std::string &key) const { return dict.find(key) != dict.end(); };
    void Insert(std::string key, double item) { dict[key] = item; };

    size_t Size() { return dict.size(); }

    double Get_double(std::string key) 
    { 
        if (Contains(key)) {
            return  dict[key];
        } else {
            return 0.00;
        }
    }

private:
    std::map<std::string, double> dict;
    // bool valid;
};

class InputConfigure
{
public:

    InputConfigure();
    InputConfigure(std::string file_name);
    ~InputConfigure();

    bool Contains(const std::string &key);
    void Insert(std::string key, std::string item);
    void ReadFile(std::string file_name);
    size_t Size();

    std::string Get(std::string key);

private:

    std::map<std::string, std::string> dict;

};

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

#endif /* UTILITY_HPP_ */