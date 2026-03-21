#include <algorithm>
#include <sstream>

#include "nicole/utils/string_utils.hpp"

namespace string_utils {
    std::string Trim(const std::string& str) {
        const std::string trim_characters = " \t\n\r";
        const std::size_t first = str.find_first_not_of(trim_characters);
        const std::size_t last = str.find_last_not_of(trim_characters);
        if (first == std::string::npos) return "";
        return str.substr(first, last - first + 1);
    }


    std::vector<std::string> Split(const std::string& str, const char delim, const bool trim_white_space) {
        std::vector<std::string> elems;
        std::stringstream ss(str);
        std::string item;

        while (std::getline(ss, item, delim)) {
            if (trim_white_space) item = Trim(item);
            if (!item.empty()) elems.push_back(item);
        }

        return elems;
    }


    bool IsInStringVector(const std::vector<std::string>& string_vector, const std::string& target) {
        return std::find(string_vector.begin(), string_vector.end(), target) != string_vector.end();
    }

    
    bool ConvertStringToBool(const std::string& strbool) {
        std::string str = Trim(strbool);
        std::string lower_str;
        std::transform(str.begin(), str.end(), std::back_insert_iterator(lower_str), ::tolower);

        if (lower_str == "true") {
            return true;
        } else if (lower_str == "false") {
            return false;
        } else {
            throw std::invalid_argument("str must be 'true' or 'false'. str = " + strbool);
        }
    }
}
