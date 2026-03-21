#ifndef INPUT_CONFIG_HPP_
#define INPUT_CONFIG_HPP_

// C++ headers
#include <map>
#include <string>

#include "nicole/nicole_defs.hpp"

namespace nicole {
    class InputConfig {
    public:
        InputConfig(const std::string& filename);

        void ReadFile(const std::string& filename);
        bool Contains(const std::string& key) const;
        void Insert(const std::string& key, const std::string& item);

        // Getter
        const std::string& GetFileName() const { return input_filename_; }
        const std::string& GetString(const std::string& key) const;
        Real GetReal(const std::string& key) const;
        int GetInt(const std::string& key) const;
        std::size_t GetSizeT(const std::string& key) const;
        bool GetBool(const std::string& key) const;
    private:
        std::map<std::string, std::string> dict_;
        std::string input_filename_;
    };
}

#endif /* INPUT_CONFIG_HPP_ */
