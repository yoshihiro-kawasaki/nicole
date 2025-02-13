/**
 * @file string_utils.cpp
 * @brief implimentation of string utils
 * @date 2025-02-11
 * @author Y. Kawasaki
 */

#include "string_utils.hpp"

namespace string_utils
{
    /**
     * @brief 文字列の前後の空白を削除する関数
     * @param str 処理対象の文字列
     * @return 空白が削除された文字列
     * @example
     * std::string str = "  hello world  \n";
     * std::string trimmed_str = string_utils::Trim(str); // "hello world"
     */
    std::string Trim(const std::string& str) {
        std::string whitespace = " \t\n\r";  // 削除対象の空白文字
        size_t first = str.find_first_not_of(whitespace);
        size_t last = str.find_last_not_of(whitespace);
        if (first == std::string::npos) return ""; // すべて空白の場合
        return str.substr(first, last - first + 1);
    }

    /**
     * @brief 文字列を特定の区切り文字で分割する関数
     * @param str 処理対象の文字列
     * @param delim 区切り文字
     * @param trimWhitespace 分割された文字列の空白をトリムするかどうか
     * @return 分割された文字列のベクトル
     * @example
     * std::string str = "apple,banana,orange";
     * std::vector<std::string> fruits = string_utils::Split(str, ",", true); // {"apple", "banana", "orange"}
     */
    std::vector<std::string> Split(const std::string &str, char delim, bool trimWhitespace) {

        std::vector<std::string> elems;
        std::stringstream ss(str);
        std::string item;

        while (std::getline(ss, item, delim)) {
            if (trimWhitespace) {
                item = Trim(item);
                if (!item.empty()) elems.push_back(item);   
            } else {
                elems.push_back(item);
            }
        }
        
        return elems;
    }

    /**
     * @brief 文字列の配列内に特定の文字列が存在するかを判定する関数
     * @param string_vector 検索対象の文字列ベクトル
     * @param target 検索する文字列
     * @return 存在する場合は true, 存在しない場合は false
     * @example
     * std::vector<std::string> fruits = {"apple", "banana", "orange"};
     * bool has_apple = string_utils::IsInStringVector(fruits, "apple"); // true
     * bool has_grape = string_utils::IsInStringVector(fruits, "grape"); // false
     */
    bool IsInStringVector(const std::vector<std::string>& string_vector, const std::string& target) {
        return std::find(string_vector.begin(), string_vector.end(), target) != string_vector.end();
    }

    /**
    * @brief 文字列をbool値に変換する関数
    * @param strbool 変換対象の文字列
    * @return 変換されたbool値
    * @throw std::invalid_argument str が "true" または "false" でない場合
    * @example
    * bool value = string_utils::ConvertStringToBool("  TRUE  "); // true
    */
    bool ConvertStringToBool(const std::string& strbool) {

        std::string str = Trim(strbool);
        std::string lower_str;
        std::transform(str.begin(), str.end(), std::back_inserter(lower_str), ::tolower); // 小文字に変換

        if (lower_str == "true") {
            return true;
        } else if (lower_str == "false") {
            return false;
        } else {
            throw std::invalid_argument("str must be 'true' or 'false'. str = " + strbool);
        }

    }
}