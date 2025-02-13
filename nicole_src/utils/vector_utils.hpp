/**
 * @file vector_utils.hpp
 * @brief Utility functions for vector
 * @date 2025-02-11
 * @author Y. Kawasaki
 */
#ifndef VECTOR_UTILS_HPP
#define VECTOR_UTILS_HPP

#include <vector>

namespace vector_utils
{
    /**
     * @brief ベクトル内の最大値を取得
     * @tparam T ベクトルの要素型
     * @param vec ベクトル
     * @return ベクトル内の最大値
     * @throw std::runtime_error ベクトルが空の場合
     * @example
     * std::vector<int> v = {1, 3, 2};
     * int max_value = vector_utils::GetMaxValue(v); // 3
     */
    template <typename T>
    T GetMaxValue(const std::vector<T>& vec) {
        if (vec.empty()) {
            throw std::runtime_error("Vector is empty.");
        }
        return *std::max_element(vec.begin(), vec.end());
    }

    /**
     * @brief ベクトル内のindiceが示すインデックスの値における最大値を取得
     * @tparam T ベクトルの要素型
     * @param vec 検索対象のベクトル
     * @param indice インデックスの配列
     * @param index_start 検索開始インデックス
     * @param index_end 検索終了インデックス
     * @return 指定された範囲内の最大値
     * @throw std::runtime_error ベクトルまたはインデックス配列が空の場合
     * @throw std::out_of_range インデックス範囲が無効な場合、またはインデックス値がベクトルの範囲外の場合
     * @example
     * std::vector<int> v = {1, 3, 2, 5, 4};
     * std::vector<int> indices = {0, 2, 4};
     * int max_value = vector_utils::GetMaxValue(v, indices, 0, 2); // 5
     */
    template <typename T>
    T GetMaxValue(const std::vector<T>& vec, const std::vector<std::size_t>& indice, const std::size_t index_start, const std::size_t index_end) {

        if (vec.empty() || indice.empty()) {
            throw std::runtime_error("Vector or index array is empty.");
        }

        if (index_start < 0 || index_start > index_end) {
            throw std::out_of_range("Invalid index range.");
        }

        if (index_end >= indice.size()) {
            throw std::out_of_range("index_end is greater than indice size.");
        }

        std::size_t vec_size = vec.size();
        T max_value = vec[indice[index_start]];  // 最初の要素を基準に最大値を初期化

        for (std::size_t i = index_start; i <= index_end; ++i) {
            if (indice[i] < 0 || indice[i] >= vec_size) {
                throw std::out_of_range("Index out of bounds.");
            }
            max_value = std::max(max_value, vec[indice[i]]);
        }

        return max_value;
    }

    /**
     * @brief ベクトル内の最大値のインデックスを取得
     * @tparam T ベクトルの要素型
     * @param vec 検索対象のベクトル
     * @return ベクトル内の最大値のインデックス
     * @throw std::runtime_error ベクトルが空の場合
     * @example
     * std::vector<int> v = {1, 3, 2};
     * int max_index = vector_utils::GetMaxValueIndex(v); // 1
     */
    template <typename T>
    int GetMaxValueIndex(const std::vector<T>& vec) {
        if (vec.empty()) {
            throw std::runtime_error("Vector is empty.");
        }
        return std::distance(vec.begin(), std::max_element(vec.begin(), vec.end()));
    }

    /**
     * @fn GetMaxAbsValue
     * @brief ベクトル内の絶対値の最大値を取得
     */
    template <typename T>
    T GetMaxAbsValue(const std::vector<T>& vec) {
        if (vec.empty()) {
            throw std::runtime_error("Vector is empty.");
        }
        return *std::max_element(vec.begin(), vec.end(), [](T a, T b) { return std::abs(a) < std::abs(b); });
    }

    /**
     * @fn GetMaxAbsValueIndex
     * @brief ベクトル内の絶対値の最大値のインデックスを取得
     */
    template <typename T>
    int GetMaxAbsValueIndex(const std::vector<T>& vec) {
        if (vec.empty()) {
            throw std::runtime_error("Vector is empty.");
        }
        return std::distance(vec.begin(), std::max_element(vec.begin(), vec.end(), [](T a, T b) { return std::abs(a) < std::abs(b); }));
    }

    /**
     * @fn GetMinValue
     * @brief ベクトル内の最小値の取得
     */
    template <typename T>
    T GetMinValue(const std::vector<T>& vec) {
        if (vec.empty()) {
            throw std::runtime_error("Vector is empty.");
        }
        return *std::min_element(vec.begin(), vec.end());
    }

    /**
     * @class GetMinValue
     * @brief ベクトル内のindiceが示すインデックスの値における最小値を取得
     */
    template <typename T>
    T GetMinValue(const std::vector<T>& vec, const std::vector<int>& indice, const std::size_t index_start, const std::size_t index_end) 
    {
        if (vec.empty() || indice.empty()) {
            throw std::runtime_error("Vector or index array is empty.");
        }

        if (index_start < 0 || index_end >= indice.size() || index_start > index_end) {
            throw std::out_of_range("Invalid index range.");
        }

        T min_value = vec[indice[index_start]];  // 最初の要素を基準に最大値を初期化
        int vec_size = static_cast<int>(vec.size());

        for (int i = index_start; i <= index_end; ++i) {
            if (indice[i] < 0 || indice[i] >= vec_size) {
                throw std::out_of_range("Index out of bounds.");
            }
            min_value = std::min(min_value, vec[indice[i]]);
        }

        return min_value;
    }

    /**
     * @fn GetMinValueIndex
     * @brief ベクトル内の最小値のインデックスの取得
     */
    template <typename T>
    int GetMinValueIndex(const std::vector<T>& vec) {
        if (vec.empty()) {
            throw std::runtime_error("Vector is empty.");
        }
        return std::distance(vec.begin(), std::min_element(vec.begin(), vec.end()));
    }

    /**
     * @fn GetMinAbsValue
     * @brief ベクトル内の絶対値の最小値を取得
     */
    template <typename T>
    T GetMinAbsValue(const std::vector<T>& vec) {
        if (vec.empty()) {
            throw std::runtime_error("Vector is empty.");
        }
        
        return *std::min_element(vec.begin(), vec.end(), [](T a, T b) { return std::abs(a) < std::abs(b); });
    }

    /**
     * @fn GetMinAbsValueIndex
     * @brief ベクトル内の絶対値の最小値のインデックスを取得
     */
    template <typename T>
    int GetMinAbsValueIndex(const std::vector<T>& vec) {
        if (vec.empty()) {
            throw std::runtime_error("Vector is empty.");
        }
        return std::distance(vec.begin(), std::min_element(vec.begin(), vec.end(), [](T a, T b) { return std::abs(a) < std::abs(b); }));
    }


}

#endif /* VECTOR_UTILS_HPP */