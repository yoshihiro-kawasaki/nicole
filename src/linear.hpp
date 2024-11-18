#ifndef LINEAR_HPP_
#define LINEAR_HPP_

#include <array>
#include <cmath>
#include <memory>
#include <vector>

class Linear
{
private:
public:

    static bool AbsCompare(double a, double b);

    size_t IndexOfLargestElementInVector(const std::vector<double> &dx, const size_t n, const size_t offset);

    void ScalesVectorByConstant(const double da, std::vector<double> &dx, const size_t n, const size_t offset);

    double DotProduct(const std::vector<double> &a, const std::vector<double> &b, const size_t n,
                 const size_t offsetA, const size_t offsetB);

    void ConstantTimesVectorPlusVector(const double da, const std::vector<double> &dx, std::vector<double> &dy,
                const size_t n, const size_t offsetX,
                const size_t offsetY);

    void SolveByGaussianElimination(const std::vector<std::vector<double>> &a, const size_t n, std::vector<int> &ipvt,
               std::vector<double> &b, const size_t job);

    void FactorMatrixByGaussianElimination(std::vector<std::vector<double>> &a, const size_t n, std::vector<int> &ipvt,
               size_t *const info);
};

#endif // LINEAR_HPP_