#include "linear.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <numeric>
#include <vector>

bool Linear::AbsCompare(double a, double b)
{
    return (std::abs(a) < std::abs(b));
}


/* Purpose : Find largest component of double vector dx */
size_t Linear::IndexOfLargestElementInVector(const std::vector<double> &dx, const size_t n, const size_t offset=0)
{
    double v = 0, vmax = 0;
    size_t idmax = 0;
    size_t i;

    for(i = 0; i < n; i++) {
        v = std::abs(dx[i+offset]);
        if( v > vmax ) {
            vmax = v;
            idmax = i;
        }
    }

    return idmax;
}

void Linear::ScalesVectorByConstant(const double da, std::vector<double> &dx, const size_t n, const size_t offset)
{
    std::transform( dx.begin()+offset, dx.end(), dx.begin()+offset, [&da](double x) -> double { return da*x; } );
}

double Linear::DotProduct(const std::vector<double> &a, const std::vector<double> &b, const size_t n,
                 const size_t offsetA = 0, const size_t offsetB = 0)
{
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) sum += a[i+offsetA] * b[i+offsetB];
    return sum;
}


void Linear::ConstantTimesVectorPlusVector(const double da, const std::vector<double> &dx, std::vector<double> &dy,
                const size_t n, const size_t offsetX = 0, const size_t offsetY = 0)
{
    for (size_t i = 0; i < n; i++) {
        dy[i + offsetY] = da * dx[i + offsetX] + dy[i + offsetY];
    }
}


void Linear::SolveByGaussianElimination(const std::vector<std::vector<double> > &a, const size_t n, std::vector<int> &ipvt,
               std::vector<double> &b, const size_t job)
{
    int k, j;
    double t;

    /*
       Job = 0, solve a * x = b.
    */
    if (job == 0)
    {
        /*
           First solve L * y = b.
        */
        for (k = 0; k < n; k++)
        {
            t = DotProduct(a[k], b, k);
            b[k] = (b[k] - t) / a[k][k];
        }

        /*
           Now solve U * x = y.
        */
        for (k = n - 2; k >= 0; k--) {

            b[k] = b[k] + DotProduct(a[k], b, n - k - 1, k + 1, k + 1);
            // std::cout << "b = " <<  b[k] << std::endl;
            j = ipvt[k];
            // std::cout << "j = " << j << std::endl;
            if (j != k)
            {
                t = b[j];
                b[j] = b[k];
                b[k] = t;
            }

        }

        return;
    }
    /*
       Job = nonzero, solve Transpose(a) * x = b.

       First solve Transpose(U) * y = b.
    */
    for (k = 0; k < n - 1; k++) {

        j = ipvt[k];
        t = b[j];
        if (j != k)
        {
            b[j] = b[k];
            b[k] = t;
        }
        ConstantTimesVectorPlusVector(t, a[k], b, n - k, k, k);
        
    }
    /*
       Now solve Transpose(L) * x = y.
    */
    for (k = n; k >= 1; k--)
    {
        b[k] = b[k] / a[k][k];
        t = -b[k];
        ConstantTimesVectorPlusVector(t, a[k], b, k - 1);
    }

}

void Linear::FactorMatrixByGaussianElimination(std::vector<std::vector<double> > &a, const size_t n, std::vector<int> &ipvt,
                  size_t *const info)
{
    size_t j = 0, k = 0, i = 0;
    double t = 0.0;

    /* Gaussian elimination with partial pivoting.   */

    *info = 0;
    for (k = 0; k < n - 1; k++) {

        /*
           Find j = pivot index.  Note that a[k]+k-1 is the address of
           the 0-th element of the row vector whose 1st element is a[k][k].
        */
        j = IndexOfLargestElementInVector(a[k], n - k, k) + k;
        ipvt[k] = j;
        // std::cout << "j = " << j << std::endl;

        /*
           Zero pivot implies this row already triangularized.
        */
        if (a[k][j] == 0.0) {
            *info = k;
            continue;
        }

        /*
           Interchange if necessary.
        */

        if (j != k)
        {
            t = a[k][j];
            a[k][j] = a[k][k];
            a[k][k] = t;
        }

        /*
           Compute multipliers.
        */
        t = -1.0 / a[k][k];
        ScalesVectorByConstant(t, a[k], n - k - 1, k + 1);
        // ScalesVectorByConstant(t, a[k], )
        /*
           Column elimination with row indexing.
        */
        for (i = k + 1; i < n; i++)
        {
            t = a[i][j];
            if (j != k)
            {
                a[i][j] = a[i][k];
                a[i][k] = t;
            }
            ConstantTimesVectorPlusVector(t, a[k], a[i], n - k - 1, k + 1, k + 1);
        }

    } /* end k-loop  */

    ipvt[n-1] = n-1;
    if (a[n-1][n-1] == 0.0) *info = n-1;

}