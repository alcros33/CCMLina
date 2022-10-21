#include "functions.hpp"

namespace ccm
{

Matrix linspace(double from, double to, size_t N, bool endpoint)
{
    Matrix R(N, 1);
    double eps = endpoint ? (to - from) / (N - 1) : (to - from) / N;
    for (size_t i = 0; i < N; i++)
        R(i, 0) = from + i * eps;
    return R;
}

} // namespace ccm