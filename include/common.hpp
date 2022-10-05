#pragma once
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <vector>
namespace ccm
{
using size_t = long long;
using usize_t = unsigned long long;

class Matrix;
class SubMatrix;
class DiagMatrix;
class BandMatrix;

#ifndef NDEBUG
#define CCM_ASSERT(expr, msg)                                                                      \
    if (!expr)                                                                                     \
    {                                                                                              \
        std::cerr << msg << " at " << __FILE__ << ":" << __LINE__ << " in " << __PRETTY_FUNCTION__ \
                  << std::endl;                                                                    \
        exit(-1);                                                                                  \
    }
#else
#define CCM_ASSERT(expr, msg) ((void)0)
#endif

#define CCM_ASSERT_SAME_SIZE(first, second) CCM_ASSERT(first.same_size(second), "Mismatched sizes")

} // namespace ccm