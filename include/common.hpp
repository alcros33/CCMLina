#pragma once
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <type_traits>

namespace ccm
{
using size_t = std::size_t;
using ssize_t = std::make_signed_t<size_t>;

class Matrix;
class SubMatrix;
class DiagMatrix;
class BandMatrix;

#ifndef NDEBUG
#define CCM_ASSERT(expr, msg)                                                                      \
    if (!expr)                                                                                     \
    {                                                                                              \
        std::cerr << msg << " in " << __PRETTY_FUNCTION__ << std::endl;                            \
        exit(1);                                                                                   \
    }
#else
#define CCM_ASSERT(expr, msg) ((void)0)
#endif

#define CCM_ASSERT_SAME_SIZE(A, B)                                                                 \
    CCM_ASSERT(A.same_size(B), "Mismatched sizes " << A.size() << " " << B.size())

} // namespace ccm