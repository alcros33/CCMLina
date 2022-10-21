#pragma once
#include "matrix.hpp"

namespace ccm
{
Matrix linspace(double from, double to, size_t N = 1000, bool endpoint = true);
} // namespace ccm
