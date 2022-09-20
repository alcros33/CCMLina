#pragma once
#include "matrix.hpp"
namespace ccm
{
void _solve_gauss_seidel(
  const Matrix& A, const Matrix& b, Matrix& xout, double eps = 1E-5, usize_t maxiters = 1000);
void _solve_jacobi(
  const Matrix& A, const Matrix& b, Matrix& xout, double eps = 1E-5, usize_t maxiters = 1000);
}