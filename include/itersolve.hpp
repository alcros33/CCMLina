#pragma once
#include "matrix.hpp"
namespace ccm
{
void _solve_gauss_seidel(
  const Matrix& A, const Matrix& b, Matrix& xout, double eps = 1E-5, usize_t maxiters = 1000);
void _solve_jacobi(
  const Matrix& A, const Matrix& b, Matrix& xout, double eps = 1E-5, usize_t maxiters = 1000);
// gradient conjugate method
Matrix solve_grad(const Matrix& A, const Matrix& b, double eps = 1E-7);
// gradient conjugate method, preconditioned with jacobi
Matrix solve_grad_precond(const Matrix& A, const Matrix& b, double eps = 1E-7);
}