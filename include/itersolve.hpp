#pragma once
#include "matrix.hpp"
namespace ccm
{
void _gauss_seidel_solve(const Matrix& A, const Matrix& b, Matrix& xout);
void _jacobi_solve(const Matrix& A, const Matrix& b, Matrix& xout);
}