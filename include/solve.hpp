#pragma once
#include "asolve.hpp"
#include "itersolve.hpp"

namespace ccm
{
enum class SolveMethod
{
    Gauss,
    Crout,
    Dolittle,
    LDU,
    Cholesky,
    GaussSeidel,
    Jacobi,
};

// Solves the linear system Ax=b
void _solve(const Matrix& A, const Matrix& b, Matrix& xout, SolveMethod method);
// Solves the linear system Ax=b
Matrix solve(const Matrix& A, const Matrix& b, SolveMethod method);

} // namespace ccm
