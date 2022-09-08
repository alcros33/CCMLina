#pragma once
#include "factorization.hpp"
namespace ccm
{
// solves the linear system asuming an upper triangular matrix
void _solve_upper(const Matrix& U, const Matrix& b, Matrix& xout);
// solves the linear system asuming an upper triangular matrix
Matrix solve_upper(const Matrix& U, const Matrix& b);

// solves the linear system asuming a lower triangular matrix
void _solve_lower(const Matrix& U, const Matrix& b, Matrix& xout);
// solves the linear system asuming a lower triangular matrix
Matrix solve_lower(const Matrix& U, const Matrix& b);

// solves the linear system asuming a diagonal matrix stored as a row
// TODO CHANGE TO DiagMatrix
void _solve_diag(const Matrix& D, const Matrix& b, Matrix& xout);

// solves the linear system asuming a diagonal matrix stored as a row
// TODO CHANGE TO DiagMatrix
Matrix solve_diag(const Matrix& D, const Matrix& b);

// solves the linear using gaussian elimination
void _solve_gauss(const Matrix& A, const Matrix& b, Matrix& xout);

// solves the linear using gaussian elimination
Matrix solve_gauss(const Matrix& A, const Matrix& b);

// Solves the linear system using results from LU factorization
void _solve_lu(const Matrix& L, const Matrix& U, const Matrix& b, Matrix& xout);

// solves the linear using result cholesky factorization
void _solve_cholesky(Matrix& L, const Matrix& b, Matrix& xout);

} // namespace ccm
