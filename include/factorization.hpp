#pragma once
#include "matrix.hpp"
#include <tuple>

namespace ccm
{
// Return two Matrices L,U (lower and upper triangular) such that A=LU and Us diagonal are Ones
void _factor_crout(const Matrix& A, Matrix& Lout, Matrix& Uout);

// Return two Matrices L,U (lower and upper triangular) such that A=LU and Us diagonal are Ones
auto factor_crout(const Matrix& A) -> std::tuple<Matrix, Matrix>;

// Return two Matrices L,U (lower and upper triangular) such that A=LU and Ls diagonal are Ones
void _factor_dolittle(const Matrix& A, Matrix& Lout, Matrix& Uout);

// Return two Matrices L,U (lower and upper triangular) such that A=LU and Ls diagonal are Ones
auto factor_dolittle(const Matrix& A) -> std::tuple<Matrix, Matrix>;

// Return three Matrices L,D,U (lower triangular, diagonal, upper triangular) such that A=LDU and Ls
// and Us diagonal are Ones
void _factor_ldu(const Matrix& A, Matrix& Lout, Matrix& Dout, Matrix& Uout);

// Return three Matrices L,D,U (lower triangular, diagonal, upper triangular) such that A=LDU and Ls
// and Us diagonal are Ones
auto factor_ldu(const Matrix& A) -> std::tuple<Matrix, Matrix, Matrix>;

// Return a Matrix L (lower triangular) such that A=LL^t. A must be symetric
void _factor_cholesky(const Matrix& A, Matrix& Lout);

// Return a Matrix L (lower triangular) such that A=LL^t. A must be symetric
Matrix factor_cholesky(const Matrix& A);

} // namespace ccm