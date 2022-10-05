#include "matrix.hpp"
#include "solve.hpp"
#include <tuple>
namespace ccm
{

// Computes an eigen value and its vector given initial values x0 and mu
std::tuple<double, Matrix> eig_rayleigh(const Matrix& A, Matrix& x0, double mu0, double eps = 1E-7);
// Computes the n biggest eigen values (and asociated vectors).
// Requires symetric matrix A=At
std::tuple<Matrix, Matrix> eig_power(const Matrix& A, usize_t n, double eps = 1E-5);
// Computes the n smallest eigen values (and asociated vectors).
// Requires symetric matrix A=At
std::tuple<Matrix, Matrix>
eig_inv_power(const Matrix& A, usize_t n, SolveMethod method, double eps = 1E-5);
// Computes all the eigen values (and asociated vectors).
// Requires symetric matrix A=At
std::tuple<Matrix, Matrix> eig_jacobi(const Matrix& A, double eps = 1E-5);
// Subspace iteration method for m biggest eigen values (and asociated vectors).
// Requires symetric matrix A=At
std::tuple<Matrix, Matrix> eig_subspace_pow(
  const Matrix& A, usize_t m, int iters, double eps = 1E-7, double jacobi_eps = 1E-3);
// Subspace iteration method for m smallest eigen values (and asociated vectors).
// Requires symetric matrix A=At
std::tuple<Matrix, Matrix> eig_subspace_inv(
  const Matrix& A, usize_t m, int inv_iters, double eps = 1E-7, double jacobi_eps = 1E-3);
} // namespace ccm
