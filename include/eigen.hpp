#include "matrix.hpp"
#include "solve.hpp"
#include <tuple>
namespace ccm
{
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
} // namespace ccm
