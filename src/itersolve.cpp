#include "itersolve.hpp"
namespace ccm
{
void _solve_gauss_seidel(
  const Matrix& A, const Matrix& b, Matrix& xout, double eps, usize_t maxiters)
{
    int itt = 0;
    auto n = A.rows();
    double rel_err, buff_sum;
    Matrix x_new(xout);
    // Pointers to iterate over x and A
    double* xp;
    const double* Ap;
    do
    {
        for (usize_t i = 0; i < n; i++)
        {
            buff_sum = 0;
            xp = x_new.data();
            Ap = A.data() + i * n; // Select ith row
            for (usize_t j = 0; j < i; j++)
                buff_sum += Ap[j] * xp[j];
            for (usize_t j = i + 1; j < n; j++)
                buff_sum += Ap[j] * xp[j];
            x_new(i, 0) = (b(i, 0) - buff_sum) / A(i, i);
        }
        rel_err = sqrt(xout.sdist(x_new) / x_new.snorm());
        xout.swap(x_new);
    } while ((itt++) < maxiters && rel_err > eps);
}

void _solve_jacobi(const Matrix& A, const Matrix& b, Matrix& xout, double eps, usize_t maxiters)
{
    int itt = 0;
    auto n = A.rows();
    double rel_err, buff_sum;
    Matrix x_new(xout);
    // Pointers to iterate over x and A
    double* xp;
    const double* Ap;
    do
    {
#pragma omp parallel for private(buff_sum, xp, Ap)
        for (usize_t i = 0; i < n; i++)
        {
            buff_sum = 0;
            xp = xout.data();
            Ap = A.data() + i * n; // Select ith row
            for (usize_t j = 0; j < i; j++)
                buff_sum += Ap[j] * xp[j];
            for (usize_t j = i + 1; j < n; j++)
                buff_sum += Ap[j] * xp[j];
            x_new(i, 0) = (b(i, 0) - buff_sum) / A(i, i);
        }
        rel_err = sqrt(xout.sdist(x_new) / x_new.snorm());
        xout.swap(x_new);
    } while ((itt++) < maxiters && rel_err > eps);
}
} // namespace ccm
