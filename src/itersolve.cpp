#include "itersolve.hpp"
namespace ccm
{
void _solve_gauss_seidel(
  const Matrix& A, const Matrix& b, Matrix& xout, double eps, size_t maxiters)
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
        for (size_t i = 0; i < n; i++)
        {
            buff_sum = 0;
            xp = x_new.data();
            Ap = A.data() + i * n; // Select ith row
            for (size_t j = 0; j < i; j++)
                buff_sum += Ap[j] * xp[j];
            for (size_t j = i + 1; j < n; j++)
                buff_sum += Ap[j] * xp[j];
            x_new(i, 0) = (b(i, 0) - buff_sum) / A(i, i);
        }
        rel_err = sqrt(xout.sdist(x_new) / x_new.snorm());
        xout.swap(x_new);
    } while ((itt++) < maxiters && rel_err > eps);
}

void _solve_jacobi(const Matrix& A, const Matrix& b, Matrix& xout, double eps, size_t maxiters)
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
        for (size_t i = 0; i < n; i++)
        {
            buff_sum = 0;
            xp = xout.data();
            Ap = A.data() + i * n; // Select ith row
            for (size_t j = 0; j < i; j++)
                buff_sum += Ap[j] * xp[j];
            for (size_t j = i + 1; j < n; j++)
                buff_sum += Ap[j] * xp[j];
            x_new(i, 0) = (b(i, 0) - buff_sum) / A(i, i);
        }
        rel_err = sqrt(xout.sdist(x_new) / x_new.snorm());
        xout.swap(x_new);
    } while ((itt++) < maxiters && rel_err > eps);
}

Matrix solve_grad(const Matrix& A, const Matrix& b, double eps)
{
    Matrix x(A.rows(), 1, 0);
    Matrix r = -b;
    Matrix P = -r;
    Matrix W = P;
    double alpha, beta;
    double rsnorm = r.snorm(), newrsnorm;
    while (sqrt(rsnorm) > eps)
    {
        A._matmul(P, W);
        alpha = rsnorm / P.dot(W);
        x += P * alpha;
        r += W * alpha;
        newrsnorm = r.snorm();
        beta = newrsnorm / rsnorm;
        rsnorm = newrsnorm;
        P *= beta;
        P -= r;
    }

    return x;
}

void diag_matmul(const Matrix& A, const Matrix& y, Matrix& xout)
{
    for (size_t i = 0; i < A.rows(); i++)
        xout(i, 0) = y(i, 0) * A(i, i);
}

Matrix solve_grad_precond(const Matrix& A, const Matrix& b, double eps)
{

    Matrix x(A.rows(), 1, 0);
    Matrix r = -b;
    Matrix Minv(A.rows(), A.cols());
    for (size_t i = 0; i < A.rows(); i++)
        Minv(i, i) = 1 / A(i, i);

    Matrix y = r;
    diag_matmul(Minv, r, y);
    Matrix P = -y, W = P;
    double alpha, beta;
    double rdoty = r.dot(y), newrdoty;
    while (true)
    {
        A._matmul(P, W);
        alpha = rdoty / P.dot(W);
        x += P * alpha;
        r += W * alpha;
        if (r.norm() < eps)
            break;
        diag_matmul(Minv, r, y);
        newrdoty = r.dot(y);
        beta = newrdoty / rdoty;
        rdoty = newrdoty;
        P *= beta;
        P -= y;
    }

    return x;
}

} // namespace ccm
