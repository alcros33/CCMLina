#include "asolve.hpp"

namespace ccm
{

void _solve_upper(const Matrix& U, const Matrix& b, Matrix& xout)
{
    CCM_ASSERT(U.is_square(), "Only for square matrices");
    CCM_ASSERT((U.rows() == b.rows()), "Mismatched sizes");
    CCM_ASSERT(b.is_vector(), "b is not a column vector");
    CCM_ASSERT_SAME_SIZE(b, xout);
    auto n = U.rows();
    double* xptr;
    const double* Uptr;
    double sum;
    for (size_t i = n - 1; i >= 0; i--)
    {
        if (U(i, i) == 0)
        {
            std::cout << "Singular Matrix" << std::endl;
            return;
        }
        sum = 0;
        Uptr = U.data() + i * n; // select ith row
        xptr = xout.data(); // al vector data
        for (usize_t j = i + 1; j < n; j++)
        {
            sum += Uptr[j] * xptr[j];
        }
        xout(i, 0) = (b(i, 0) - sum) / U(i, i);
    }
}

// solves the linear system asuming a lower triangular matrix
void _solve_lower(const Matrix& L, const Matrix& b, Matrix& xout)
{
    CCM_ASSERT(L.is_square(), "Only for square matrices");
    CCM_ASSERT((L.rows() == b.rows()), "Mismatched sizes");
    CCM_ASSERT(b.is_vector(), "b is not a column vector");
    CCM_ASSERT_SAME_SIZE(b, xout);
    auto n = L.rows();
    double* xptr;
    const double* Lptr;
    double sum;
    for (usize_t i = 0; i < n; i++)
    {
        if (L(i, i) == 0)
        {
            std::cout << "Singular Matrix" << std::endl;
            return;
        }
        sum = 0;
        Lptr = L.data() + i * n; // select ith row
        xptr = xout.data(); // All vector data
        for (usize_t j = 0; j < i; j++)
        {
            sum += Lptr[j] * xptr[j];
        }
        xout(i, 0) = (b(i, 0) - sum) / L(i, i);
    }
}

void _solve_diag(const Matrix& D, const Matrix& b, Matrix& xout)
{
    CCM_ASSERT(D.is_square(), "Only for square matrices");
    CCM_ASSERT((D.rows() == b.rows()), "Mismatched sizes");
    CCM_ASSERT(b.is_vector(), "b is not a column vector");
    CCM_ASSERT_SAME_SIZE(b, xout);
    for (size_t i = 0; i < D.rows(); i++)
    {
        if (D(i, i) == 0)
        {
            std::cout << "Singular Matrix" << std::endl;
            return;
        };
        xout(i, 0) = b(i, 0) / D(i, i);
    }
}

void _solve_gauss(const Matrix& A, const Matrix& b, Matrix& xout)
{
    CCM_ASSERT(A.is_square(), "Only for square matrices");
    CCM_ASSERT((A.rows() == b.rows()), "Mismatched sizes");
    CCM_ASSERT(b.is_vector(), "b is not a column vector");
    CCM_ASSERT_SAME_SIZE(b, xout);

    // Gaussian Elimination
    ccm::Matrix Am(A), bm(b);
    ccm::Matrix row_tmp(1, Am.cols());
    auto n = Am.rows();
    ccm::usize_t piv_idx;
    double factor;
    double piv_val;
    for (ccm::usize_t i = 0; i < n - 1; i++)
    {
        // Row pivoting
        piv_idx = i;
        piv_val = fabs(Am(i, i));
        for (ccm::usize_t row = i + 1; row < n; row++)
        {
            if (piv_val < fabs(Am(row, i)))
            {
                piv_idx = row;
                piv_val = fabs(Am(row, i));
            }
        }
        if (piv_idx != i)
        {
            row_tmp = Am.slice(i, i + 1, 0, n);
            Am.slice(i, i + 1, 0, n) = Am.slice(piv_idx, piv_idx + 1, 0, n);
            Am.slice(piv_idx, piv_idx + 1, 0, n) = row_tmp;
            std::swap(bm(i, 0), bm(piv_idx, 0));
        }

        for (ccm::usize_t row = i + 1; row < n; row++)
        {
            if (Am(i, i) == 0)
            {
                std::cout << "Singular Matrix" << std::endl;
                return;
            }
            factor = Am(row, i) / Am(i, i);
            Am.slice(row, row + 1, i, n) -= Am.slice(i, i + 1, i, n) * factor;
            bm(row, 0) -= bm(i, 0) * factor;
        }
    }

    _solve_upper(Am, bm, xout);
}

void _solve_lu(const Matrix& L, const Matrix& U, const Matrix& b, Matrix& xout)
{
    auto y = solve_lower(L, b);
    _solve_upper(U, y, xout);
}

void _solve_cholesky(Matrix& L, const Matrix& b, Matrix& xout)
{
    auto y = solve_lower(L, b);
    _solve_upper(L.itranspose(), y, xout);
    L.itranspose();
}

Matrix solve_upper(const Matrix& U, const Matrix& b)
{
    Matrix x(U.rows(), 1);
    _solve_upper(U, b, x);
    return x;
}

Matrix solve_lower(const Matrix& L, const Matrix& b)
{
    Matrix x(L.rows(), 1);
    _solve_lower(L, b, x);
    return x;
}

Matrix solve_diag(const Matrix& D, const Matrix& b)
{
    Matrix x(D.rows(), 1);
    _solve_diag(D, b, x);
    return x;
}

Matrix solve_gauss(const Matrix& A, const Matrix& b)
{
    Matrix x(A.rows(), 1);
    _solve_gauss(A, b, x);
    return x;
}

Matrix solve_QR(const Matrix& Q, const Matrix& R, const Matrix& b)
{
    auto qtb = Q.transpose().matmul(b);
    return solve_upper(R, qtb);
}

} // namespace ccm