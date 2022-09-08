#include "factorization.hpp"

namespace ccm
{

void _factor_dolittle(const Matrix& A, Matrix& Lout, Matrix& Uout)
{
    CCM_ASSERT(A.is_square(), "Only for square matrices");
    CCM_ASSERT_SAME_SIZE(A, Lout);
    CCM_ASSERT_SAME_SIZE(A, Uout);
    for (ccm::usize_t i = 0; i < A.rows(); i++)
    {
        for (ccm::usize_t j = 0; j < i; j++)
        {
            Uout(j, i) = A(j, i);
            for (ccm::usize_t k = 0; k < j; k++)
            {
                Uout(j, i) -= Lout(j, k) * Uout(k, i);
            }
        }
        for (ccm::usize_t j = 0; j < i; j++)
        {
            Lout(i, j) = A(i, j);
            for (ccm::usize_t k = 0; k < j; k++)
            {
                Lout(i, j) -= Lout(i, k) * Uout(k, j);
            }
            Lout(i, j) /= Uout(j, j);
        }
        Uout(i, i) = A(i, i);
        for (ccm::usize_t k = 0; k < i; k++)
        {
            Uout(i, i) -= Lout(i, k) * Uout(k, i);
        }
        if (Uout(i, i) == 0)
        {
            std::cout << "Singular Matrix" << std::endl;
            break;
        }
    }
}

void _factor_crout(const Matrix& A, Matrix& Lout, Matrix& Uout)
{
    CCM_ASSERT(A.is_square(), "Only for square matrices");
    CCM_ASSERT_SAME_SIZE(A, Lout);
    CCM_ASSERT_SAME_SIZE(A, Uout);

    for (ccm::usize_t i = 0; i < A.rows(); i++)
    {
        for (ccm::usize_t j = 0; j < i; j++)
        {
            Lout(i, j) = A(i, j);
            for (ccm::usize_t k = 0; k < j; k++)
            {
                Lout(i, j) -= Lout(i, k) * Uout(k, j);
            }
        }
        for (ccm::usize_t j = 0; j < i; j++)
        {
            Uout(j, i) = A(j, i);
            for (ccm::usize_t k = 0; k < j; k++)
            {
                Uout(j, i) -= Lout(j, k) * Uout(k, i);
            }
            Uout(j, i) /= Lout(j, j);
        }
        Lout(i, i) = A(i, i);
        for (ccm::usize_t k = 0; k < i; k++)
        {
            Lout(i, i) -= Lout(i, k) * Uout(k, i);
        }
        if (Lout(i, i) == 0)
        {
            std::cout << "Singular Matrix" << std::endl;
            break;
        }
    }
}

void _factor_ldu(const Matrix& A, Matrix& Lout, Matrix& Dout, Matrix& Uout)
{
    CCM_ASSERT(A.is_square(), "Only for square matrices");
    CCM_ASSERT_SAME_SIZE(A, Lout);
    CCM_ASSERT_SAME_SIZE(A, Dout);
    CCM_ASSERT_SAME_SIZE(A, Uout);

    for (ccm::usize_t i = 0; i < A.rows(); i++)
    {
        // Row
        for (ccm::usize_t j = 0; j < i; j++)
        {
            Lout(i, j) = A(i, j);
            for (ccm::usize_t k = 0; k < j; k++)
            {
                Lout(i, j) -= Lout(i, k) * Dout(k, k) * Uout(k, j);
            }
            Lout(i, j) /= Dout(j, j);
        }
        // Column
        for (ccm::usize_t j = 0; j < i; j++)
        {
            Uout(j, i) = A(j, i);
            for (ccm::usize_t k = 0; k < j; k++)
            {
                Uout(j, i) -= Lout(j, k) * Dout(k, k) * Uout(k, i);
            }
            Uout(j, i) /= Dout(j, j);
        }
        // Diagonal
        Dout(i, i) = A(i, i);
        for (ccm::usize_t k = 0; k < i; k++)
        {
            Dout(i, i) -= Lout(i, k) * Dout(k, k) * Uout(k, i);
        }
        if (Dout(i, i) == 0)
        {
            std::cout << "Singular Matrix" << std::endl;
            break;
        }
    }
}

void _factor_cholesky(const Matrix& A, Matrix& Lout)
{
    CCM_ASSERT(A.is_square(), "Only for square matrices");
    CCM_ASSERT_SAME_SIZE(A, Lout);
    auto n = A.rows();
    double *Lp, *LTp;
    double sum = 0;

    for (usize_t i = 0; i < n; i++)
    {
        for (usize_t j = 0; j <= i; j++)
        {

            sum = 0;
            Lp = Lout.begin() + i * n; // Select ith row
            LTp = Lout.begin() + j * n; // Seletct jth row
            for (usize_t k = 0; k < j; k++) // iter over columns until jth column
            {
                sum += Lp[k] * LTp[k];
            }

            if (i == j)
                Lout(i, i) = sqrt(A(i, i) - sum);
            else
                Lout(i, j) = (1.0 / Lout(j, j) * (A(i, j) - sum));
        }
    }
}

auto factor_crout(const Matrix& A) -> std::tuple<Matrix, Matrix>
{
    Matrix L(A.rows(), A.cols()), U(A.rows(), A.cols());
    _factor_crout(A, L, U);
    return std::make_tuple(L, U);
}

auto factor_dolittle(const Matrix& A) -> std::tuple<Matrix, Matrix>
{
    Matrix L(A.rows(), A.cols()), U(A.rows(), A.cols());
    _factor_dolittle(A, L, U);
    return std::make_tuple(L, U);
}

auto factor_ldu(const Matrix& A) -> std::tuple<Matrix, Matrix, Matrix>
{
    Matrix L(A.rows(), A.cols()), D(A.rows(), A.cols()), U(A.rows(), A.cols());
    _factor_ldu(A, L, D, U);
    return std::make_tuple(L, D, U);
}

Matrix factor_cholesky(const Matrix& A)
{
    Matrix L(A.rows(), A.cols());
    _factor_cholesky(A, L);
    return L;
}

} // namespace ccm