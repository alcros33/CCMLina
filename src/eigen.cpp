#include "eigen.hpp"

namespace ccm
{

std::tuple<double, Matrix> eig_rayleigh(const Matrix& A, Matrix& x0, double mu0, double eps)
{
    auto LP = std::make_tuple(mu0, Matrix(x0));
    auto& lam = std::get<0>(LP);
    auto& phi = std::get<1>(LP);
    Matrix AmuI = A;
    double lam_old;
    phi.normalize();
    do
    {
        for (size_t i = 0; i < A.rows(); i++)
            AmuI(i, i) = A(i, i) - lam;
        x0 = solve_grad_precond(AmuI, phi);
        lam_old = lam;
        lam += 1 / phi.dot(x0);
        x0.normalize();
        phi.swap(x0);
    } while (std::abs(lam - lam_old) > eps);
    if (phi(0, 0) < 0)
        phi *= -1;
    return LP;
}

std::tuple<Matrix, Matrix> eig_power(const Matrix& A, size_t n, double eps)
{
    CCM_ASSERT((A.is_square()), "Only for square matrices");

    Matrix lam(n, 1, 0);
    std::vector<Matrix> VS;
    std::vector<double> ws(n);
    VS.reserve(n);

    auto sz = A.rows();
    Matrix tmp(sz, 1);
    double last_lam;

    for (size_t i = 0; i < n; i++)
    {
        VS.emplace_back(Matrix::randu(sz, 1));
        VS[i].normalize();
        tmp = VS[i];
        do
        {
            last_lam = lam(i, 0);
            for (size_t v_idx = 0; v_idx < i; v_idx++)
                ws[v_idx] = tmp.dot(VS[v_idx]);
            for (size_t v_idx = 0; v_idx < i; v_idx++)
                tmp -= (ws[v_idx] * VS[v_idx]);
            tmp.normalize();
            A._matmul(tmp, VS[i]);
            lam(i, 0) = VS[i].dot(tmp);
            VS[i].normalize();
            tmp.swap(VS[i]);
        } while (std::abs(lam(i, 0) - last_lam) > eps);
    }

    Matrix phis(sz, n);
    for (size_t i = 0; i < n; i++)
        phis.slice(0, sz, i, i + 1) = VS[i];

    return std::make_tuple(lam, phis);
}

std::tuple<Matrix, Matrix> eig_inv_power(const Matrix& A, size_t n, SolveMethod method, double eps)
{
    CCM_ASSERT((A.is_square()), "Only for square matrices");

    Matrix lam(n, 1, 0);
    std::vector<Matrix> VS;
    std::vector<double> ws(n);
    VS.reserve(n);
    ccm::Matrix L, U;
    switch (method)
    {
    case SolveMethod::Crout:
        L.resize(n, n);
        U.resize(n, n);
        _factor_crout(A, L, U);
        break;
    case SolveMethod::Dolittle:
        L.resize(n, n);
        U.resize(n, n);
        _factor_dolittle(A, L, U);
        break;
    case SolveMethod::Cholesky:
        L.resize(n, n);
        _factor_cholesky(A, L);
        break;
    default:
        break;
    }

    auto sz = A.rows();
    Matrix tmp(sz, 1);
    double last_lam;

    for (size_t i = 0; i < n; i++)
    {
        VS.emplace_back(Matrix::randu(sz, 1));
        VS[i].normalize();
        tmp = VS[i];
        do
        {
            last_lam = lam(i, 0);
            for (size_t v_idx = 0; v_idx < i; v_idx++)
                ws[v_idx] = tmp.dot(VS[v_idx]);
            for (size_t v_idx = 0; v_idx < i; v_idx++)
                tmp -= (ws[v_idx] * VS[v_idx]);
            tmp.normalize();

            switch (method)
            {
            case SolveMethod::Crout:
            case SolveMethod::Dolittle:
                _solve_lu(L, U, tmp, VS[i]);
                break;
            case SolveMethod::Cholesky:
                _solve_cholesky(L, tmp, VS[i]);
                break;
            case SolveMethod::Gauss:
                _solve_gauss(A, tmp, VS[i]);
                break;
            case SolveMethod::GaussSeidel:
                _solve_gauss_seidel(A, tmp, VS[i], eps);
                break;
            case SolveMethod::Jacobi:
                _solve_jacobi(A, tmp, VS[i], eps);
                break;
            default:
                break;
            }

            lam(i, 0) = 1.0 / VS[i].dot(tmp);
            VS[i].normalize();
            tmp.swap(VS[i]);
        } while (std::abs(lam(i, 0) - last_lam) > eps);
    }

    Matrix phis(sz, n);
    for (size_t i = 0; i < n; i++)
        phis.slice(0, sz, i, i + 1) = VS[i];

    return std::make_tuple(lam, phis);
}

std::tuple<size_t, size_t> matrix_max_offdiag(const Matrix& A)
{
    std::tuple<size_t, size_t> idx{0, 1};
    double m_v = 0;
    for (size_t i = 0; i < A.rows(); i++)
    {
        for (size_t j = 0; j < i; j++)
        {
            if (m_v < std::abs(A(i, j)))
            {
                idx = {i, j};
                m_v = std::abs(A(i, j));
            }
        }
    }

    return idx;
}

std::tuple<Matrix, Matrix> eig_jacobi(const ccm::Matrix& A, double eps)
{
    CCM_ASSERT((A.is_square()), "Only for square matrices");
    auto n = A.rows();
    Matrix Lam(A);
    auto Phi = Matrix::eye(n);
    Matrix Lam_i(n, 1), Lam_j(n, 1);
    Matrix Phi_i(n, 1);

    ccm::size_t i, j;
    std::tie(i, j) = matrix_max_offdiag(A);
    double theta, sv, cv;
    while (std::abs(Lam(i, j)) > eps)
    {
        theta = 0.5 * std::atan2(2 * Lam(i, j), Lam(j, j) - Lam(i, i));
        sv = std::sin(theta);
        cv = std::cos(theta);

        for (size_t it = 0; it < n; it++)
        {
            Lam_i(it, 0) = Lam(i, it);
            Lam_j(it, 0) = Lam(j, it);
            Phi_i(it, 0) = Phi(i, it);
        }

        Lam.col(i) = Lam_i * cv - Lam_j * sv;
        Lam.col(j) = Lam_j * cv + Lam_i * sv;
        for (ccm::size_t k = 0; k < n; k++)
        {
            Lam(j, k) = Lam(k, j);
            Lam(i, k) = Lam(k, i);
        }

        Lam(i, i) = cv * cv * Lam_i(i, 0) - 2 * sv * cv * Lam_i(j, 0) + sv * sv * Lam_j(j, 0);
        Lam(j, j) = sv * sv * Lam_i(i, 0) + 2 * sv * cv * Lam_i(j, 0) + cv * cv * Lam_j(j, 0);
        Lam(i, j) = Lam(j, i) = 0;

        Phi.col(i) = Phi_i * cv - Phi.col(j) * sv;
        Phi.col(j) = Phi.col(j) * cv + Phi_i * sv;

        std::tie(i, j) = matrix_max_offdiag(Lam);
    }
    return std::make_tuple(Lam, Phi);
}

std::tuple<Matrix, Matrix>
eig_subspace_pow(const Matrix& A, size_t m, int iters, double eps, double jacobi_eps)
{
    auto LP = std::make_tuple(Matrix(m, 1), Matrix::randu(A.rows(), m));
    CCM_ASSERT((A.is_square()), "Only for square matrices");
    CCM_ASSERT((A.rows() > m), "m must be less than the size of a");
    auto n = A.rows();

    auto& Phi = std::get<1>(LP);
    auto B = Matrix::eye(m);
    Matrix APhi(m, n);
    Matrix Z(m, m);
    Matrix Phitmp = Phi;
    Matrix tmp(n, 1), tmpnew(1, n);

    std::vector<double> ws(m);
    size_t i, j;
    while (true)
    {
        // power method for each Phi_i
        for (size_t i = 0; i < m; i++)
        {
            tmp = Phi.col(i);
            tmp.normalize();
            for (int itts = 0; itts < iters; itts++)
            {
                for (size_t v_idx = 0; v_idx < i; v_idx++)
                {
                    ws[v_idx] = tmp.dot(Phi.col(v_idx));
                    tmp -= (ws[v_idx] * Phi.col(v_idx));
                }
                tmp.itranspose()._tmatmul(A, tmpnew);
                tmp.swap(tmpnew);
                tmp.itranspose().normalize();
            }
            Phi.col(i) = tmp;
        }

        Phi.itranspose();
        Phi._tmatmul(A, APhi); // Since A is symetric
        APhi._tmatmul(Phi, B); // Phi is still transposed
        Phi.itranspose();

        std::tie(i, j) = matrix_max_offdiag(B);
        if (std::abs(B(i, j)) < eps)
            break;

        Z = std::get<1>(eig_jacobi(B, jacobi_eps));
        Phi._matmul(Z, Phitmp);
        Phi.swap(Phitmp);
    }
    auto& Lam = std::get<0>(LP);
    for (size_t i = 0; i < m; i++)
    {
        if (Phi(0, i) < 0)
            Phi.col(i) *= -1;
        Lam(i, 0) = B(i, i);
    }
    return LP;
}

std::tuple<Matrix, Matrix>
eig_subspace_inv(const Matrix& A, size_t m, int inv_iters, double eps, double jacobi_eps)
{
    auto LP = std::make_tuple(Matrix(m, 1), Matrix::randu(A.rows(), m));
    CCM_ASSERT((A.is_square()), "Only for square matrices");
    CCM_ASSERT((A.rows() > m), "m must be less than the size of a");
    auto n = A.rows();

    auto& Phi = std::get<1>(LP);
    auto L = factor_cholesky(A);
    Matrix APhi(m, n);
    Matrix Z(m, m);
    Matrix tmp(n, 1);
    Matrix Phitmp = Phi;

    auto B = Matrix::eye(m);
    Matrix B_old(B);

    std::vector<double> ws(m);
    size_t i, j;
    while (true)
    {
        // Inverse power method for each Phi_i
        for (size_t i = 0; i < m; i++)
        {
            tmp = Phi.col(i);
            tmp.normalize();

            for (int itts = 0; itts < inv_iters; itts++)
            {
                for (size_t v_idx = 0; v_idx < i; v_idx++)
                {
                    ws[v_idx] = tmp.dot(Phi.col(v_idx));
                    tmp -= (ws[v_idx] * Phi.col(v_idx));
                }
                _solve_cholesky(L, tmp, tmp);
                tmp.normalize();
            }
            Phi.col(i) = tmp;
        }

        Phi.itranspose();
        Phi._tmatmul(A, APhi); // Since A is symetric
        APhi._tmatmul(Phi, B); // Phi is still transposed
        Phi.itranspose();

        std::tie(i, j) = matrix_max_offdiag(B);
        if (std::abs(B(i, j)) < eps)
            break;

        Z = std::get<1>(eig_jacobi(B, jacobi_eps));

        Phi._matmul(Z, Phitmp);
        Phi.swap(Phitmp);
    }
    auto& Lam = std::get<0>(LP);
    for (size_t i = 0; i < m; i++)
    {
        if (Phi(0, i) < 0)
            Phi.col(i) *= -1;
        Lam(i, 0) = B(i, i);
    }

    return LP;
}

} // namespace ccm