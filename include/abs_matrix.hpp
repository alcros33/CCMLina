#pragma once
#include "common.hpp"

namespace ccm
{

#define CCM_ABS_EWISE_OP_IMPL(OP)                                                                  \
    template <class Other>                                                                         \
    EvalRetType operator OP(const AbsMatrix<Other, EvalRetType>& rhs) const                        \
    {                                                                                              \
        CCM_ASSERT_SAME_SIZE((*this), rhs);                                                        \
        EvalRetType res(underlying());                                                             \
        res OP## = rhs;                                                                            \
        return res;                                                                                \
    }                                                                                              \
    template <class Other>                                                                         \
    UnderlyingT& operator OP##=(const AbsMatrix<Other, EvalRetType>& rhs)                          \
    {                                                                                              \
        CCM_ASSERT_SAME_SIZE((*this), rhs);                                                        \
        auto it1 = underlying().begin();                                                           \
        auto it2 = rhs.underlying().cbegin();                                                      \
        auto count = underlying().end() - it1;                                                     \
        for (size_t i = 0; i < count; i++)                                                         \
        {                                                                                          \
            it1[i] OP## = it2[i];                                                                  \
        }                                                                                          \
        return underlying();                                                                       \
    }

#define CCM_ABS_EWISE_SINGLE_OP_IMPL(OP)                                                           \
    EvalRetType operator OP(double rhs) const                                                      \
    {                                                                                              \
        EvalRetType res(underlying());                                                             \
        res OP## = rhs;                                                                            \
        return res;                                                                                \
    }                                                                                              \
    UnderlyingT& operator OP##=(double rhs)                                                        \
    {                                                                                              \
        for (auto it = underlying().begin(); it != underlying().end(); it++)                       \
        {                                                                                          \
            (*it) OP## = rhs;                                                                      \
        }                                                                                          \
        return underlying();                                                                       \
    }                                                                                              \
    friend EvalRetType operator OP(double lhs, const AbsMatrix& rhs)                               \
    {                                                                                              \
        EvalRetType res(rhs.underlying());                                                         \
        for (auto it = res.begin(); it != res.end(); it++)                                         \
        {                                                                                          \
            (*it) = lhs OP(*it);                                                                   \
        }                                                                                          \
        return res;                                                                                \
    }

template <class UnderlyingT, class EvalRetType>
class AbsMatrix
{
public:
    explicit AbsMatrix() {}
    AbsMatrix(const AbsMatrix& other) = delete;
    AbsMatrix(AbsMatrix&& other) = delete;
    AbsMatrix& operator=(const AbsMatrix& other) = delete;
    UnderlyingT& operator=(double rhs)
    {
        for (auto it = underlying().fbegin(); it != underlying().fend(); it++)
        {
            (*it) = rhs;
        }
        return underlying();
    }

    UnderlyingT& underlying() { return static_cast<UnderlyingT&>(*this); }
    const UnderlyingT& underlying() const { return static_cast<const UnderlyingT&>(*this); }

    usize_t rows() const { return m_nrow; }
    usize_t cols() const { return m_ncol; }

    bool is_square() const { return m_nrow == m_ncol; }
    bool is_vector() const { return m_ncol == 1; }

    // access to raw pointer
    double* data() { return m_data; }
    const double* data() const { return m_data; }

    // Squared L2 Norm of Matrix
    double snorm() const
    {
        double res = 0;
        for (auto it = underlying().cfbegin(); it != underlying().cfend(); it++)
        {
            res += (*it) * (*it);
        }
        return res;
    }
    // L2 Norm of matrix
    double norm() const { return sqrt(snorm()); }

    // Divides all elements by the norm, so that the norm becomes 1
    UnderlyingT& normalize()
    {
        double _norm = this->norm();
        CCM_ASSERT((_norm != 0), "Matrix is 0, cannot normalize");
        this->operator/=(_norm);
        return underlying();
    }

    // Squared euclidean distance to other Matrix of same size
    template <class Other>
    double sdist(const AbsMatrix<Other, EvalRetType>& B) const
    {
        CCM_ASSERT_SAME_SIZE((*this), B);
        double res = 0;
        double tmp;
        auto it1 = underlying().cbegin();
        auto it2 = B.underlying().cbegin();
        for (usize_t i = 0; i < m_ncol * m_nrow; i++)
        {
            tmp = it1[i] - it2[i];
            res += tmp * tmp;
        }
        return res;
    }

    // For tests and debugging
    template <class Other>
    bool operator==(const AbsMatrix<Other, EvalRetType>& B) const
    {
        CCM_ASSERT_SAME_SIZE((*this), B);
        auto it1 = underlying().cbegin();
        auto it2 = B.underlying().cbegin();
        for (usize_t i = 0; i < m_ncol * m_nrow; i++)
        {
            if (it1[i] != it2[i])
                return false;
        }
        return true;
    }

    // Euclidean distance to other Matrix of same size
    template <class Other>
    double dist(const AbsMatrix<Other, EvalRetType>& B) const
    {
        return sqrt(sdist<Other>(B));
    }

    double sum() const
    {
        double res = 0;
        for (auto it = underlying().cfbegin(); it != underlying().cfend(); it++)
        {
            res += (*it);
        }
        return res;
    }
    double mean() const { return sum() / (m_nrow * m_ncol); }
    double item() const { return underlying().data()[0]; }

    // element access
    double& operator()(size_t r, size_t c)
    {
        assert(0 <= r && r < m_nrow && 0 <= c && c < m_ncol);
        return underlying().data()[underlying().row_col_to_idx(r, c)];
    }
    // element access
    const double& operator()(size_t r, size_t c) const
    {
        assert(0 <= r && r < m_nrow && 0 <= c && c < m_ncol);
        return underlying().data()[underlying().row_col_to_idx(r, c)];
    }

    EvalRetType transpose() const
    {
        EvalRetType res(underlying());
        res.itranspose();
        return res;
    }
    template <class Other>
    bool same_size(const AbsMatrix<Other, EvalRetType>& B) const
    {
        return m_ncol == B.cols() && m_nrow == B.rows();
    }

    // Operators

    CCM_ABS_EWISE_OP_IMPL(+)
    CCM_ABS_EWISE_OP_IMPL(-)
    CCM_ABS_EWISE_OP_IMPL(*)
    CCM_ABS_EWISE_OP_IMPL(/)

    CCM_ABS_EWISE_SINGLE_OP_IMPL(+)
    CCM_ABS_EWISE_SINGLE_OP_IMPL(-)
    CCM_ABS_EWISE_SINGLE_OP_IMPL(*)
    CCM_ABS_EWISE_SINGLE_OP_IMPL(/)

    template <class Other>
    EvalRetType matmul(const AbsMatrix<Other, EvalRetType>& B) const
    {
        CCM_ASSERT((m_ncol == B.rows()), "Mismatched dimensions for matmul");
        EvalRetType res(m_nrow, B.cols(), 0);
        _matmul(B, res);
        return res;
    }

    // Stores matmul result in res output parameter
    template <class Other>
    void _matmul(const AbsMatrix<Other, EvalRetType>& B, EvalRetType& ResOut) const
    {
        // Transpose B first because iterating over cols is faster
        auto BT = B.transpose();
        double *Rp, *Bp;
        auto Ap = underlying().cbegin();
#pragma omp parallel for collapse(2) private(Rp, Ap, Bp)
        for (usize_t i = 0; i < ResOut.rows(); i++)
        {
            for (usize_t j = 0; j < ResOut.cols(); j++)
            {
                Rp = ResOut.data() + i * ResOut.cols() + j; // select ith row jth column from result
                Ap = underlying().cbegin() + i * m_ncol; // Select ith row from A
                Bp = BT.data() + j * BT.cols(); // select jth row from BTransposed
                for (usize_t k = 0; k < m_ncol; k++)
                {
                    (*Rp) += Ap[k] * Bp[k];
                }
            }
        }
    }

protected:
    usize_t m_nrow{0}, m_ncol{0};
    double* m_data{nullptr};
};

template <class UnderlyingT, class EvalRetType>
std::ostream& operator<<(std::ostream& os, const AbsMatrix<UnderlyingT, EvalRetType>& M)
{
    auto old_prec = os.precision();
    os << M.rows() << " " << M.cols() << "\n";
    os << std::setprecision(12);
    for (usize_t i = 0; i < M.rows(); i++)
    {
        for (usize_t j = 0; j < M.cols(); j++)
        {
            os << M(i, j);
            if (j != (M.cols() - 1))
                os << " ";
        }
        if (i != (M.rows() - 1))
            os << "\n";
    }
    os << std::setprecision(old_prec);
    return os;
}
} // namespace ccm
