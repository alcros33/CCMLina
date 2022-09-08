#include "submatrix.hpp"
#include "matrix.hpp"
namespace ccm
{

SubMatrix::SubMatrix(Matrix& Mat, usize_t offrow, usize_t offcol, usize_t rows, usize_t cols)
    : parent()
    , m_const(false)
    , m_offset_row(offrow)
    , m_offset_col(offcol)
    , m_og_cols(Mat.cols())
    , m_cdata(Mat.data())
{
    m_data = Mat.data();
    m_nrow = rows;
    m_ncol = cols;
}

SubMatrix::SubMatrix(const Matrix& Mat, usize_t offrow, usize_t offcol, usize_t rows, usize_t cols)
    : parent()
    , m_const(true)
    , m_offset_row(offrow)
    , m_offset_col(offcol)
    , m_og_cols(Mat.cols())
    , m_cdata(Mat.data())
{
    m_data = nullptr;
    m_nrow = rows;
    m_ncol = cols;
}

usize_t SubMatrix::row_col_to_idx(size_t r, size_t c) const
{
    return (r + m_offset_row) * m_og_cols + c + m_offset_col;
}
SliceIter SubMatrix::begin()
{
    CCM_ASSERT(!m_const, "Trying to get non-const iterator from const SubMatrix");
    return SliceIter(this);
}
SliceIter SubMatrix::end()
{
    CCM_ASSERT(!m_const, "Trying to get non-const iterator from const SubMatrix");
    SliceIter res(this);
    res.m_idx_col = 0, res.m_idx_row = m_nrow;
    return res;
}
ConstSliceIter SubMatrix::cbegin() const { return ConstSliceIter(this); }
ConstSliceIter SubMatrix::cend() const
{
    ConstSliceIter res(this);
    res.m_idx_col = 0, res.m_idx_row = m_nrow;
    return res;
}
SliceIter SubMatrix::fbegin() { return begin(); }
SliceIter SubMatrix::fend() { return end(); }
ConstSliceIter SubMatrix::cfbegin() const { return cbegin(); }
ConstSliceIter SubMatrix::cfend() const { return cend(); }

double* SubMatrix::data()
{
    CCM_ASSERT(!m_const, "Trying to get non-const pointer from const SubMatrix");
    return m_data;
}
const double* SubMatrix::data() const { return m_cdata; }

#define SUBMATRIX_ASSIGN_IMPL(OtherMat)                                                            \
    SubMatrix& SubMatrix::operator=(const OtherMat& rhs)                                           \
    {                                                                                              \
        CCM_ASSERT_SAME_SIZE((*this), rhs);                                                        \
        CCM_ASSERT(!m_const, "Trying to assign on const SubMatrix");                               \
        auto it1 = begin();                                                                        \
        auto it2 = rhs.cbegin();                                                                   \
        for (size_t i = 0; i < m_ncol * m_nrow; i++)                                               \
        {                                                                                          \
            it1[i] = it2[i];                                                                       \
        }                                                                                          \
        return *this;                                                                              \
    }

SUBMATRIX_ASSIGN_IMPL(SubMatrix)
SUBMATRIX_ASSIGN_IMPL(Matrix)

namespace inner
{
    template <bool Const, class UT>
    implSliceIterator<Const, UT>::implSliceIterator(implSliceIterator::mat_ptr M) : m_M(M)
    {
        m_data = m_M->data() + m_M->row_col_to_idx(m_idx_row, m_idx_col);
    }

    template <bool Const, class UT>
    UT& implSliceIterator<Const, UT>::underlying()
    {
        return static_cast<UT&>(*this);
    }
    template <bool Const, class UT>
    UT& implSliceIterator<Const, UT>::operator+=(int incr)
    {
        auto dr = lldiv(m_idx_col + incr, m_M->cols());
        m_idx_row += dr.quot;
        m_idx_col = dr.rem;
        m_data = m_M->data() + m_M->row_col_to_idx(m_idx_row, m_idx_col);
        return underlying();
    }
    template <bool Const, class UT>
    UT implSliceIterator<Const, UT>::operator+(int incr)
    {
        UT res(underlying());
        res += incr;
        return res;
    }
    template <bool Const, class UT>
    UT implSliceIterator<Const, UT>::operator++(int)
    {
        UT res(underlying());
        this->operator++();
        return res;
    }
    template <bool Const, class UT>
    UT& implSliceIterator<Const, UT>::operator++()
    {
        if (m_idx_col == m_M->cols() - 1)
        {
            m_idx_row++;
            m_idx_col = 0;
            m_data = m_M->data() + m_M->row_col_to_idx(m_idx_row, m_idx_col);
        }
        else
        {
            m_idx_col++;
            m_data++;
        }
        return underlying();
    }
    template <bool Const, class UT>
    auto implSliceIterator<Const, UT>::operator-(const implSliceIterator& other) const ->
      typename implSliceIterator<Const, UT>::difference_type
    {
        auto tot = m_idx_row * m_M->cols() + m_idx_col;
        auto otot = other.m_idx_row * other.m_M->cols() + other.m_idx_col;
        return tot - otot;
    }
    template <bool Const, class UT>
    auto implSliceIterator<Const, UT>::operator[](usize_t incr) ->
      typename implSliceIterator<Const, UT>::reference
    {
        auto dr = lldiv(m_idx_col + incr, m_M->cols());
        return m_M->data()[m_M->row_col_to_idx(m_idx_row + dr.quot, dr.rem)];
    }
    template <bool Const, class UT>
    auto implSliceIterator<Const, UT>::operator*() ->
      typename implSliceIterator<Const, UT>::reference
    {
        return *m_data;
    }
    template <bool Const, class UT>
    typename implSliceIterator<Const, UT>::pointer implSliceIterator<Const, UT>::operator->()
    {
        return m_data;
    }
    template <bool Const, class UT>
    bool implSliceIterator<Const, UT>::operator==(const implSliceIterator& b)
    {
        return m_idx_row == b.m_idx_row && m_idx_col == b.m_idx_col;
    }
    template <bool Const, class UT>
    bool implSliceIterator<Const, UT>::operator!=(const implSliceIterator& b)
    {
        return !(this->operator==(b));
    }
    template class implSliceIterator<true, ConstSliceIter>;
    template class implSliceIterator<false, SliceIter>;
} // namespace inner

} // namespace ccm