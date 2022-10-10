#pragma once
#include "matrix_trait.hpp"
namespace ccm
{

namespace inner
{
    template <bool Const, class UT>
    class implSliceIterator
    {
    public:
        using iterator_category = std::random_access_iterator_tag;
        using difference_type = size_t;
        using reference = typename std::conditional_t<Const, const double&, double&>;
        using pointer = typename std::conditional_t<Const, const double*, double*>;
        using value_type = typename std::decay<reference>::type;
        using mat_ptr = typename std::conditional_t<Const, const SubMatrix*, SubMatrix*>;

        implSliceIterator() {}

        implSliceIterator(mat_ptr M);

        UT& underlying();
        UT& operator+=(int incr);
        UT operator+(int incr);
        UT operator++(int);
        UT& operator++();

        difference_type operator-(const implSliceIterator& other) const;

        reference operator[](size_t incr);

        reference operator*();
        pointer operator->();

        bool operator==(const implSliceIterator& b);
        bool operator!=(const implSliceIterator& b);

    protected:
        mat_ptr m_M;
        size_t m_idx_row{0}, m_idx_col{0};
        pointer m_data;
    };
} // namespace inner

class ConstSliceIter : public inner::implSliceIterator<true, ConstSliceIter>
{
    using parent = inner::implSliceIterator<true, ConstSliceIter>;
    using parent::parent;
    friend class SubMatrix;
};
class SliceIter : public inner::implSliceIterator<false, SliceIter>
{
    using parent = inner::implSliceIterator<false, SliceIter>;
    using parent::parent;
    friend class SubMatrix;
};

// A reference view into a matrix
class SubMatrix : public MatrixTrait<SubMatrix, Matrix>
{
public:
    using parent = MatrixTrait<SubMatrix, Matrix>;
    using parent::operator=;

    SubMatrix() = delete;
    SubMatrix(const SubMatrix& other) = default;
    SubMatrix(SubMatrix&& other) = default;

    SubMatrix(Matrix& Mat, size_t offrow, size_t offcol, size_t rows, size_t cols);
    SubMatrix(const Matrix& Mat, size_t offrow, size_t offcol, size_t rows, size_t cols);

    SubMatrix& operator=(const SubMatrix& other);
    SubMatrix& operator=(const Matrix& other);

    size_t row_col_to_idx(size_t r, size_t c) const;

    double* data();
    const double* data() const;

    friend class ConstSliceIter;
    friend class SliceIter;

    // Iterator of all elements
    SliceIter begin();
    // Iterator of all elements
    SliceIter end();
    // constant Iterator of all elements
    ConstSliceIter cbegin() const;
    // constant Iterator of all elements
    ConstSliceIter cend() const;

    // Iterator of all non-zero elements
    SliceIter fbegin();
    // Iterator of all non-zero elements
    SliceIter fend();
    // constant Iterator of all non-zero elements
    ConstSliceIter cfbegin() const;
    // constant Iterator of all non-zero elements
    ConstSliceIter cfend() const;

protected:
    bool m_const;
    size_t m_offset_row, m_offset_col;
    size_t m_og_cols;
    const double* m_cdata;
};

} // namespace ccm