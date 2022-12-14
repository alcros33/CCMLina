#pragma once
#include "matrix_trait.hpp"
#include "submatrix.hpp"
#include <vector>

namespace ccm
{

class Matrix : public MatrixTrait<Matrix, Matrix>
{
public:
    using parent = MatrixTrait<Matrix, Matrix>;
    using parent::operator=;

    Matrix(std::initializer_list<double> L);
    Matrix(std::initializer_list<std::initializer_list<double>> L);

    explicit Matrix();
    Matrix(size_t n, size_t m, double fill = 0);
    Matrix(Matrix&& other);
    void swap(Matrix&);

    static Matrix randu(size_t n, size_t m);
    static Matrix randu(size_t n, size_t m, double from, double to);
    static Matrix eye(size_t n);

    Matrix(const Matrix& other);

    template <class MatT>
    Matrix(const MatT& other) : parent()
    {
        m_nrow = other.rows(), m_ncol = other.cols();
        m_data = new double[m_nrow * m_ncol];
        this->operator=(other);
    }
    Matrix& operator=(const Matrix& rhs);
    Matrix& operator=(const SubMatrix& rhs);
    // Matrix& operator=(const DiagMatrix& rhs);
    // Matrix& operator=(const BandMatrix& rhs);
    Matrix& operator=(Matrix&& rhs);

    ~Matrix();

    // Resizes matrix, re-allocating if needed
    // Data is lost in the process
    Matrix& resize(size_t n, size_t m);

    // Iterator of all elements
    double* begin();
    // constant Iterator of all elements
    const double* cbegin() const;
    // Iterator of all elements
    double* end();
    // constant Iterator of all elements
    const double* cend() const;

    // Iterator of all non-zero elements
    double* fbegin();
    // constant Iterator of all non-zero elements
    const double* cfbegin() const;
    // Iterator of all non-zero elements
    double* fend();
    // constant Iterator of all non-zero elements
    const double* cfend() const;

    size_t row_col_to_idx(size_t r, size_t c) const;

    SubMatrix slice(size_t from_row, size_t to_row, size_t from_col, size_t to_col);
    const SubMatrix slice(size_t from_row, size_t to_row, size_t from_col, size_t to_col) const;

    SubMatrix row(size_t r);
    const SubMatrix row(size_t r) const;

    SubMatrix col(size_t c);
    const SubMatrix col(size_t c) const;

    // Transposes the matrix in-place
    Matrix& itranspose();

    // IO
    friend std::istream& operator>>(std::istream& is, Matrix& M);

protected:
    void _itranspose_cycle();
    void _itranspose_square();
};

} // namespace ccm