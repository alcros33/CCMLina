#include "matrix.hpp"
#include "submatrix.hpp"
#include <random>
#include <vector>

namespace ccm
{

Matrix::Matrix() {}

Matrix::Matrix(usize_t n, usize_t m, double fill)
{
    m_nrow = n, m_ncol = m;
    m_data = new double[n * m];
    this->operator=(fill);
}

Matrix::Matrix(std::initializer_list<std::initializer_list<double>> L) : parent()
{
    m_nrow = L.size();
    m_ncol = L.begin()->size();
    m_data = new double[m_nrow * m_ncol];
    auto curr = L.begin();
    auto curr_iter = curr->begin();
    for (usize_t i = 0; i < m_nrow * m_ncol; i++)
    {
        m_data[i] = *curr_iter;
        if (++curr_iter == curr->end())
        {
            ++curr;
            curr_iter = curr->begin();
        }
    }
}

Matrix::Matrix(std::initializer_list<double> L) : parent()
{
    m_nrow = L.size();
    m_ncol = 1;
    m_data = new double[m_nrow];
    auto it = L.begin();
    for (usize_t i = 0; i < m_nrow; i++)
    {
        m_data[i] = it[i];
    }
}

Matrix::Matrix(const Matrix& other) : parent()
{
    m_nrow = other.m_nrow, m_ncol = other.m_ncol;
    m_data = new double[other.m_nrow * other.m_ncol];
    this->operator=(other);
}

#define MATRIX_ASSIGN_IMPL(OtherMat)                                                               \
    Matrix& Matrix::operator=(const OtherMat& rhs)                                                 \
    {                                                                                              \
        CCM_ASSERT_SAME_SIZE((*this), rhs);                                                        \
        auto it1 = begin();                                                                        \
        auto it2 = rhs.cbegin();                                                                   \
        for (size_t i = 0; i < m_ncol * m_nrow; i++)                                               \
        {                                                                                          \
            it1[i] = it2[i];                                                                       \
        }                                                                                          \
        return *this;                                                                              \
    }

MATRIX_ASSIGN_IMPL(Matrix)
MATRIX_ASSIGN_IMPL(SubMatrix)
// ASSIGN_MAT_IMPL(DiagMatrix)
// ASSIGN_MAT_IMPL(BandMatrix)

Matrix::Matrix(Matrix&& other)
{
    m_nrow = other.m_nrow, m_ncol = other.m_ncol;
    std::swap(m_data, other.m_data);
}

void Matrix::swap(Matrix& b)
{
    std::swap(m_data, b.m_data);
    std::swap(m_ncol, b.m_ncol);
    std::swap(m_nrow, b.m_nrow);
}

Matrix::~Matrix() { delete[] m_data; }

Matrix& Matrix::resize(usize_t n, usize_t m)
{
    if (m * n > m_nrow * m_ncol)
    {
        delete[] m_data;
        m_data = new double[m * n];
    }
    m_nrow = n, m_ncol = m;
    return *this;
}

Matrix Matrix::randomu(usize_t n, usize_t m)
{
    Matrix R(n, m);
    std::random_device rd; // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);
    for (auto it = R.begin(); it != R.end(); it++)
    {
        (*it) = dis(gen);
    }

    return R;
}

Matrix Matrix::eye(usize_t n)
{
    Matrix R(n, n);
    for (usize_t i = 0; i < n; i++)
        R(i, i) = 1;
    return R;
}

// iterators
double* Matrix::begin() { return m_data; }
const double* Matrix::cbegin() const { return m_data; }
double* Matrix::end() { return m_data + m_ncol * m_nrow; }
const double* Matrix::cend() const { return m_data + m_ncol * m_nrow; }
double* Matrix::fbegin() { return m_data; }
const double* Matrix::cfbegin() const { return m_data; }
double* Matrix::fend() { return m_data + m_ncol * m_nrow; }
const double* Matrix::cfend() const { return m_data + m_ncol * m_nrow; }

usize_t Matrix::row_col_to_idx(size_t r, size_t c) const { return r * m_ncol + c; }

Matrix& Matrix::itranspose()
{
    if (m_ncol != 1 && m_nrow != 1)
    {
        if (is_square())
            _itranspose_square();
        else
            _itranspose_cycle();
    }
    std::swap(m_nrow, m_ncol);
    return *this;
}

void Matrix::_itranspose_cycle()
{
    auto last = m_data + (m_nrow * m_ncol);
    const int mn1 = (last - m_data - 1);
    std::vector<bool> visited(last - m_data);
    auto cycle = m_data;
    while (++cycle != last)
    {
        if (visited[cycle - m_data])
            continue;
        int a = cycle - m_data;
        do
        {
            a = a == mn1 ? mn1 : (m_nrow * a) % mn1;
            std::swap(*(m_data + a), *cycle);
            visited[a] = true;
        } while ((m_data + a) != cycle);
    }
}

void Matrix::_itranspose_square()
{
    for (ccm::usize_t i = 1; i < m_nrow; i++)
    {
        for (ccm::usize_t j = 0; j < i; j++)
        {
            std::swap(this->operator()(i, j), this->operator()(j, i));
        }
    }
}

SubMatrix Matrix::slice(usize_t from_row, usize_t to_row, usize_t from_col, usize_t to_col)
{
    return SubMatrix(
      static_cast<Matrix&>(*this), from_row, from_col, (to_row - from_row), (to_col - from_col));
}
const SubMatrix
Matrix::slice(usize_t from_row, usize_t to_row, usize_t from_col, usize_t to_col) const
{
    return SubMatrix(*this, from_row, from_col, (to_row - from_row), (to_col - from_col));
}

SubMatrix Matrix::row(usize_t r)
{
    CCM_ASSERT((r < m_nrow), "Row out of range");
    return slice(r, r + 1, 0, m_ncol);
}
const SubMatrix Matrix::row(usize_t r) const
{
    CCM_ASSERT((r < m_nrow), "Row out of range");
    return slice(r, r + 1, 0, m_ncol);
}

SubMatrix Matrix::col(usize_t c)
{
    CCM_ASSERT((c < m_ncol), "Column out of range");
    return slice(0, m_nrow, c, c + 1);
}
const SubMatrix Matrix::col(usize_t c) const
{
    CCM_ASSERT((c < m_ncol), "Column out of range");
    return slice(0, m_nrow, c, c + 1);
}

// IO

std::istream& operator>>(std::istream& is, Matrix& M)
{
    auto old_prec = is.precision();
    is >> M.m_nrow >> M.m_ncol;
    is >> std::setprecision(12);
    delete[] M.m_data;
    M.m_data = new double[M.m_ncol * M.m_nrow];
    for (usize_t i = 0; i < M.m_ncol * M.m_nrow; i++)
        is >> M.m_data[i];
    is >> std::setprecision(old_prec);
    return is;
}

} // namespace ccm