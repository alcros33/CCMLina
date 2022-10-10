#pragma once
#include "matrix.hpp"
namespace ccm
{
class BandMatrix : public MatrixTrait<BandMatrix, Matrix>
{
public:
    double* begin();
    const double* cbegin() const;
    double* end();
    const double* cend() const;

    size_t row_col_to_idx(size_t r, size_t c) const;

    friend std::istream& operator>>(std::istream& is, BandMatrix& M);
};
} // namespace ccm
