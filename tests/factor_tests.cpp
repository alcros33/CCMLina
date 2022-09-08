#include "linalg.hpp"
#include <gtest/gtest.h>

inline bool mat_near(const ccm::Matrix& A, const ccm::Matrix& B, double eps = 1E-5)
{
    return A.dist(B) < eps;
}

#define TEST_TOL 1E-5

#define EXPECT_MAT_NEAR(A, B) EXPECT_PRED3(mat_near, A, B, TEST_TOL)

TEST(FactorTest, CholeskyTest)
{
    ccm::Matrix A{{4, -1, 1}, {-1, 4.25, 2.75}, {1, 2.75, 3.5}};
    ccm::Matrix R{{2, 0, 0}, {-0.5, 2, 0}, {0.5, 1.5, 1}};
    auto L = ccm::factor_cholesky(A);
    EXPECT_MAT_NEAR(L, R);
}