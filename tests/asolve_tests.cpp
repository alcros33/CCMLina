#include "linalg.hpp"
#include <gtest/gtest.h>

#define TEST_TOL 1E-5

inline bool mat_near(const ccm::Matrix& A, const ccm::Matrix& B) { return A.dist(B) < TEST_TOL; }

#define EXPECT_MAT_NEAR(A, B) EXPECT_PRED2(mat_near, A, B)

TEST(ASolveTest, LowerTest)
{
    ccm::Matrix L{{1, 0, 0, 0}, {2, 1, 0, 0}, {3, 4, 1, 0}, {-1, -3, 0, 1}};
    ccm::Matrix b{8, 7, 14, -7};
    ccm::Matrix R{8, -9, 26, -26};
    auto x = ccm::solve_lower(L, b);
    EXPECT_MAT_NEAR(x, R);
}

TEST(ASolveTest, UpperTest)
{
    ccm::Matrix U{{1, 1, 0, 3}, {0, -1, -1, -5}, {0, 0, 3, 13}, {0, 0, 0, -13}};
    ccm::Matrix b{8, -9, 26, -26};
    ccm::Matrix R{3, -1, 0, 2};
    auto x = ccm::solve_upper(U, b);
    EXPECT_MAT_NEAR(x, R);
}

TEST(ASolveTest, GaussElimTest)
{
    ccm::Matrix A{{1, -1, 2, -1}, {2, -2, 3, -3}, {1, 1, 1, 0}, {1, -1, 4, 3}};
    ccm::Matrix b{-8, -20, -2, 4};
    ccm::Matrix R{-7, 3, 2, 2};
    auto x = ccm::solve_gauss(A, b);
    EXPECT_MAT_NEAR(x, R);
}