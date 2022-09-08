#include "linalg.hpp"
#include <gtest/gtest.h>

TEST(MatrixTest, AssignTest)
{

    ccm::Matrix M(3, 3, 0);
    M = 10;
    ccm::Matrix R(3, 3, 10);
    // Expect equality.
    EXPECT_EQ(M, R);
}

TEST(MatrixTest, AddititonTest)
{
    ccm::Matrix M(3, 3, 3);
    M += 10;
    ccm::Matrix R(3, 3, 13);
    // Expect equality.
    EXPECT_EQ(M, R);
}

TEST(MatrixTest, SubstractTest)
{
    ccm::Matrix M(3, 3, 3);
    M -= 10;
    ccm::Matrix R(3, 3, -7);
    // Expect equality.
    EXPECT_EQ(M, R);
}

TEST(MatrixTest, MatmulTest)
{
    ccm::Matrix A{{-5, 5, -4}, {-2, 2, 2}};
    ccm::Matrix B{{-3, -6}, {2, 3}, {6, 1}};
    ccm::Matrix Res{{1, 41}, {22, 20}};
    // Expect equality.
    EXPECT_EQ(A.matmul(B), Res);
}