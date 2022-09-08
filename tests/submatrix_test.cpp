#include "linalg.hpp"
#include <gtest/gtest.h>

TEST(SubMatrixTest, AssignTest)
{
    ccm::Matrix M(5, 5, 5);
    auto sm = M.slice(1, 3, 1, 3);
    ccm::Matrix A(2, 2, 3);
    sm = A;
    ccm::Matrix R{
      {5, 5, 5, 5, 5}, {5, 3, 3, 5, 5}, {5, 3, 3, 5, 5}, {5, 5, 5, 5, 5}, {5, 5, 5, 5, 5}};
    ASSERT_EQ(sm, A);
    ASSERT_EQ(sm, R.slice(1, 3, 1, 3));
    ASSERT_EQ(M, R);
}

TEST(SubMatrixTest, AssignTestNum)
{
    ccm::Matrix M(5, 5, 5);
    auto sm = M.slice(1, 3, 1, 3);
    sm = 3;
    ccm::Matrix R{
      {5, 5, 5, 5, 5}, {5, 3, 3, 5, 5}, {5, 3, 3, 5, 5}, {5, 5, 5, 5, 5}, {5, 5, 5, 5, 5}};
    ASSERT_EQ(sm, R.slice(1, 3, 1, 3));
    ASSERT_EQ(M, R);
}

TEST(SubMatrixTest, AditionTest)
{
    ccm::Matrix M(5, 5, 5);
    auto sm = M.slice(1, 3, 1, 3);
    ccm::Matrix A(2, 2, 4);
    sm += A;
    ccm::Matrix R{
      {5, 5, 5, 5, 5}, {5, 9, 9, 5, 5}, {5, 9, 9, 5, 5}, {5, 5, 5, 5, 5}, {5, 5, 5, 5, 5}};
    ASSERT_EQ(sm, R.slice(1, 3, 1, 3));
    ASSERT_EQ(M, R);
}

TEST(SubMatrixTest, AditionTestNum)
{
    ccm::Matrix M(5, 5, 5);
    auto sm = M.slice(1, 3, 1, 3);
    sm += 3;
    ccm::Matrix R{
      {5, 5, 5, 5, 5}, {5, 8, 8, 5, 5}, {5, 8, 8, 5, 5}, {5, 5, 5, 5, 5}, {5, 5, 5, 5, 5}};
    ASSERT_EQ(M, R);
}

TEST(SubMatrixTest, CopyAssignTest)
{
    ccm::Matrix M(5, 5, 5);
    auto sm = M.slice(1, 3, 1, 3);
    ccm::Matrix A(sm);
    ASSERT_EQ(A, sm);
}

TEST(SubMatrixTest, MatmulTestA)
{
    ccm::Matrix A{{0, -5, 5, -4, 0}, {0, -2, 2, 2, 0}};
    ccm::Matrix B{{-3, -6}, {2, 3}, {6, 1}};
    ccm::Matrix Res{{1, 41}, {22, 20}};
    // Expect equality.
    auto R = A.slice(0, 2, 1, 4).matmul(B);
    EXPECT_EQ(R, Res);
}
TEST(SubMatrixTest, MatmulTestB)
{
    ccm::Matrix A{{-5, 5, -4}, {-2, 2, 2}};
    ccm::Matrix B{{0, 0}, {-3, -6}, {2, 3}, {6, 1}};
    ccm::Matrix Res{{1, 41}, {22, 20}};
    // Expect equality.
    auto R = A.matmul(B.slice(1, 4, 0, 2));
    EXPECT_EQ(R, Res);
}