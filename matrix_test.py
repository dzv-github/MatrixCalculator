# matrix_test.py

import unittest
from matrix import Matrix 
import math 

class TestMatrixOperations(unittest.TestCase):
    
    def setUp(self):
        self.m = Matrix()

    # ----------------------------------------------------------------------
    # A. 구조 및 유효성 검사 (Validation & Structure)
    # ----------------------------------------------------------------------

    def test_a1_exception_ragged_and_type_mismatch(self):
        # 1. 들쭉날쭉한 행렬 (Ragged matrix)
        with self.assertRaisesRegex(ValueError, "same number of columns"):
            self.m.add([[1, 2], [3]], [[1, 1], [1, 1]])
            
        # 2. 행렬 내부에 리스트가 아닌 타입 포함
        with self.assertRaisesRegex(ValueError, "must be lists"):
             self.m.transpose([[1, 2], 'a'])
             
        # 3. 빈 입력
        with self.assertRaisesRegex(ValueError, "non-empty list of lists"):
             self.m.transpose([])

    def test_a2_dimension_mismatch(self):
        # 1. 덧셈 차원 불일치
        with self.assertRaisesRegex(ValueError, "Dimension mismatch for matrix addition"):
            self.m.add([[1, 2]], [[3, 4], [5, 6]])
            
        # 2. 곱셈 차원 불일치
        with self.assertRaisesRegex(ValueError, "Inner dimensions must match"):
            self.m.multiply([[1, 2]], [[1], [2], [3]])
            
        # 3. 비정방 행렬 역행렬/행렬식
        A_non_square = [[1, 2, 3], [4, 5, 6]]
        with self.assertRaisesRegex(ValueError, "Matrix must be square"):
            self.m.inverse(A_non_square)
        with self.assertRaisesRegex(ValueError, "Matrix must be square"):
            self.m.determinant(A_non_square)

    # ----------------------------------------------------------------------
    # B. 기본 산술 및 스칼라 연산 (Basic Arithmetic & Scalar)
    # ----------------------------------------------------------------------

    def test_b1_basic_add_multiply_transpose(self):
        A = [[1, 2], [3, 4]]
        B = [[5, 6], [7, 8]]
        
        # 덧셈
        self.assertEqual(self.m.add(A, B), [[6, 8], [10, 12]])
        # 곱셈
        self.assertEqual(self.m.multiply(A, B), [[19, 22], [43, 50]])
        # 전치
        self.assertEqual(self.m.transpose(A), [[1, 3], [2, 4]])
        # 1xn 행렬 전치
        A_1xn = [[1, 2, 3]]
        self.assertEqual(self.m.transpose(A_1xn), [[1], [2], [3]])

    def test_b2_scalar_multiplication_and_zero_matrix(self):
        A = [[1, 2], [3, 4]]
        scalar = 2.5
        
        # 스칼라 * 행렬
        self.assertEqual(self.m.multiply(scalar, A), [[2.5, 5.0], [7.5, 10.0]])
        # 0 * 행렬 (경계 케이스)
        self.assertEqual(self.m.multiply(0, A), [[0, 0], [0, 0]])


    # ----------------------------------------------------------------------
    # C. 선형 시스템 해법 (Solve System)
    # ----------------------------------------------------------------------
    
    def test_c1_solve_system_unique_solution(self):
        # 고유 해: x=1.6, y=1.8
        A = [[4, 2], [1, 3]]
        b = [[10], [7]]
        x = self.m.solve_system(A, b, least_square=False) 
        
        self.assertAlmostEqual(x[0][0],1.6) 
        self.assertAlmostEqual(x[1][0], 1.8) 

    def test_c2_solve_system_no_solution_and_infinite(self):
        # 1. 해 없음
        A_no = [[1, 1], [1, 1]] 
        b_no = [[1], [2]]       
        self.assertEqual(self.m.solve_system(A_no, b_no, least_square=False), "NO SOLUTION")

        # 2. 무한대 해
        A_inf = [[1, 1], [2, 2]] 
        b_inf = [[1], [2]]
        self.assertEqual(self.m.solve_system(A_inf, b_inf, least_square=False), "INFINITY SOLUTION")

    def test_c3_solve_system_least_squares(self):
        # 🚩 수정됨: 수학적 정답은 x=1, y=1 입니다. (4/3 -> 1.0)
        A = [[1, 0], [1, 1], [0, 1]] 
        b = [[1], [2], [1]] 
        x = self.m.solve_system(A, b, least_square=True)
        
        self.assertAlmostEqual(x[0][0], 1.0) 
        self.assertAlmostEqual(x[1][0], 1.0) 


    # ----------------------------------------------------------------------
    # D. 특이점 및 정밀도 테스트 (Singularity & Precision)
    # ----------------------------------------------------------------------

    def test_d1_determinant_enhanced(self):
        # 🚩 수정됨: 3x3 행렬식의 수학적 정답은 5.0 입니다. (7.0 -> 5.0)
        A_3x3 = [[3, 2, 0], [1, 4, 5], [0, 1, 2]]
        self.assertAlmostEqual(self.m.determinant(A_3x3), 5.0)

        # 2. 치환이 필요한 행렬 (Determinant -2)
        A_swap = [[0, 1], [2, 3]]
        self.assertAlmostEqual(self.m.determinant(A_swap), -2.0)
        
        # 3. 상/하삼각 행렬 (대각선 곱만 확인)
        A_triangular = [[2, 0, 0], [4, 5, 0], [6, 7, 8]]
        self.assertAlmostEqual(self.m.determinant(A_triangular), 80.0)
        
    def test_d2_inverse_singular_matrix(self):
        # 1. 특이 행렬 (Determinant 0)
        A_singular = [[1, 1], [1, 1]]
        with self.assertRaisesRegex(RuntimeError, "Singular matrix"):
            self.m.inverse(A_singular)
            
        # 2. Zero Pivot 검출
        A_zero_pivot = [[0, 1], [0, 2]]
        with self.assertRaisesRegex(RuntimeError, "Singular matrix|ZERO PIVOT"):
            self.m.inverse(A_zero_pivot)

    def test_d3_is_close_precision(self):
        # 1. 톨러런스 범위 내 비교 (Success)
        a = 1.0
        b = 1.0 + 5e-10 
        self.assertTrue(self.m._is_close(a, b))

        # 2. 톨러런스 범위 초과 비교 (Fail)
        c = 1.0 + 2e-09 
        self.assertFalse(self.m._is_close(a, c))

if __name__ == '__main__':
    unittest.main(argv=['first-arg-is-ignored'], exit=False)