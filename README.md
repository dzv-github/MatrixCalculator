# MatrixCalculator

identity_matrix(size): creates $size \\times size$ identity matrix and return result as lists of lists  
  
zero_matrix(row,column): creates $row \\times column$ zero matrix and return result as list of lists  
  
determinant(target_matrix): calculates $det(target \\_ matrix)$ and return as an int or float  
  
transpose(target_matrix): creates transposed matrix and return as list of lists  
  
multiply(matrix_a,matrix_b): calculates $matrix \\_ a \\times matrix \\_ b$ and return as list of lists  
  
add(matrix_a,matrix_b): calculates $matrix \\_ a + matrix \\_ b$ and return as list of lists  
  
inverse(inp_matrix): creates inverse matrix and return as list of lists  
  
lu_decomposition(inp_matrix, include_permutation_count=False, PLU=True): do LU decomposition and return l_matrix and u_matrix as list of lists. if PLU is true, do PLU decomposition and return p_matrix as list of lists in addition. if include_permutation_count is true, return number of row substitutions at pivoting as int in addition  
  
solve_system(coefficient_matrix, constant_matrix, least_square=True): solve the $coefficient \\_ matrix \\times [x] = constant \\_ matrix$ equation and return solution of it as list of lists. but if solution is infinity or no solution, return as string. if least_square is true, do least squares method and return a result of it as list of lists instead of return as string.
