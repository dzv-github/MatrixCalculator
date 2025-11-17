# MatrixCalculator

identity_matrix(size): create $size \\times size$ identity matrix and return result as lists in list  
  
zero_matrix(row,column): create $row \\times column$ zero matrix and return result as lists in list  
  
determinant(target_matrix): calculate $det(target \\_ matrix)$ and return as int or float  
  
transpose(target_matrix): create transposed matrix and return as lists in list  
  
multiply(matrix_a,matrix_b): calculate $matrix \\_ a \\times matrix \\_ b$ and return as lists in list  
  
add(matrix_a,matrix_b): calculate $matrix \\_a + matrix \\_ b$ and return as lists in list  
  
inverse(inp_matrix): create inverse matrix and return as lists in list  
  
lu_decomposition(inp_matrix, include_permutation_count=False, PLU=True): do LU decomposition and return l_matrix and u_matrix as lists in list. if PLU is true, do PLU decomposition and return p_matrix as lists in                                                                            list in addition. if include_permutation_count is true, return number of row substitutions at pivoting as int in addition  
  
solve_system(coefficient_matrix, constant_matrix, least_square=True): solve the $coefficient \\_ matrix \\times [x] = constant \\_ matrix$ equation and return solution of it as lists in list. but if solution is infinity or non,                                                                        return as string. if least_square is true, do least squares method and return a result of it as lists of list instead of return as string                                                                              sometimes.
