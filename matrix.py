class Matrix:

    def __init__(self):
        pass

    def _validate_matrix_structure(self, inp_matrix):
        if not isinstance(inp_matrix,list) or not inp_matrix:
            raise ValueError("Input must be a non-empty list of lists.")
        
        if not all(isinstance(row,list)for row in inp_matrix):
            raise ValueError("All rows in the matrix must be lists.")

        if not all(len(row)==len(inp_matrix[0]) for row in inp_matrix):
            raise ValueError("All rows in the matrix must have the same number of columns.")
        return len(inp_matrix),len(inp_matrix[0])

    def _deep_copy(self, target_matrix):
        return [row[:] for row in target_matrix]

    def _is_close(self,a,b,rel_tol=1e-09,abs_tol=0.0):
        return abs(a-b)<=max(rel_tol*max((abs(a),abs(b))),abs_tol)

    def identity_matrix(self,size):
        return [[1 if i==j else 0 for i in range(size)]for j in range(size)]
    
    def zero_matrix(self, row,column):
        return [[0]* column for j in range(row)]
    
    def determinant(self,target_matrix):
        #return scalar
        self._validate_matrix_structure(target_matrix)
        if len(target_matrix)!=len(target_matrix[0]):
            raise ValueError("DETERMINANT:Matrix must be square.")

        if(len(target_matrix)==1):
            return target_matrix[0][0]

        l_matrix,u_matrix,p_matrix,number_of_substitutions=self.lu_decomposition(target_matrix,include_permutation_count=True)

        result_determinant=-1 if number_of_substitutions%2==1 else 1

        for row in range(len(u_matrix)):
            result_determinant*=u_matrix[row][row]

        return result_determinant


    def transpose(self,target_matrix):
        #return matrix
        self._validate_matrix_structure(target_matrix)
        result_matrix=[[target_matrix[i][j] for i in range(len(target_matrix))] for j in range(len(target_matrix[0]))]

        return result_matrix


    def multiply(self, matrix_a, matrix_b):
        #return matrix
        if isinstance(matrix_a,(int,float))or isinstance(matrix_b,(int,float)):
            if isinstance(matrix_a,(int,float)) and isinstance(matrix_b,(int,float)):
                return [[matrix_a*matrix_b]]
            else:
                if isinstance(matrix_a,(int,float)):
                    return [[i*matrix_a for i in j]for j in matrix_b]
                else:
                    return [[i*matrix_b for i in j]for j in matrix_a]
                
        self._validate_matrix_structure(matrix_a)

        self._validate_matrix_structure(matrix_b)
            
        if len(matrix_a[0])!=len(matrix_b):
            raise ValueError("MULTIPLY: Inner dimensions must match for matrix multiplication.")
        
        result_matrix=[[0 for j in range(len(matrix_b[0]))]for i in range(len(matrix_a))]

        for row_a in range(len(matrix_a)):
            for column_b in range(len(matrix_b[0])):
                for column_a in range(len(matrix_a[row_a])):
                    result_matrix[row_a][column_b]+=matrix_a[row_a][column_a]*matrix_b[column_a][column_b]

        return result_matrix

    def add(self,matrix_a,matrix_b):
        #return matrix
        self._validate_matrix_structure(matrix_a)

        self._validate_matrix_structure(matrix_b)
        if len(matrix_a)!= len(matrix_b) or len(matrix_a[0]) != len(matrix_b[0]):
            raise ValueError("ADD: Dimension mismatch for matrix addition.")

        result_matrix=[[0 for i in range(len(matrix_a[0]))]for j in range(len(matrix_a))]

        for row in range(len(matrix_a)):
            for column in range(len(matrix_a[row])):
                result_matrix[row][column]=matrix_a[row][column]+matrix_b[row][column]

        return result_matrix

    def inverse(self,inp_matrix):
        #return matrix
        self._validate_matrix_structure(inp_matrix)
        try:
            if len(inp_matrix) != len(inp_matrix[0]):
                raise ValueError("INVERSE: Matrix must be square.")
        
            if self._is_close(self.determinant(inp_matrix), 0.0):
                raise RuntimeError("INVERSE: Singular matrix (Determinant is zero).")

            id_matrix=[[1 if i==j else 0 for j in range(len(inp_matrix[0]))]for i in range(len(inp_matrix))]
            
            target_matrix=self._deep_copy(inp_matrix)
            augmented_matrix=self._deep_copy(id_matrix)

            #forward

            for row in range(len(target_matrix)):

                max_row=row
                max_val=target_matrix[row][row]

                for target_row in range(row,len(target_matrix)):
                    if abs(max_val)<abs(target_matrix[target_row][row]):
                        max_row=target_row
                        max_val=target_matrix[target_row][row]
                    
                target_matrix[row],target_matrix[max_row]=target_matrix[max_row],target_matrix[row]
                augmented_matrix[row],augmented_matrix[max_row]=augmented_matrix[max_row],augmented_matrix[row]

                factor=target_matrix[row][row]

                for column in range(len(target_matrix[row])):
                    target_matrix[row][column]=target_matrix[row][column]/factor
                    augmented_matrix[row][column]=augmented_matrix[row][column]/factor

                for target_row in range(row+1,len(target_matrix)):
                    factor=target_matrix[target_row][row]

                    for column in range(len(target_matrix[target_row])):
                        target_matrix[target_row][column]=target_matrix[target_row][column]-factor*target_matrix[row][column]
                        augmented_matrix[target_row][column]=augmented_matrix[target_row][column]-factor*augmented_matrix[row][column]

                if self._is_close(target_matrix[row][row], 0.0):
                    raise RuntimeError("INVERSE: ZERO PIVOT DETECTED")
                
            #backward
            for row in range(len(target_matrix)-1,-1,-1):
                for target_row in range(row):
                    factor=target_matrix[target_row][row] 

                    for column in range(len(target_matrix[target_row])):
                        target_matrix[target_row][column]=target_matrix[target_row][column]-factor*target_matrix[row][column]
                        augmented_matrix[target_row][column]=augmented_matrix[target_row][column]-factor*augmented_matrix[row][column]

            return augmented_matrix
        except ValueError as e:
            raise e
        except SystemExit:
            raise
        except RuntimeError as e:
            raise e
        except:
            raise RuntimeError("INVERSE: Unexpected error during inversion.")
        

    def lu_decomposition(self,inp_matrix,include_permutation_count=False,PLU=True):
        self._validate_matrix_structure(inp_matrix)

        #return matrix
        try:
            if not PLU and len(inp_matrix)!=len(inp_matrix[0]):
                raise ValueError("LU_DECOMPOSITION: Matrix must be square for standard LU decomposition (PLU=False).")
            target_matrix=self._deep_copy(inp_matrix)
            number_of_substitutions=0

            l_matrix=[[1 if i==j else 0 for j in range(len(target_matrix))]for i in range(len(target_matrix))]
            u_matrix=target_matrix
            p_matrix=[[1 if i==j else 0 for j in range(len(target_matrix))]for i in range(len(target_matrix))]

            for row in range(min(len(u_matrix),len(u_matrix[0]))):
                max_row=row
                max_val=target_matrix[row][row]

                if PLU==True:
                    for target_row in range(row,len(target_matrix)):
                        if abs(max_val)<abs(target_matrix[target_row][row]):
                            max_row=target_row
                            max_val=target_matrix[target_row][row]
                        
                    if row!= max_row: number_of_substitutions+=1
                    target_matrix[row],target_matrix[max_row]=target_matrix[max_row],target_matrix[row]
                    p_matrix[row],p_matrix[max_row]=p_matrix[max_row],p_matrix[row]
                    for column in range(row):
                        l_matrix[row][column],l_matrix[max_row][column]=l_matrix[max_row][column],l_matrix[row][column]

                    if self._is_close(target_matrix[row][row],0.0):
                        continue
                
                for target_row in range(row+1,len(u_matrix)):
                    factor=u_matrix[target_row][row]/u_matrix[row][row] if u_matrix[row][row]!=0 else 0
                    for column in range(len(u_matrix[target_row])):
                        u_matrix[target_row][column]=u_matrix[target_row][column]-factor*u_matrix[row][column]
                    
                    l_matrix[target_row][row]=factor
            
            if PLU==True:
                if include_permutation_count==False:
                    return l_matrix,u_matrix,p_matrix
                else:
                    return l_matrix,u_matrix,p_matrix,number_of_substitutions
            else:
                if include_permutation_count==False:
                    return l_matrix,u_matrix
                else:
                    return l_matrix,u_matrix,0

        except ValueError as e:
            raise e
        except SystemExit:
            raise
        except RuntimeError as e:
            raise e
        except:
            raise RuntimeError("LU_DECOMPOSITION: Unexpected error during decomposition.")


    def solve_system(self,coefficient_matrix,constant_matrix,least_square=True):
        self._validate_matrix_structure(coefficient_matrix)

        self._validate_matrix_structure(constant_matrix)

        #return matrix
        try:
            if least_square==True:
                l_matrix,u_matrix,p_matrix=self.lu_decomposition(self.multiply(self.transpose(coefficient_matrix),coefficient_matrix))
                swapped_constant_matrix=self.multiply(p_matrix,self.multiply(self.transpose(coefficient_matrix),constant_matrix))

            else:
                if len(coefficient_matrix) != len(constant_matrix) or len(coefficient_matrix) != len(coefficient_matrix[0]):
                    raise ValueError("SOLVE_SYSTEM: Coefficient matrix must be square and dimensions must match constant matrix")
                else:
                    l_matrix,u_matrix,p_matrix=self.lu_decomposition(coefficient_matrix)
                    swapped_constant_matrix=self.multiply(p_matrix,constant_matrix)
                    

            y_matrix=[0 for i in range(len(u_matrix))]
            x_matrix=[0 for i in range(len(l_matrix))]

            for row in range(0,len(l_matrix)):
                y_sum=0
                for column in range(row):
                    y_sum+=y_matrix[column]*l_matrix[row][column]

                
                y_matrix[row]=swapped_constant_matrix[row][0]-y_sum

            for row in range(len(u_matrix)-1,-1,-1):
                x_sum=0
                for column in range(row+1,len(u_matrix[row])):
                    x_sum+=u_matrix[row][column]*x_matrix[column]

                if self._is_close(u_matrix[row][row],0.0):
                    if not self._is_close(y_matrix[row]-x_sum,0.0):
                        return "NO SOLUTION"
                    else:
                        return "INFINITY SOLUTION"
                x_matrix[row]=(y_matrix[row]-x_sum)/u_matrix[row][row]

            return [[x] for x in x_matrix]
        
        except ValueError as e:
            raise e
        except SystemExit:
            raise 
        except RuntimeError as e:
            raise e
        except:
            raise RuntimeError("SOLVE_SYSTEM: Unexpected error during system solution.")
        