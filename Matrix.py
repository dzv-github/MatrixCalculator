class Matrix:

    def __init__(self):
        pass

    def __MatrixdeepCopy(self, targetMatrix):
        return [row[:] for row in targetMatrix]

    def __MatrixIsClose(self,a,b,rel_tol=1e-09,abs_tol=0.0):
        return abs(a-b)<=max(rel_tol*max((abs(a),abs(b))),abs_tol)

    def __IntSwap(self,a,b):
        tmp=a
        a=b
        b=tmp
        return a,b


    def rawSwap(self,RawA,RawB):
        tmpMatrix=RawA
        RawA=RawB
        RawB=tmpMatrix
        return RawA,RawB

    def Det(self,targetMatrix):
        if len(targetMatrix)!=len(targetMatrix[0]):
            raise Exception("DET:THIS IS NOT SQUARE MATRIX")

        if(len(targetMatrix)==1):
            return targetMatrix[0][0]

        Lmatrix,Umatrix,Pmatrix,NumberOfSubstitutions=self.LU(targetMatrix,returnNumberOfSubstitution=True)

        Det=-1 if NumberOfSubstitutions%2==1 else 1

        for raw in range(len(Umatrix)):
            Det*=Umatrix[raw][raw]

        return Det


    def Transposition(self,targetMatrix):
        resultMatrix=[[targetMatrix[i][j] for i in range(len(targetMatrix))] for j in range(len(targetMatrix[0]))]

        return resultMatrix


    def Multiply(self, matrixA, matrixB):
        if len(matrixA[0])!=len(matrixB):
            raise Exception("MULTIPLY: MULTIPLY DIMENSION MISMATCH")
        
        resultMatrix=[[0 for j in range(len(matrixB[0]))]for i in range(len(matrixA))]

        for rawA in range(len(matrixA)):
            for columnB in range(len(matrixB[0])):
                for columnA in range(len(matrixA[rawA])):
                    resultMatrix[rawA][columnB]+=matrixA[rawA][columnA]*matrixB[columnA][columnB]

        return resultMatrix

    def Sum(self,matrixA,matrixB):
        if len(matrixA)!= len(matrixB) or len(matrixA[0]) != len(matrixB[0]):
            raise Exception("SUM: DIMENSION MISMATCH")

        resultMatrix=[[0 for i in range(len(matrixA[0]))]for j in range(len(matrixA))]

        for raw in range(len(matrixA)):
            for column in range(len(matrixA[raw])):
                resultMatrix[raw][column]=matrixA[raw][column]+matrixB[raw][column]

        return resultMatrix

    def Inverse(self,INPmatrix):
        if len(INPmatrix) != len(INPmatrix[0]):
            raise Exception("INVERSE: MATRIX IS NOT SQUARE MATRIX")


        try:
            if self.__MatrixIsClose(self.Det(INPmatrix), 0.0):
                raise Exception("INVERSE: SINGULAR MATRIX (DET IS ZERO)")
        except SystemExit:
            raise
        except:
            raise Exception("INVERSE: FAILED DETERMINANT CHECK")

        IDmatrix=[[1 if i==j else 0 for j in range(len(INPmatrix[0]))]for i in range(len(INPmatrix))]

        targetMatrix=self.__MatrixdeepCopy(INPmatrix)
        augmentedMatrix=self.__MatrixdeepCopy(IDmatrix)

        #forward

        for raw in range(len(targetMatrix)):

            maxraw=raw
            maxval=targetMatrix[raw][raw]

            for targetRaw in range(raw,len(targetMatrix)):
                if abs(maxval)<abs(targetMatrix[targetRaw][raw]):
                    maxraw=targetRaw
                    maxval=targetMatrix[targetRaw][raw]
                
            targetMatrix[raw],targetMatrix[maxraw]=self.rawSwap(targetMatrix[raw],targetMatrix[maxraw])
            augmentedMatrix[raw],augmentedMatrix[maxraw]=self.rawSwap(augmentedMatrix[raw],augmentedMatrix[maxraw])

            factor=targetMatrix[raw][raw]

            for column in range(len(targetMatrix[raw])):
                targetMatrix[raw][column]=targetMatrix[raw][column]/factor
                augmentedMatrix[raw][column]=augmentedMatrix[raw][column]/factor

            for targetRaw in range(raw+1,len(targetMatrix)):
                factor=targetMatrix[targetRaw][raw]

                for column in range(len(targetMatrix[targetRaw])):
                    targetMatrix[targetRaw][column]=targetMatrix[targetRaw][column]-factor*targetMatrix[raw][column]
                    augmentedMatrix[targetRaw][column]=augmentedMatrix[targetRaw][column]-factor*augmentedMatrix[raw][column]

            if self.__MatrixIsClose(targetMatrix[raw][raw], 0.0):
                raise Exception("INVERSE: ZERO PIVOT DETECTED")
            
        #backward
        for raw in range(len(targetMatrix)-1,-1,-1):
            for targetRaw in range(raw):
                factor=targetMatrix[targetRaw][raw] 

                for column in range(len(targetMatrix[targetRaw])):
                    targetMatrix[targetRaw][column]=targetMatrix[targetRaw][column]-factor*targetMatrix[raw][column]
                    augmentedMatrix[targetRaw][column]=augmentedMatrix[targetRaw][column]-factor*augmentedMatrix[raw][column]

        return augmentedMatrix

    def LU(self,INPmatrix,returnNumberOfSubstitution=False,PLU=True):
        try:
            if len(INPmatrix)==0:
                raise Exception("LU: EMPTY MATRIX")
            targetMatrix=self.__MatrixdeepCopy(INPmatrix)
            NumberOfSubstitution=0

            Lmatrix=[[1 if i==j else 0 for j in range(len(targetMatrix))]for i in range(len(targetMatrix))]
            Umatrix=targetMatrix
            Pmatrix=[[1 if i==j else 0 for j in range(len(targetMatrix))]for i in range(len(targetMatrix))]

            for raw in range(len(Umatrix)):
                maxraw=raw
                maxval=targetMatrix[raw][raw]

                if PLU==True:
                    for targetRaw in range(raw,len(targetMatrix)):
                        if abs(maxval)<abs(targetMatrix[targetRaw][raw]):
                            maxraw=targetRaw
                            maxval=targetMatrix[targetRaw][raw]
                        
                    if raw!= maxraw: NumberOfSubstitution+=1
                    targetMatrix[raw],targetMatrix[maxraw]=self.rawSwap(targetMatrix[raw],targetMatrix[maxraw])
                    Pmatrix[raw],Pmatrix[maxraw]=self.rawSwap(Pmatrix[raw],Pmatrix[maxraw])
                    for column in range(raw):
                        Lmatrix[raw][column],Lmatrix[maxraw][column]=self.__IntSwap(Lmatrix[raw][column],Lmatrix[maxraw][column])

                    if self.__MatrixIsClose(targetMatrix[raw][raw],0.0):
                        continue
                
                for targetraw in range(raw+1,len(Umatrix)):
                    factor=Umatrix[targetraw][raw]/Umatrix[raw][raw] if Umatrix[raw][raw]!=0 else 0
                    for column in range(len(Umatrix[targetraw])):
                        Umatrix[targetraw][column]=Umatrix[targetraw][column]-factor*Umatrix[raw][column]
                    
                    Lmatrix[targetraw][raw]=factor
            
            if PLU==True:
                if returnNumberOfSubstitution==False:
                    return Lmatrix,Umatrix,Pmatrix
                else:
                    return Lmatrix,Umatrix,Pmatrix,NumberOfSubstitution
            else:
                if returnNumberOfSubstitution==False:
                    return Lmatrix,Umatrix
                else:
                    return Lmatrix,Umatrix,0

        except:
            raise Exception("LU : ERROR IN LU DECOMPOSITIONED")


    def EquationSol(self,CoefficientMatrix,ConstantMatrix,LeastSquare=True):
        try:
            if LeastSquare==True:

                Lmatrix,Umatrix,Pmatrix=self.LU(self.Multiply(self.Transposition(CoefficientMatrix),CoefficientMatrix))
                SwappedConstantMatrix=self.Multiply(Pmatrix,self.Multiply(self.Transposition(CoefficientMatrix),ConstantMatrix))
            else:
                if len(CoefficientMatrix) != len(ConstantMatrix) or len(CoefficientMatrix) != len(CoefficientMatrix[0]):
                    raise Exception("ERROR")
                else:
                    Lmatrix,Umatrix,Pmatrix=self.LU(CoefficientMatrix)
                    SwappedConstantMatrix=self.Multiply(Pmatrix,ConstantMatrix)
                    

            ymatrix=[0 for i in range(len(Umatrix))]
            xmatrix=[0 for i in range(len(Lmatrix))]

            for raw in range(0,len(Lmatrix)):
                ysum=0
                for column in range(raw):
                    ysum+=ymatrix[column]*Lmatrix[raw][column]

                
                ymatrix[raw]=SwappedConstantMatrix[raw][0]-ysum

            for raw in range(len(Umatrix)-1,-1,-1):
                xsum=0
                for column in range(raw+1,len(Umatrix[raw])):
                    xsum+=Umatrix[raw][column]*xmatrix[column]

                if self.__MatrixIsClose(Umatrix[raw][raw],0.0):
                    if not self.__MatrixIsClose(ymatrix[raw]-xsum,0.0):
                        return "NO SOLUTION"
                    else:
                        return "INFINITY SOLUTION"
                xmatrix[raw]=(ymatrix[raw]-xsum)/Umatrix[raw][raw]

            return [[x] for x in xmatrix]
        
        except SystemExit:
            raise # Multiply 또는 LU 실패 오류 전파
        except:
            raise Exception("EQUATIONSOL: ERROR")
