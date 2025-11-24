def SimpleCalculateMatrix(matrixA, SizeN, Vectorb):
    """
     just simple gaos solution
     get upper traingle and then calculate with the solution we get
    """

    # copy matrix and vector
    matrix = [list(map(float, row)) for row in matrixA]
    vector = list(map(float, Vectorb))

    temp = [0.0] * SizeN  # solution will be here

    # forward part - upper triangle
    for i in range(SizeN):

        # if the pivot is 0 we need to change the line so we can continue
        if abs(matrix[i][i]) == 0:
            found_swap = False
            for k in range(i + 1, SizeN):
                #find pyvot that isnt zero and exchange between the lines
                if abs(matrix[k][i]) != 0:
                    matrix[i], matrix[k] = matrix[k], matrix[i]
                    vector[i], vector[k] = vector[k], vector[i]
                    found_swap = True
                    break

        # Forward Propagation
        for k in range(i + 1, SizeN):
            pivot =  matrix[i][i]
            l = matrix[k][i] / pivot

            # update matrix
            for j in range(i, SizeN):
                matrix[k][j] -= l * matrix[i][j]

            # update vector
            vector[k] -= l * vector[i]

    # backward propagation
    #do the same as goes and just do all the equations without updateing the matrix
    # last line
    temp[SizeN - 1] = vector[SizeN - 1] / matrix[SizeN - 1][SizeN - 1]
    # 2.the rest
    for i in range(SizeN - 2, -1, -1):
        sum1 = 0
        for j in range(i + 1, SizeN):
            sum1 += matrix[i][j] * temp[j]

        temp[i] = (vector[i] - sum1) / matrix[i][i]
    return temp
