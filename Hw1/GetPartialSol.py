def GetPartialSol(matrixA, SizeN, Vectorb):
    """
    1. find max{a[i][j]} in row i
    2. swap rows.
    3. do the same as SimpleCalculateMatrix but only one itteration and then do again from one until we finish
    """

    # creat copies
    matrix = [list(map(float, row)) for row in matrixA]
    vector = list(map(float, Vectorb))

    temp = [0.0] * SizeN


    # forward part - upper triangle
    for i in range(SizeN):

        # PP
        # Find max line
        row = i
        maximum = 0
        for k in range(i, SizeN):
            if abs(matrix[k][i]) > maximum:
                maximum = abs(matrix[k][i])
                row = k

        # error if the entire line is zero
        if maximum == 0:
            return None

            # exchange lines  and vector if needed
        if row != i:
            matrix[i], matrix[row] = matrix[row], matrix[i]
            vector[i], vector[row] = vector[row], vector[i]

        # do single iterration of elimination - only on this line

        # Forward Propagation
        for k in range(i + 1, SizeN):
            pivot = matrix[i][i]
            l = matrix[k][i] / pivot

            # update matrix
            for j in range(i, SizeN):
                matrix[k][j] -= l * matrix[i][j]

            # update vector
            vector[k] -= l * vector[i]

    # Backward Propagation

    # Last line
    if abs(matrix[SizeN - 1][SizeN - 1]) < 1e-10: return None
    temp[SizeN - 1] = vector[SizeN - 1] / matrix[SizeN - 1][SizeN - 1]

    for i in range(SizeN - 2, -1, -1):
        sum1 = 0
        for j in range(i + 1, SizeN):
            sum1 += matrix[i][j] * temp[j]

        temp[i] = (vector[i] - sum1) / matrix[i][i]

    return temp