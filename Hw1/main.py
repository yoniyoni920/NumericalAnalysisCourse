from GetPartialSol import GetPartialSol
from SimpleCalculateMatrix import SimpleCalculateMatrix
def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.
def main(matrix_a, size_n, vector_b):
    if is_single_sol(matrix_a):
        simple_answer = SimpleCalculateMatrix(matrix_a, size_n, vector_b)
        partial_sol = GetPartialSol(matrix_a, size_n, vector_b)
        print(get_closeness_vector(simple_answer, partial_sol))
def is_single_sol(matrix_a):
    """
    checks if the matrix has single solution
    """
    value = False
    determinate = get_determinant(matrix_a)
    print(f"Calculated Determinant: {determinate}")
    # if det is not zero with small margin return true
    if abs(determinate) > 1e-10:
        return True
    return False
# will return the minor basicly a smaller matrix without  0 line  and j column
def get_minor(matrix, j):
    return [row[:j] + row[j+1:] for row in (matrix[:0] + matrix[0+1:])]


def get_determinant(matrix_a):
    """
    calculate the determinant in a recursive way
    """
    size_n = len(matrix_a)
    # if matrix is sizeN == 1 then matrix size is 1 on 1
    if size_n == 1:
        return matrix_a[0][0]
    det = 0
    # go over all the first row
    for i in range(size_n):
         # change the mark so that even times will be + and uneven -
        sign = (-1) ** i
        current = matrix_a[0][i]
         #get the minor
        inner_matrix = get_minor(matrix_a, i)
        inner_det = get_determinant(inner_matrix)
        det += sign * current * inner_det

    return det

def get_closeness_vector(simple_answer, partial_sol):
    if not simple_answer  or not partial_sol:
        print("error empty vector")
        exit(1)
    error_vector = []
    #calculate the error vector before printing:
    if len(simple_answer) != len(partial_sol):
        print("vector size mismatch")
        exit(1)
    i=0
    for simpleElement,partialElement in zip(simple_answer, partial_sol):
        if simpleElement == 0:
            error_vector.append(0)  # protect from zero devision
        else:
            error_vector.append(abs((simpleElement - partialElement)/simpleElement))
    return str(error_vector)


if __name__ == '__main__':
    print("=== simple solution  ===")

    # solve able solution  [1, 1, 1, 1]
    matrix_solvable = [
        [4, 1, 1, 1],
        [1, 4, 1, 1],
        [1, 1, 4, 1],
        [1, 1, 1, 4]
    ]
    vector_solvable = [7, 7, 7, 7]
    size_n = 4

    print(f"Matrix A: {matrix_solvable}")
    print(f"Vector b: {vector_solvable}")
    main(matrix_solvable, size_n, vector_solvable)

    print("\n" + "=" * 50 + "\n")

    print("=== single solution more complex matrix where we can see the PP algorithem effects ===")

    matrix_singular1 = [
        [0, 2, 0, 1],
        [2, 2, 3, 2],
        [4, -3, 0, 1],
        [6, 1, -6, -5]
    ]
    vector_solvable1 = [0, -2, -7, 6]


    print(f"Matrix A: {matrix_singular1}")
    print(f"Vector b: {vector_solvable1}")
    main(matrix_singular1, size_n, vector_solvable1)


    print("=== Not single solution  ===")
    # det = 0 second line is double the first one
    matrix_singular = [
        [1, 2, 3, 4],
        [2, 4, 6, 8],
        [5, 6, 7, 8],
        [9, 10, 11, 12]
    ]
    vector_singular = [10, 20, 30, 40]

    print(f"Matrix A: {matrix_singular}")
    print(f"Vector b: {vector_singular}")
    main(matrix_singular, size_n, vector_singular)
