import numpy as np
import re
import sklearn.preprocessing
from scipy.sparse import csr_matrix
from scipy.sparse import dok_matrix
from scipy.sparse import random

np.set_printoptions(suppress=True)

INFLATION_OPERATOR = 2
EXPANSION_OPERATOR = 2
THRESHOLD_A = 0.01
THRESHOLD_B = 4.0

def expand(matrix):
    return matrix.dot(matrix)

def normalize(matrix):
    return sklearn.preprocessing.normalize(matrix, norm="l1", axis=0)

def inflate(matrix, inflation_op):
    return normalize(matrix.power(inflation_op))

def selfLoop(matrix, dimension):
    for i in range(dimension):
        matrix[i][i] = 1

    return matrix


"""
Based of off numpy.allclose

atol: absolute tolerance parameter
rtol: relative tolerance parameter
"""
def convergence(matrix1, matrix2, rtol=1e-05, atol=1e-08):
    allclose = (np.abs(matrix1 - matrix2)) - (rtol * np.abs(matrix2))
    return allclose.max() <= atol


def get_diag(matrix, dimensions):
    nonzero_diag = list()
    for i in range(dimensions):
        x = float("{:.5f}".format(float(matrix.item(i, i))))
        if x != 0:
            nonzero_diag.append(i)

    return nonzero_diag

def get_clusters(matrix, dimension):
    #attractors = get_diag(matrix.todense(), dimension)
    attractors = matrix.diagonal().nonzero()[0].tolist()

    clusters = list()
    for row in attractors:
        elements = matrix.getrow(row).nonzero()[1].tolist()
        if (len(elements) > 1):
                clusters.append(elements)

    return clusters

def threshold(matrix, dimensions):
    max_values = matrix.max(axis=0).toarray().tolist()[0]
    sum_ctr = list()

    for col in range(dimensions):
        sum = 0
        rows = matrix.getcol(col).nonzero()[0].tolist()
        for row in rows:
            sum = sum + (matrix[row, col]**2)

        sum_ctr.append(sum) 

    col_threshold = list()
    
    for col in range(dimensions):
        ctr = sum_ctr[col]
        col_threshold.append( (THRESHOLD_A*ctr) * (1 - (THRESHOLD_B * (max_values[col] - ctr))))

    return col_threshold


def prune(matrix, dimensions):
    threshold_values = threshold(matrix, dimensions)

    for col in range(dimensions):
        rows = matrix.getcol(col).nonzero()[0].tolist()

        for row in rows:
            if matrix[row, col] <= threshold_values[col]:
                matrix[row, col] = 0
 

    return normalize(matrix)


#THIS IS TO BE CALLED
def mcl_algorithm(matrix, dimensions):
    selfLoop(matrix, dimensions)

    sp_matrix_orig = csr_matrix(np.matrix(matrix))
    
    normalize(sp_matrix_orig)

    sp_matrix_eval = inflate( expand(sp_matrix_orig.copy()), INFLATION_OPERATOR)

    while not convergence(sp_matrix_eval, sp_matrix_orig):
        sp_matrix_orig = sp_matrix_eval
        sp_matrix_eval = prune( inflate(expand(sp_matrix_orig), INFLATION_OPERATOR), dimensions)

    clusters = get_clusters(sp_matrix_eval, dimensions)
    print(sp_matrix_eval.todense())
    
    #clusters = get_clusters(sp_matrix_eval, dimensions)
    print(len(clusters))

    return clusters


def split_lines(data):
    new_data = list()
    length = len(data)
    for el in data[:length-1]:
        add = el.split('\t')
        add[0] = int(add[0])
        add[1] = int(add[1])
        new_data.append(add)

   # print(new_data)
    for e in new_data:
        if len(e) != 2:
            print (e)
        if type(e[0]) is not int or type(e[1]) is not int:
            print("Not Int")
            print(e)
    return new_data

def form_matrix():
    with open("dataset/numbered_data", 'r') as file:
        data = file.read().split('\n')
        dimensions = int (re.findall(r'\d+', data[0])[0])

    new_data = split_lines(data[1:])
    
    matrix = [0] * dimensions
    for i in range(dimensions):
        matrix[i] = [0] * dimensions

    for el in new_data:
        if matrix[el[1]][el[0]] != 1:
            matrix[el[0]][el[1]] = 1

    return (matrix, dimensions)
