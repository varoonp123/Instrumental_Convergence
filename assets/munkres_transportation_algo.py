# @Author: Varoon Pazhyanur <varoon>
# @Date:   05-11-2017
# @Filename: munkres_assignment_algo.py
# @Last modified by:   varoon
# @Last modified time: Jan-11-2018
import random as r
import copy
import itertools    
import numpy as np      

'''
Algorithm in O(n^3) to solve the assignment problem as described in
"Algorithms for the Assignment and Transportation Problems" by James Munkres in the Journal of the Society of Industrial and Applied Math March 1957.

Problem: Given n people, n jobs and a rating matrix (r_{i,j}) describing a rating of person i in job j, assign each person to exactly 1 job to maximize the total rating. ie Find n independent (do not share a row or column) elements in the matrix with maximal sum. This is equivalent to minimizing that sum with a nonnegative matrix. This script runs through the algorithm with a random nxn matrix and gets the min sum. It also compares the output to a large number of options generated in the naive method.
'''

#get a set of ind starred zeroes and cover their cols. This is always possible because when this runs, each row and each col must have atleast one 0. 
def star_ind_zeros_and_cover(m, star,prime, cov_row, cov_col):
    #Note: This can probably be done more efficiently (but probably not with numpy/vector operations). 
    for i in range(0, m.shape[0]):
        for j in range(0, m.shape[1]):
            if m[i,j] == 0 and np.sum(star[:,j]) + np.sum(star[i,:]) == 0:
                star[i,j] = True
                cov_col[:,j] = True
    return [m, star, prime, cov_row, cov_col]

# Check if  matrix m has a noncovered zero. 
def has_noncovered_zero(m,cov_row, cov_col):
    return np.any((m == 0) & np.invert(cov_row) & np.invert(cov_col))

# These numbered steps correspond to steps in Munkres' paper. After some initial reductions, step one initiates, and there is no fixed order after that! Steps 1,2,3 can call eachother. One must show via an argument from the paper that this eventually terminates. 
def step_1(m,star,prime,cov_row, cov_col):
    if np.all(np.sum(cov_col, axis=1) == len(cov_col)):
        return [m, star, prime, cov_row, cov_col]
    while has_noncovered_zero(m,cov_row,cov_col):
        #prime the first uncovered zero. If same row has a starred zero, cover the row. 
        indices = np.where((m==0) & np.invert(cov_row) & np.invert(cov_col))
        prime_i = indices[0][0]
        prime_j = indices[1][0]
        prime[prime_i, prime_j] = True
        if not np.any(star[prime_i,:] & (m[prime_i,:] == 0)):
            [m, star, prime, cov_row, cov_col] = step_2(m, star, prime, cov_row, cov_col)
        else:
            cov_row[prime_i,:] = True
            cov_col[:, prime_j] = False
    # Proceed to step 3 when ALL zeroes are covered.
    return step_3(m,star, prime, cov_row, cov_col)

def step_2(m, star, prime, cov_row, cov_col):
    #construct a sequence of alternating starred and primed zeroes. Start with Z_0 = uncovered prime 0. This first line gets thr first such zero. Note that Z contains indices. 2 cols. 
    #Z is a 2-tuple of arrays ([rows],[cols])
    Z = np.where((m == 0) & np.invert(cov_col) & np.invert(cov_row) & prime)

    if len(Z[0]) != 1:
        print("ERROR FOUND TOO MANY UNCOVERED 0' IN STEP 2")

    continue_bool = True
    while continue_bool:
        z_i = Z[0][-1]
        z_j = Z[1][-1]
        #if just added a prime 0, add a starred 0 in the same col if it exists
        star_zero = np.where((m[:,z_j] == 0) & star[:,z_j])
        if prime[z_i, z_j] and np.any( (m[:,z_j] == 0) & star[:,z_j]):
            np.append(Z[0], z_i)
            np.append(Z[1], star_zero[0][-1] )
        if star[z_i, z_j]:
            prime_zero = np.where((m[z_i,:] == 0) & prime[z_i,:])
            np.append(Z[0], prime_zero[0][-1])
            np.append(Z[1], z_j)
        else:
            continue_bool = False
            
    star[(m == 0) & star] = False
    star[(m == 0) & prime] = True

    #Finally Erase all primes, uncover each row, and cover every column with a 0*.
    prime = np.zeros((n,n), dtype=bool)
    cov_row =  np.zeros((n,n), dtype=bool)
    col_has_star_zero = np.sum((m ==0) & star, axis=0)
    cov_col = np.tile(col_has_star_zero,(len(col_has_star_zero,1)))
    return [m, star, prime, cov_row, cov_col]

#When starting this step, ALL zeroes of m are covered. Each 0* is covered by one line.
def step_3(m,star,prime,cov_row,cov_col):
    #Get the smallest noncovered element of m, add it to each covered row and subtract it from each uncovered column. 
    tmp = np.where(np.invert(cov_row) & np.invert(cov_col))
    h = np.amin(m[tmp])
    if np.sum(cov_col[1,:]) == len(cov_col):
        return [m,star,prime,cov_row,cov_col]
    m = m + h*cov_row - h * np.invert(np.array(cov_col,dtype=bool))
    return step_1(m,star,prime,cov_row,cov_col)

def assignment(m):
    #These matrices keep track of whether each element is starred or primed. They start out as false. 
    starred =np.zeros((n,n), dtype=bool)
    primed = np.zeros((n,n), dtype=bool)

    # boolean matrices initialized to false. 
    covered_cols = np.zeros((n,n), dtype=bool)
    covered_rows =  np.zeros((n,n), dtype=bool)

    # Subtract min of each col from that col. Do the same for each row. Now every row and col has atleast one zero.
    m = m - np.tile(np.amin(m,1).reshape(len(m),1), (1,len(m)))
    m = m - np.tile(np.amin(m,0).reshape(1,len(m)), (len(m),1))

    #star ind zeros and cover their cols. The starred zeroes are ind
    [x,starred, primed,covered_rows,covered_cols] = star_ind_zeros_and_cover(m,starred, primed,covered_rows,covered_cols)

    [x,starred, primed,covered_rows,covered_cols] = step_1(x, starred, primed, covered_rows, covered_cols)
    return starred

#INITIALIZE cost, starred, primed,matrices and boolean vecs for covering.

n = 4

ORIGINAL_MATRIX = np.random.randint(0,high=10,size=(n,n))
result = assignment(ORIGINAL_MATRIX)
print(ORIGINAL_MATRIX)

MIN = 999999
counter = 0
for perm in itertools.permutations(range(0,n)):
    counter = counter + 1
    tmp = 0
    if counter > 1000000:
        break
    for ind in range(0,len(perm)):
        tmp = tmp + ORIGINAL_MATRIX[ind, perm[ind]]
    MIN = min(tmp, MIN)
print("The min calculated sum is {}".format(MIN))
print("The sum from the algorithm is {}.".format(np.sum(ORIGINAL_MATRIX[result])))
