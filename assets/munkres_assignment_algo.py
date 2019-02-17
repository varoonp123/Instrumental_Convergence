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

Problem: Given n people, n jobs and a rating matrix (r_{i,j}) describing a rating of person i in job j, assign each person to exactly 1 job to maximize the total rating.

ie Find n independent (do not share a row or column) elements in the matrix with maximal sum. This is equivalent to minimizing that sum with a nonnegative matrix.

This script runs through the algorithm with a random nxn matrix and gets the min sum. It also compares the output to a large number of options generated in the naive method.
'''

#takes a matrix m and a row index i and subtracts the min element in that row from each elem in the row

def subtr_min_from_row(m,i):
    m[i,:]= m[i,:] - np.amin(m,1)[i]

#takes a matrix m and a col index j and subtracts the min element in that col from each elem in the col

def subtr_min_from_col(m,j):
    m[:,j]= m[:,j] - np.amin(m,0)[j]
    return m

#Subtract min element of each from from that row. Repeat for each col. Return resulting matrix.
# TODO: I bet this can be done better with a vector operation, but I couldn't figure it out. 
def subtr_min_from_ea_row_col(m):

    for i in range(0, m.shape[0]):
        m = subtr_min_from_row(m, i)
    for j in range(0, m.shape[1]):
        m = subtr_min_from_col(m, j)
    return m

#get a set of ind starred zeroes and cover their cols. This is always possible because when this runs, each row and each col must have atleast one 0. 
def star_ind_zeros_and_cover(m, star,prime, cov_row, cov_col):
    #Note: This can probably be done more efficiently (but probably not with numpy/vector operations). 
    for i in range(0, m.shape[0]):
        for j in range(0, m.shape[1]):
            if m[i,j] == 0 and np.sum(star[:,j]) + np.sum(star[i,:]) == 0:
                star[i,j] = True
                cov_col[j] = True
    return [m, star, prime, cov_row, cov_col]

# Check if  matrix m has a noncovered zero. 
def has_noncovered_zero(m,cov_row, cov_col):
    return np.sum((m == 0)*np.matmult(np.invert(cov_row).reshape(cov_row.size,1), np.invert(cov_col).reshape(1,cov_col.size))) > 0

# These numbered steps correspond to steps in Munkres' paper. After some initial reductions, step one initiates, and there is no fixed order after that! Steps 1,2,3 can call eachother. One must show via an argument from the paper that this eventually terminates. 
def step_1(m,star,prime,cov_row, cov_col):
    while has_noncovered_zero(m,cov_row,cov_col):
        #prime the first uncovered zero. If same row has a starred zero, cover the row.
        #note: Some of this CAN be done with vector/numpy operations, but this is an algorithm. It matters that we just find the first occurance of something, so 
        count = 0
        for i in range(0,m.shape[0]):
            for j in range(0,m.shape[1]):
                if count == 0 and m[i,j] == 0 and not cov_row[i] and not cov_col[j]:
                    count = count+1
                    prime[i,j] = True
                    if prime[i,j] and star[i,j]:
                    #If has starred zero in same row(starred zeros are ind by construction), redo coverings.
                    row_has_starred_zero = np.any(m[i,:] == 0 and star[i,:])
                    #CONTINUE HERE!
                    
                    for k in range(0,m.shape[0]):
                        if m[i,k] == 0 and star[i,k]:
                            row_has_starred_zero = True
                            cov_row[i] = True
                            cov_col[k] = False
                    if not row_has_starred_zero:
                        [m, star, prime, cov_row, cov_col] = step_2(m, star, prime, cov_row, cov_col)
                        if np.sum(cov_col) == len(cov_col):
                            return [m, star, prime, cov_row, cov_col]

    # Proceed to step 3 when ALL zeroes are covered. May need to repeat steps 1 and 2 several times.
    return step_3(m,star, prime, cov_row, cov_col)

def step_2(m, star, prime, cov_row, cov_col):
    #construct a sequence of alternating starred and primed zeroes. Start with Z_0 = uncovered prime 0. This first line gets thr first such zero. Note that Z contains indices. 2 cols. 
    Z = np.array([[x,y] for x in range(0, m.shape[0]) for y in m.shape[1] if m[x,y] ==0 and prime[x,y] and not cov_col[y] and not cov_row[x] ])[0,:]

    if Z.shape[0] != 1:
        print("ERROR FOUND TOO MANY UNCOVERED 0' IN STEP 2")

    continue_bool = True
    while continue_bool:
        prev_z = Z[len(Z) - 1,:]
        has_star_zero = False
        #if just added a prime 0, add a starred 0 in the same col if it exists
        for x in range(0,len(m)):
            if m[x][prev_z[1]]==0 and star[x][prev_z[1]] and continue_bool:
                Z.append([x,prev_z[1]])
                has_star_zero = True
        if has_star_zero:
            #Add a prime zero in the same row (It Must exist)
            prev_z = Z[len(Z) - 1]
            for y in range(0, len(m[0])):
                if m[prev_z[0]][y] == 0 and prime[prev_z[0]][y] and continue_bool:
                    Z.append([prev_z[0], y])
        else:
            continue_bool = False

    #After constructing sequence, unstar each starred 0 and star each prime 0 of the sequence, resulting in ind starred zeroes
    for n in range(0,len(Z)):
        if n % 2 == 0:
            star[Z[n][0]][Z[n][1]] = True
        else:
            star[Z[n][0]][Z[n][1]] = False

    #Finally Erase all primes, uncover each row, and cover every column with a 0*.
    prime = np.zeros((n,n), dtype=bool)
    cover_row =  np.zeros(n, dtype=bool)
    cov_col = np.sum(np.any(m ==0 and star == 1), axis = 1)


    for j in range(0,len(cov_col)):
        col_has_star_zero = False
        for i in range(0,len(m)):
            if m[i][j] == 0 and star[i][j]:
                col_has_star_zero = True
        cov_col[j] = col_has_star_zero
    return [m, star, prime, cov_row, cov_col]

#When starting this step, ALL zeroes of m are covered. Each 0* is covered by one line.
def step_3(m,star,prime,cov_row,cov_col):
    #Get the smallest noncovered element of m
    h = np.amin(m[np.matmult(np.invert(cov_row).reshape(cov_row.size,1), np.invert(cov_col).reshape(1,cov_col.size))])

    if np.sum(cov_col) == np.size(cov_col):
        return [m,star,prime,cov_row,cov_col]

    #add h to each covered row and subtract h from each uncovered col

    m = h*np.matmul(cov_row.reshape(cov_row.size, 1), np.ones(1, cov_row.size) - h * np.matmul(np.ones(cov_row.size,1), cov_col.reshape(1,cov_col.size))
    return step_1(m,star,prime,cov_row,cov_col)

def assignment(m):
    #These matrices keep track of whether each element is starred or primed. They start out as false. 
    starred =np.zeros((n,n), dtype=bool)
    primed = np.zeros((n,n), dtype=bool)

    # row boolean vectors initialized to false. 
    covered_cols = np.zeros(n, dtype=bool)
    covered_rows =  np.zeros(n, dtype=bool)

    x = copy.deepcopy(ORIGINAL_MATRIX)

    #Subtract min from each row and col. Now every row and col has atleast one zero.
    x = subtr_min_from_ea_row_col(x)

    #star ind zeros and cover their cols. The starred zeroes are ind
    [x,starred, primed,covered_rows,covered_cols] = star_ind_zeros_and_cover(x,starred, primed,covered_rows,covered_cols)

    [x,starred, primed,covered_rows,covered_cols] = step_1(x, starred, primed, covered_rows, covered_cols)
    return primed
#INITIALIZE cost, starred, primed,matrices and boolean vecs for covering.

n = 1000

ORIGINAL_MATRIX = np.random.randint(0,high=10,size=(n,n))
