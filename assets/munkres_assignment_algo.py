# @Author: Varoon Pazhyanur <varoon>
# @Date:   5 November 2017
# @Filename: munkres_assignment_algo.py
# @Last modified by:   varoon
# @Last modified time: 23 February 2019
import random as r
import numpy as np

'''
Algorithm in O(n^3) to solve the assignment problem as described in
"Algorithms for the Assignment and Transportation Problems" by James Munkres in the Journal of the Society of Industrial and Applied Math March 1957.

Problem: Given n people, n jobs and a rating matrix (r_{i,j}) describing a rating of person i in job j, assign each person to exactly 1 job to maximize the total rating. ie Find n independent (do not share a row or column) elements in the matrix with maximal sum. This is equivalent to minimizing that sum with a nonnegative matrix. This script runs through the algorithm with a random nxn matrix and gets the min sum. It also compares the output to a large number of options generated in the naive method.

'''
# The inputs and outputs for these functions are all essentially the 'state' of the marked matrix. 
# Get a set of independent starred zeroes and cover their cols. This is always possible because when this runs, each row and each col must have atleast one 0. 
def star_ind_zeros_and_cover(m, star,prime, cov_row, cov_col):
    #Note: This can probably be done more efficiently (but probably not with numpy/vector operations). TODO: This is basically the nonattacking rook problem from chess. 
    for i in range(0, m.shape[0]):
        for j in range(0, m.shape[1]):
            if m[i,j] == 0 and (np.sum(star[:,j]) + np.sum(star[i,:]) == 0):
                star[i,j] = True
                cov_col[:,j] = True
    return [m, star, prime, cov_row, cov_col]

# These numbered steps correspond to steps in Munkres' paper. After some initial reductions, step one initiates, and there is no fixed order after that! Steps 1,2,3 can call eachother. One must show via an argument from the paper that this eventually terminates. 
def step_1(m,star,prime,cov_row, cov_col):
    if np.sum(star)==len(star):
        return [m, star, prime, cov_row, cov_col]
    while np.any((m == 0) & np.invert(np.array(cov_row, dtype=bool)) & np.invert(np.array(cov_col, dtype=bool))):#while has noncovered zero 
        indices = np.where((m==0) & np.invert(np.array(cov_row,dtype=bool)) &np.invert(np.array(cov_col,dtype=bool)))
        prime_i = indices[0][0]
        prime_j = indices[1][0]
        prime[prime_i, prime_j] = True
        if not np.any(star[prime_i,:] & (m[prime_i,:] == 0)):
            [m, star, prime, cov_row, cov_col] = step_2(m, star, prime, cov_row, cov_col) 
            if np.sum(star)==len(star):
                return [m, star, prime, cov_row, cov_col]
        else:
            cov_row[prime_i,:] = True
            cov_col[:, np.where(star[prime_i,:] & (m[prime_i,:] == 0))[-1]] = False
                   
    return step_3(m,star, prime, cov_row, cov_col)

def step_2(m, star, prime, cov_row, cov_col):
    #construct a sequence of alternating starred and primed zeroes. Start with Z_0 = uncovered prime 0. Each col of Z corresponds to an index pair. 
    Z = np.array(np.where((m == 0) & np.invert(np.array(cov_col,dtype=bool)) & np.invert(cov_row) & prime))

    continue_bool = True
    while continue_bool:
        z_i = Z[0,-1]
        z_j = Z[1,-1]
        #Add a starred 0 in the same col  as the primed zero if it exists
        star_zero = np.where((m[:,z_j] == 0) & star[:,z_j])
        if np.any( (m[:,z_j] == 0) & star[:,z_j]):
            Z = np.append(Z, [[star_zero[0][-1]],[z_j]],1)
            z_i = Z[0,-1]
            z_j = Z[1,-1]
            prime_zero = np.where((m[z_i,:] == 0) & prime[z_i,:])
            Z = np.append(Z, [[z_i],[prime_zero[0][-1]]],1) 
        else:
            continue_bool = False   
    for j in range(Z.shape[1]):
        star[Z[0,j],Z[1,j]] = j%2==0

    #Finally Erase all primes, uncover each row, and cover every column with a 0*.
    n = m.shape[0]
    prime = np.zeros((n,n), dtype=bool)
    cov_row =  np.zeros((n,n), dtype=bool)
    cov_col[:,np.where(np.sum((m ==0) & star, axis=0)==1)[0]] = True
    return [m, star, prime, cov_row, cov_col]

#When starting this step, ALL zeroes of m are covered. Each 0* is covered by one line.

def step_3(m,star,prime,cov_row,cov_col):
    if np.sum(star)==len(star):
        return [m, star, prime, cov_row, cov_col]
    #Get the smallest noncovered element of m, add it to each covered row and subtract it from each uncovered column. 
    h = np.amin(m[np.invert(np.array(cov_row,dtype=bool)) & np.invert(np.array(cov_col,dtype=bool))])
    if np.sum(cov_col[0,:]) == len(cov_col):
        return [m,star,prime,cov_row,cov_col]
    t = h * cov_row
    s = h * np.invert(np.array(cov_col,dtype=bool))
    m = m + t-s
    return step_1(m,star,prime,cov_row,cov_col)

#INPUT: m - a square np array. 
#OUTPUT: A square boolean array that has a 1 at (i,j) if and only if person i is assigned to job j
def assignment(m):
    #Using the terms from Munkres' paper, these matrices keep track of whether each element is starred or primed and whether the row/col of the element is covered. They start out as false. 
    n = m.shape[0]
    star =np.zeros((n,n), dtype=bool)
    prime = np.zeros((n,n), dtype=bool)
    cov_col = np.zeros((n,n), dtype=bool)
    cov_row =  np.zeros((n,n), dtype=bool)

    # Subtract min of each col from that col. Do the same for each row. Now every row and col has atleast one zero.
    m = m - np.tile(np.amin(m,1).reshape(len(m),1), (1,len(m)))
    m = m - np.tile(np.amin(m,0).reshape(1,len(m)), (len(m),1))

    [x,star, prime,cov_row,cov_col] = star_ind_zeros_and_cover(m,star, prime,cov_row,cov_col)

    [x,star, prime,cov_row,cov_col] = step_1(x, star, prime, cov_row, cov_col)
    return star

if __name__ == "__main__":
    k = 200

    r = np.random.randint(0,20,size=(k,k))
    np.savetxt('assignment_large.txt', r, fmt = '%d')
    np.savetxt('assignment_large_solution.txt', assignment(r), fmt = '%d')
