# @Author: Varoon Pazhyanur <varoon>
# @Date:   05-11-2017
# @Filename: munkres_assignment_algo.py
# @Last modified by:   varoon
# @Last modified time: Jan-11-2018
import random as r
import copy
import itertools
'''
Algorithm in O(n^3) to solve the assignment problem as described in
"Algorithms for the Assignment and Transportation Problems" by James Munkres
in the Journal of the Society of Industrial and Applied Math March 1957.

Problem: Given n people, n jobs and a rating matrix (r_{i,j}) describing a rating of person i
in job j, assign each person to exactly 1 job to maximize the total rating.

ie Find n independent (do not share a row or column) elements in the matrix with
maximal sum. Equivelent to minimizing that sum with a nonnegative matrix.

Runs through the algorithm with a random nxn matrix and gets the min sum. Also compares to up to
the max sum of 1000000 sets of independent elements.
'''

def print_mat_alt(m):

    for i in range(0,len(m)):
        for j in range(0,len(m[0])):
            print("{0: <8}".format(str(m[i][j])), end = "\t")
        print()

def print_mat(m,star,prime,cov_row, cov_col):

    for i in range(0,len(m)):
        for j in range(0,len(m[0])):
            st = ""
            pr = ""
            row = ""
            col = ""
            if star[i][j]:
                st = "*"
            if prime[i][j]:
                pr = "\'"
            if cov_row[i]:
                col = "R"
            if cov_col[j]:
                row = " C"
            res = str(m[i][j])+st+pr+row+col
            print("{0: <8}".format(res), end = "\t")
        print()

#takes a matrix m and a row index i and subtracts the min element in that row from each elem in the row
def subtr_min_from_row(m,i):
    MIN = 99999999
    for j in range(0,len(m[0])):
        MIN = min(MIN, m[i][j])
    for k in range(0,len(m[0])):
        m[i][k] = m[i][k] - MIN
    return m

#takes a matrix m and a col index j and subtracts the min element in that col from each elem in the col

def subtr_min_from_col(m,j):
    MIN = 99999999
    for i in range(0,len(m)):
        MIN = min(MIN, m[i][j])
    for k in range(0,len(m)):
        m[k][j] = m[k][j] - MIN
    return m

#Subtract min element of each from from that row. Repeat for each col. Return resulting matrix.

def subtr_min_from_ea_row_col(m):

    for i in range(0, len(m)):
        m = subtr_min_from_row(m, i)

    for j in range(0, len(m[0])):
        m = subtr_min_from_col(m, j)

    return m

#get a set of ind starred zeroes and cover their cols
def star_ind_zeros_and_cover(m, star,prime, cov_row, cov_col):
    for i in range(0, len(m)):
        for j in range(0, len(m[0])):

            if m[i][j] == 0 and (sum([int(star[i][k]) for k in range(0, len(star[0]))]) + sum([int(star[k][j]) for k in range(0, len(star))]))==0:
                star[i][j] = True
                cov_col[j] = True
    print("Starred and covered independent zeros. There are {} covered rows, {} covered cols, {} primes, {} stars, "
          "and the sum is {}".format(sum(cov_row), sum(cov_col), sum(sum(prime, [])), sum(sum(star, [])), sum(sum(m, []))))
    print_mat(m, star, prime, cov_row, cov_col)
    print()
    return [m, star, prime, cov_row, cov_col]

#Check of matrix has a noncovered zero.
def has_noncovered_zero(m,cov_row, cov_col):
    for i in range(0,len(m)):
        for j in range(0,len(m[0])):
            if m[i][j] == 0 and not cov_row[i] and not cov_col[j]:
                return True
    return False

def step_1(m,star,prime,cov_row, cov_col):
    global num_steps_1
    num_steps_1 = num_steps_1 + 1
    while has_noncovered_zero(m,cov_row,cov_col):
        #prime the first uncovered zero. If same row has a starred zero, cover the row.
        count = 0
        for i in range(0,len(m)):
            for j in range(0,len(m[0])):
                if count == 0 and m[i][j] == 0 and not cov_row[i] and not cov_col[j]:
                    count = count+1
                    prime[i][j] = True
                    if prime[i][j] and star[i][j]:
                        print("ERROR STARRED PRIMED 0")
                        print_mat(m, star, prime, cov_row, cov_col)
                    #If has starred zero in same row(starred zeros are ind by construction), redo coverings.
                    row_has_starred_zero = False
                    for k in range(0,len(m[0])):
                        if m[i][k] == 0 and star[i][k]:
                            row_has_starred_zero = True
                            cov_row[i] = True
                            cov_col[k] = False
                            print("End of Step 1. There are {} covered rows, {} covered cols, {} primes, {} stars, "
                                  "and the sum is {}".format(sum(cov_row), sum(cov_col), sum(sum(prime, [])), sum(sum(star, [])),sum(sum(m, []))))
                            print_mat(m, star, prime, cov_row, cov_col)
                            print()
                    if not row_has_starred_zero:
                        [m, star, prime, cov_row, cov_col] = step_2(m, star, prime, cov_row, cov_col)
                        if sum(cov_col) == len(cov_col):
                            return [m, star, prime, cov_row, cov_col]

    # Proceed to step 3 when ALL zeroes are covered. May need to repeat steps 1 and 2 several times.
    return step_3(m,star, prime, cov_row, cov_col)

def step_2(m, star, prime, cov_row, cov_col):
    global num_steps_2
    num_steps_2 = num_steps_2 + 1
    #construct a sequence of alternating starred and primed zeroes. Start with Z_0 = uncovered prime 0
    Z = []
    for x in range(0,len(m)):
        for y in range(0,len(m[0])):
            if m[x][y] == 0 and prime[x][y] and not cov_col[y] and not cov_row[x]:
                Z.append([x,y])
                break

    if len(Z) is not 1:
        print("ERROR FOUND TOO MANY UNCOVERED 0' IN STEP 2")

    continue_bool = True
    while continue_bool:
        prev_z = Z[len(Z) - 1]
        print(len(Z))
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
    prime = [[False for i in range(0,len(m))]for j in range(0, len(m[0]))]
    cov_row = [False for i in range(0, len(m))]
    for j in range(0,len(cov_col)):
        col_has_star_zero = False
        for i in range(0,len(m)):
            if m[i][j] == 0 and star[i][j]:
                col_has_star_zero = True
        cov_col[j] = col_has_star_zero
    print("End of Step 2. There are {} covered rows, {} covered cols, {} primes, {} stars, and the sum is {}".format(sum(cov_row), sum(cov_col), sum(sum(prime, [])), sum(sum(star, [])), sum(sum(m, []))))
    print_mat(m, star, prime, cov_row, cov_col)
    print()
    return [m, star, prime, cov_row, cov_col]

def step_3(m,star,prime,cov_row,cov_col):
    global num_steps_3
    num_steps_3 = num_steps_3 + 1
    #When starting this step, ALL zeroes of m are covered. Each 0* is covered by one line.
    #Get the smallest noncovered element
    #print_mat(m, star, prime, cov_row, cov_col)
    h = 99999999
    for i in range(0,len(m)):
        for j in range(0,len(m[0])):
            if not cov_row[i] and not cov_col[j] :
                h = min(h,m[i][j])
    if sum(cov_col) == len(cov_col):
        return [m,star,prime,cov_row,cov_col]
    #add h to each covered row and subtract h from each uncovered col
    for a in range(0,len(m)):
        for b in range(0,len(m[0])):
            if cov_row[a]:
                m[a][b] = m[a][b]+h
            if not cov_col[b]:
                m[a][b] = m[a][b]-h
    print("End of Step 3. There are {} covered rows, {} covered cols, {} primes, {} stars, and the sum is {}".format(sum(cov_row), sum(cov_col), sum(sum(prime,[])), sum(sum(star,[])), sum(sum(m,[])) ))
    print_mat(m, star, prime, cov_row, cov_col)
    print()
    return step_1(m,star,prime,cov_row,cov_col)

#INITIALIZE cost, starred, primed,matrices and boolean vecs for covering.
num_steps_1 = 0
num_steps_2 = 0
num_steps_3 = 0
num_rows = 10
num_cols = 10

ORIGINAL_MATRIX = [[r.randint(0,11) for i in range(0,num_cols)] for j in range(0,num_rows)]
starred = [[False for i in range(0,num_cols)]for j in range(0,num_rows)]
primed = [[False for i in range(0,num_cols)]for j in range(0,num_rows)]
covered_cols = [False for i in range(0,num_cols)]
covered_rows = [False for j in range(0,num_rows)]
x = copy.deepcopy(ORIGINAL_MATRIX)

print("Original Matrix")
print_mat(ORIGINAL_MATRIX,starred, primed,covered_rows,covered_cols)
print()

#Subtract min from each row and col. Now every row and col has atleast one zero.
x = subtr_min_from_ea_row_col(x)
print_mat(x,starred, primed,covered_rows,covered_cols)
print()

#star ind zeros and cover their cols. The starred zeroes are ind
[x,starred, primed,covered_rows,covered_cols] = star_ind_zeros_and_cover(x,starred, primed,covered_rows,covered_cols)

[x,starred, primed,covered_rows,covered_cols] = step_1(x, starred, primed, covered_rows, covered_cols)
print("FINISHED!")
print_mat(x,starred, primed,covered_rows,covered_cols)

print()
sum = 0
for i in range(0, len(x)):
    for j in range(0, len(x[0])):
        if starred[i][j]:
            sum = sum + ORIGINAL_MATRIX[i][j]

print("The sum from the algorithm is {}".format(sum))


#Now loop through possibilities to find max with brute force method (cap at 10k.
MIN = 999999
counter = 0
for perm in itertools.permutations(range(0,len(x))):
    counter = counter + 1
    tmp = 0
    if counter > 1000000:
        break
    for ind in range(0,len(perm)):
        tmp = tmp + ORIGINAL_MATRIX[ind][perm[ind]]
    MIN = min(tmp, MIN)
print("The min calculated sum is {}".format(MIN))
print("It took {}, {}, and {} loops of steps 1,2,3 respectively".format(num_steps_1, num_steps_2, num_steps_3))
