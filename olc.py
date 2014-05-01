#!/usr/bin/env python

match_award = 1
mismatch_penalty = -2
gap_penalty = -2
      
#Creates empty matrix with all zeros
def zeros(shape):
    retval = []
    for x in range(shape[0]):
        retval.append([])
        for y in range(shape[1]):
            retval[-1].append(0)
    return retval

#No substituition matrix, just simple linear gap penalty model
def match_score(alpha, beta):
    if alpha == beta:
        return match_award
    elif alpha == '-' or beta == '-':
        return gap_penalty
    else:
        return mismatch_penalty

def nw(seq1, seq2):
    global max_score
   
    m = len(seq1)
    n = len(seq2) # lengths of two sequences
    
    # Generate DP (Dynamic Programmed) table and traceback path pointer matrix
    score = zeros((n+1, m+1)) # the DP table
    for i in range(n+1) :
      score[i][0] = 0
    
    for j in range(m+1) :
      score[0][j] = 0
    
    #Traceback matrix
    pointer = zeros((n+1, m+1)) # to store the traceback path
    for i in range(n+1) :
      pointer[i][0] = 1
    for j in range(m+1) :
      pointer[0][j] = 2
    
    # Calculate DP table and mark pointers
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            score_diagonal = score[i-1][j-1] + match_score(seq1[j-1], seq2[i-1])
            score_up = score[i][j-1] + gap_penalty
	    score_left = score[i-1][j] + gap_penalty
	    score[i][j] = max(score_left, score_up, score_diagonal)
	    
	    if score[i][j] == score_diagonal :
                pointer[i][j] = 3 # 3 means trace diagonal
	    elif score[i][j] == score_up :
                pointer[i][j] = 2 # 2 means trace left
	    elif score[i][j] == score_left :
                pointer[i][j] = 1 # 1 means trace up
            
    #Finding the right-most match which represents a longest overlap
    #Note that .index() will find the index of the first item in the list that matches, 
    #so if you had several identical "max" values, the index returned would be the one for the first.
    max_i = -200
    for ii in range(n+1) :
      if score[ii][-1] >= max_i :
	max_i = score[ii][-1]
	i = ii
#    print max_i

    prei = i
    prej = j
    #Traceback, follow pointers in the traceback matrix
    align1, align2 = '', '' # initial sequences
    while 1 : 
      if pointer[i][j] == 3 :
	align1 = seq1[j-1] + align1
	align2 = seq2[i-1] + align2
	i -= 1
	j -= 1
      elif pointer[i][j] == 2 : #2 means trace left
	align2 = '-' + align2
	align1 = seq1[j-1] + align1
	j -= 1
      elif  pointer[i][j] == 1 : # 1 means trace up
	align2 = seq2[i-1] + align2
	align1 = '-' + align1
	i -= 1
#      print align1, align2, j, i
      if (i == 0 or j == 0) : break
    
    return (align1, align2, prej, j, prei, i, max_i)
#    print align1 + '\n' + align2

#print seq11
#print seq22    
#nw(seq11, seq22)
#nw(seq22,seq11)
