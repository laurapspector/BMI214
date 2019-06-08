'''
October 16, 2018
This file takes as input an input file as specified in the Project 1 guidelines
and generates an output file listing the optimal alignment score followed by
the optimal alignments with gaps indicated by '_'. The program can align
sequences using both 'global with no end gaps' or 'local alignment' rules.

Usage: python3 align.py input_file output_file
'''

import sys

#### ------ USEFUL FUNCTIONS ------- ####
def fuzzy_equals(a, b):
    """
    Checks if two floating point numbers are equivalent.

    Inputs:
        a = a number
        b = another number
    Returns:
        True or False
    """
    epsilon = 10**(-6) 
    return (abs(a - b) < epsilon)
    

#### ------- CLASSES ------- ####

class MatchMatrix(object):
    """
    Match matrix class stores the scores of matches in a data structure
    Data structure used here is a dictionary
    """
    def __init__(self):

        # Initialize a match matrix object
        self.matchMatrix = {}

    def fill_match_matrix(self, inputList):
        """
        Stores tuples of each letter pair as keys in a dictionary with the 
        score corresponding to that position in the match matrix as the value.

        Inputs:
           inputList = row-delimited list of match matrix entries parsed from input file
        """
        # Empty match matrix in case it already exists
        self.matchMatrix = {}

        for i in inputList:    # i is a single row from inputList
            i = i.strip().split()
            alpha_A_letter = i[2]   # = match matrix row
            alpha_B_letter = i[3]   # = match matrix column
            score = float(i[4])
            self.matchMatrix[ (alpha_A_letter, alpha_B_letter)] = score

    def set_score(self, a, b, score):
        """
        Updates or adds a score for a specified match

        Inputs:
           a = the character from sequence A
           b = the character from sequence B
           score = the score to set it for
        """
        self.matchMatrix[ (a,b) ] = score


    def get_score(self, a, b):
        """
        Returns the score for a particular match, where a is the
        character from sequence A and b is from sequence B.

        Inputs:
           a = the character from sequence A
           b = the character from sequence B
        Returns:
           the score of that match
        """
        return self.matchMatrix[ (a,b) ]


class ScoreMatrixEntry(object):
    '''
    Object to store an entry in the ScoreMatrix object. Each score matrix
    (M, Ix, and Iy) has a single ScoreMatrixEntry object at each (i,j) entry.
    '''

    def __init__(self, row, col, name=''):
        self.row = row
        self.col = col
        self.score = 0.0
        self.pointedObjects = [[]] # will contain the pointers for this entry
        self.name = name # identifier for the score matrix in which it exists

class ScoreMatrix(object):
    """
    Object to store a score matrix, which generated during the alignment process. 
    The score matrix consists of a 2-D array of ScoreMatrixEntry objects that are 
    updated during alignment and used to output the maximum alignment.
    """

    def __init__(self, name, nrow=None, ncol=None):
        self.name = name # identifier for the score matrix - Ix, Iy, or M
        self.score_matrix = []
        self.nrow = nrow
        self.ncol = ncol

        # Allow for creation of a score matrix of pre-defined size for 
        # testing unit tests
        if self.nrow and self.ncol != None:
            self.initialize_matrix(self.nrow, self.ncol, self.name)

    def initialize_matrix(self, nrow, ncol, name):
        '''
        Fills the ScoreMatrix object's 2-D array with ScoreMatrixEntry objects
        and updates each object with the matrix name and its indices.
        '''
        self.nrow = nrow
        self.ncol = ncol    

        # Initialize 2-D array
        for i in range(self.nrow): 
            self.score_matrix.append([None]*self.ncol)
        
        # Fill 2-D array with ScoreMatrixEntry objects
        for i in range(self.nrow):
            for j in range(self.ncol):
                self.score_matrix[i][j] = ScoreMatrixEntry(i, j)
                self.score_matrix[i][j].name = name

    def get_score(self, row, col):
        '''
        Returns the score of the object at a given position (0-indexed) 
        in the score matrix.
        Inputs:
            row = row
            col = column
        Returns:
            score
        '''
        return self.score_matrix[row][col].score
        
    def set_score(self, row, col, score):    
        '''
        Sets the score of the object at a given position (0-indexed)
        Inputs:
            row = row
            col = column
        '''     
        self.score_matrix[row][col].score = score

    def get_entry_object(self, row, col):
        '''
        Returns the ScoreMatrixEntry object at a given position in the score matrix.
        Inputs:
            row = row
            col = column
        Returns:
            ScoreMatrixEntry object
        '''
        return self.score_matrix[row][col]

    def get_pointers(self, row, col):
        """
        Returns the indices (0-indexed) of the entries that are pointed to and the 
        name of the matrix to which they belong. Works by iterating through a list
        of pointed-to ScoreMatrixEntry objects.
        Inputs:
            row = row
            col = column
        Returns:
            A list of tuples containing (name, row, column) of the entry pointed to.
        """
        pointer_list = []
        for i in self.score_matrix[row][col].pointedObjects:
            pointer_list.append( (i.name, i.row, i.col) )
        return pointer_list

    def set_pointers(self, row, col, entryObjects):
        '''
        Sets the pointers for a given ScoreMatrixEntry object.
        Inputs:
            row = row
            col = column
            entryObjects = list of ScoreMatrixEntry objects that are pointed to
        '''
        self.score_matrix[row][col].pointedObjects = entryObjects


    def print_scores(self):
        """
        Returns a nicely formatted string containing the scores in the score matrix
        that can be printed from the main methods.

        When printed, looks like:
        M=
            0.0, 0.0, 0.0, 0.0, 0.0
            0.0, 1.0, 0.0, 0.0, 0.0
            0.0, 1.0, 1.0, 1.0, 1.0
            0.0, 0.0, 1.0, 1.0, 1.0
            0.0, 0.0, 2.0, 2.0, 1.0
            0.0, 0.0, 1.0, 2.0, 3.0

        """
        score_string = self.name + "=" + '\n'

        for i in self.score_matrix:
            score_string = score_string + '\t' + str(i[0].score)
            for j in i[1:]:               
                score_string = score_string + ', ' + str(j.score)
            score_string += '\n'

        return score_string
  

    def print_pointers(self):
        """
        Returns a nicely formatted string containing the pointers for each entry 
        in the score matrix that can be printed from the main methods.
        
        Prints:
            Prints by row a comma(column)-delimited set of pointers corresponding
            to the pointers for that object in that (row,column).
        """
        matrix_pointers = ''
        pointers = []

        for i in range(self.nrow):
            matrix_pointers = matrix_pointers + 'Row ' + str(i) + ':'
            for j in range(self.ncol):
                pointers = self.get_pointers(i, j)
                for k in pointers:
                    matrix_pointers = matrix_pointers + ' ' + k[0] + '(' + str(k[1]) + ',' + str(k[2]) + ')'
                matrix_pointers += ','
            matrix_pointers += '\n'  
        return matrix_pointers

class AlignmentParameters(object):
    '''
    Object to score alignment parameters.
    '''    

    def __init__(self):
        # default values for variables that are filled in by reading
        # the input alignment file
        self.seq_a = ""
        self.seq_b = ""
        self.global_alignment = False
        self.dx = 0
        self.ex = 0
        self.dy = 0
        self.ey = 0
        self.alphabet_a = "" 
        self.alphabet_b = ""
        self.len_alphabet_a = 0
        self.len_alphabet_b = 0
        self.match_matrix = MatchMatrix()

    def load_params_from_file(self, input_file): 
        """
        Reads the parameters from an input file and stores in the object

        Input:
           input_file = specially formatted alignment input file
        """

        f=open(input_file)

        inputLines = [line.strip() for line in f if line.strip()]
        self.seq_a, self.seq_b = inputLines[0:2]

        # update the global_alignment setting (flag is 0, set as True)
        if int(inputLines[2]) == 0:
            self.global_alignment = True

        penalties = [float(i) for i in inputLines[3].split(' ')]
        self.dx, self.ex, self.dy, self.ey = penalties
        self.len_alphabet_a = int(inputLines[4])
        self.alphabet_a = inputLines[5]
        self.len_alphabet_b = int(inputLines[6])
        self.alphabet_b = inputLines[7]

        # fill the match matrix
        self.match_matrix.fill_match_matrix(inputLines[8:])
        
        f.close()


class AlignmentAlgorithm(object):
    '''
    Object that contains algorithms for filling and updating the score matrices.
    '''

    def __init__(self):

        self.m_matrix = ScoreMatrix("M")
        self.ix_matrix = ScoreMatrix("Ix")
        self.iy_matrix = ScoreMatrix("Iy")

    def populate_score_matrices(self, align_params):
        """
        Method to populate the score matrices.
        Should call update(i,j) for each entry in the score matrices

        Input:
           align_params = an AlignmentParameters object with the input
        """
        # Add a row and column for the boundary conditions
        nrow = len(align_params.seq_a) + 1
        ncol = len(align_params.seq_b) + 1

        self.m_matrix.initialize_matrix(nrow, ncol, 'M')
        self.ix_matrix.initialize_matrix(nrow, ncol, 'Ix')
        self.iy_matrix.initialize_matrix(nrow, ncol, 'Iy')

        for i in range(nrow):
            for j in range(ncol):

                # Assign pointers to the boundary indices; default score is 0.0
                if i == 0 and j == 0:
                    self.m_matrix.set_pointers(i, j, [self.m_matrix.get_entry_object(i,j)])
                    self.ix_matrix.set_pointers(i, j, [self.m_matrix.get_entry_object(i,j)])
                    self.iy_matrix.set_pointers(i, j, [self.m_matrix.get_entry_object(i,j)])
                elif i == 0 and not j == 0:
                    self.m_matrix.set_pointers(i, j, [self.m_matrix.get_entry_object(i,j-1)])
                    self.ix_matrix.set_pointers(i, j, [self.m_matrix.get_entry_object(i,j-1)])
                    self.iy_matrix.set_pointers(i, j, [self.m_matrix.get_entry_object(i,j-1)])
                elif j == 0 and not i == 0:
                    self.m_matrix.set_pointers(i, j, [self.m_matrix.get_entry_object(i-1,j)])
                    self.ix_matrix.set_pointers(i, j, [self.m_matrix.get_entry_object(i-1,j)])
                    self.iy_matrix.set_pointers(i, j, [self.m_matrix.get_entry_object(i-1,j)])
                # Fill all other entries:
                else:
                    self.update(i, j, align_params)


    def update(self, row, col, align_params):
        """
        Method to update the matrices at a given row and column index (0-indexed!).

        Input:
           row = the row index to update
           col = the column index to update
           align_params = an AlignmentParameters object with the input
        """
        self.update_m(row, col, align_params)
        self.update_ix(row, col, align_params)
        self.update_iy(row, col, align_params)

    def update_m(self, row, col, align_params):
        '''
        Method to update the entries in the M matrix
        Inputs:
            row = row to update
            col = column to update
            align_params = an AlignmentParameters object with the input
        '''
        
        match_matrix = align_params.match_matrix

        seq_a = align_params.seq_a
        seq_b = align_params.seq_b

        # Identify which letters to align using indexing for position in sequence
        # Use these to get the score for aligning those letters from match matrix
        letter_a = align_params.seq_a[row-1]
        letter_b = align_params.seq_b[col-1]
        Sij = match_matrix.get_score(letter_a, letter_b)

        # Get all scores and max score used to assign score and this entry
        all_scores = [self.m_matrix.get_score(row-1, col-1) + Sij, \
            self.ix_matrix.get_score(row-1, col-1) + Sij, \
            self.iy_matrix.get_score(row-1, col-1) + Sij]
        max_score = max(all_scores)
        
        # Condition to prevent negative scores in score matrix for local alignment
        if max_score < 0.0 and not align_params.global_alignment:
            max_score = 0.0
        self.m_matrix.set_score(row, col, max_score)
        
        # Assign pointers for the updated entry by pairing max score
        # with matrix or matrices from which it came
        pointerObjects = []
        matrix_list = [self.m_matrix, self.ix_matrix, self.iy_matrix]
        max_matrices_by_index = [i for i, j in enumerate(all_scores) if fuzzy_equals(j,max_score)]
        pointerObjects = [matrix_list[k].get_entry_object(row-1, col-1) for k in max_matrices_by_index]
        
        self.m_matrix.set_pointers(row, col, pointerObjects)
        

    def update_ix(self, row, col, align_params):
        '''
        Method to update the entries in the Ix matrix
        Inputs:
            row = row to update
            col = column to update
            align_params = an AlignmentParameters object with the input
        '''
        
        dy = align_params.dy
        ey = align_params.ey
        # Get all scores and max score used to assign score and this entry
        all_scores = [self.m_matrix.get_score(row-1, col) - dy, \
            self.ix_matrix.get_score(row-1, col) - ey]
        max_score = max(all_scores)

        # Condition to prevent negative scores in score matrix for local alignment
        if max_score < 0.0 and not align_params.global_alignment:
        	max_score = 0.0
        self.ix_matrix.set_score(row, col, max_score)

        # Assign pointers for the updated entry by pairing max score
        # with matrix or matrices from which it came
        pointerObjects = []
        matrix_list = [self.m_matrix, self.ix_matrix]
        max_matrices_by_index = [i for i, j in enumerate(all_scores) if fuzzy_equals(j,max_score)]
        pointerObjects = [matrix_list[k].get_entry_object(row-1, col) for k in max_matrices_by_index]

        self.ix_matrix.set_pointers(row, col, pointerObjects)

    def update_iy(self, row, col, align_params):
        '''
        Method to update the entries in the Iy matrix
        Inputs:
            row = row to update
            col = column to update
            align_params = an AlignmentParameters object with the input
        '''
        dx = align_params.dx
        ex = align_params.ex
        # Get all scores and max score used to assign score and this entry
        all_scores = [self.m_matrix.get_score(row, col-1) - dx, \
            self.iy_matrix.get_score(row, col-1) - ex]
        max_score = max(all_scores)

        # Condition to prevent negative scores in score matrix for local alignment
        if max_score < 0.0 and not align_params.global_alignment:
        	max_score = 0.0
        self.iy_matrix.set_score(row, col, max_score)

        # Assign pointers for the updated entry by pairing max score
        # with matrix or matrices from which it came
        pointerObjects = []
        matrix_list = [self.m_matrix, self.iy_matrix]
        max_matrices_by_index = [i for i, j in enumerate(all_scores) if fuzzy_equals(j,max_score)]
        pointerObjects = [matrix_list[k].get_entry_object(row, col-1) for k in max_matrices_by_index]

        self.iy_matrix.set_pointers(row, col, pointerObjects)

#### ------- MAIN METHODS ------- ####

def find_traceback_start(align_object, align_params):
    """
    Finds the location to start the traceback.

    Inputs:
        align_object = an AlignmentAlgorithm object with populated score matrices
        align_params = an AlignmentParameters object with the input

    Returns:
        Two variables:
        high_score = the high score
        high_score_objects = a list of ScoreMatrixEntry objects that are the traceback 
        starts (highest score in the proper location for the given alignment type).
    """
    # Identify alignment type
    global_alignment = align_params.global_alignment

    high_score_global = -1000.0  # arbitrary low sub-zero number (will be reset to 0.0 by boundary scores)
    high_score_local = 0.0  # score matrix will never have score < 0.0 for local
    high_score = 0.0

    score_matrices = [align_object.m_matrix, align_object.ix_matrix, align_object.iy_matrix]
    high_score_objects = []

    # Iterate through last column and last row of all matrices to get max score and position(s)
    if global_alignment:
        for matrix in score_matrices:
            for j in range(matrix.ncol-1):  # Avoid counting bottom right corner twice
                score = matrix.get_score(matrix.nrow-1, j)
                if score > high_score_global:
                    high_score_global = score
                    high_score_objects = [ matrix.get_entry_object(matrix.nrow-1, j) ]
                elif fuzzy_equals(score, high_score_global):
                    high_score_objects.append( matrix.get_entry_object(matrix.nrow-1, j) )
            for k in range(matrix.nrow):
                score = matrix.get_score(k, matrix.ncol-1)
                if score > high_score_global:
                    high_score_global = score
                    high_score_objects = [ matrix.get_entry_object(k, matrix.ncol-1) ]
                elif fuzzy_equals(score, high_score_global):
                    high_score_objects.append( matrix.get_entry_object(k, matrix.ncol-1) )
        high_score = high_score_global

    # Iterate through all of all matrices to get max score and position(s)
    elif not global_alignment:
        for matrix in score_matrices:
            for i in range(matrix.nrow):
                for m in range(matrix.ncol):
                    score = matrix.get_score(i, m)
                    if score > high_score_local:
                        high_score_local = score
                        high_score_objects = [ matrix.get_entry_object(i,m) ]
                    elif fuzzy_equals(score, high_score_local):
                        high_score_objects.append( matrix.get_entry_object(i,m) )
        high_score = high_score_local

    return "%.1f" % high_score, high_score_objects


def traceback(align_object, align_params, traceback_starts):
    """
    Generates a list of entries for the traceback and calls trace() function
    to perform the traceback recursively to generate the final traceback list.

    Input:
        align_object = an AlignmentAlgorithm object with populated score matrices
        align_params = an AlignmentParameters object with the input
        traceback_starts = a list of ScoreMatrixEntry objects that are the start points
    Returns:
        Two variables:
        master_list_unique = a list of final traceback lists
        print_string = a set of all tracebacks for printing which includes separate 
        tracebacks with ties for M, Ix, and Iy in the final entry at the end.
    """
    
    master_list = []
    for i in traceback_starts:
        old_list = [ (i.name, i.row, i.col) ]   # Initialize list with traceback start
        trace(i, master_list, old_list, align_params) # call trace() function recursively

    print_string = 'For tracebacks with ties for the final pointers (point to M, Ix, and Iy), only one alignment is ultimately printed.\n'
    for i in master_list:
        print_string = print_string + i[0][0]+ '(' + str(i[0][1]) + ',' + str(i[0][2]) +')'
        for j in i[1:]:
            print_string = print_string + '->' + j[0]+ '(' + str(j[1]) + ',' + str(j[2]) +')'
        print_string += '\n'

    # remove the last entry in the traceback in order to get rid of ties
    # for M, Ix, and Iy and find unique tracebacks
    
    for i in master_list:
        del i[-1]
    master_list_set = set([tuple(j) for j in master_list])
    master_list_unique = list(master_list_set)

    return master_list_unique, print_string

def trace(entry_object, master_list, old_list, align_params):
    '''
    Performs a traceback.
    Inputs:
        entry_object = an ScoreMatrixEntry object
        master_list = list of traceback entries
        old_list = an initialized list containing the traceback start point
        align_params = an AlignmentParameters object with the input
    Returns:
        Updates a final traceback list of lists with tuples of (name, row, column) 
        info for each object in a traceback.
    '''
    # terminate recursion if reach object score of 0.0 for local
    if fuzzy_equals(entry_object.score, 0.0) and not align_params.global_alignment:  # local alignment; entry object is the current object
        master_list.append(old_list)
        return
    # terminate recursion if reach first row or column for global
    elif align_params.global_alignment and (entry_object.row == 0 or entry_object.col == 0):   # global alignment
        master_list.append(old_list)
        return
    for i in entry_object.pointedObjects:
        info = (i.name, i.row, i.col)   # save info for each pointed-to object
        new_list = old_list[:]
        new_list.append(info)
        trace(i, master_list, new_list, align_params)


def write_alignment(output_file, align_params, traceback_list, high_score):
    '''
    Writes optimal alignment score and all alignments to an output file.
    Inputs:
        output_file = output file name
        align_params = an AlignmentParameters object with the input
        traceback_list = final list traceback list of tuples with (name, row, column)
        high_score = optimal alignment score
    Returns:
        Writes to output file formatted as specified in Project 1 guidelines
    '''
    seq_a = align_params.seq_a
    seq_b = align_params.seq_b

    alignment_pairs = []

    for i in traceback_list: # i is a single traceback
        traceback_string_A = ''
        traceback_string_B = ''
        if i == []: # In case the alignment fails
            alignment_pairs.append( (traceback_string_A, traceback_string_B) )
            break
        # convert the traceback to letters
        for entry in i: # entry is a tuple
            if entry[0] == 'M':
                traceback_string_A = seq_a[entry[1]-1] + traceback_string_A
                traceback_string_B = seq_b[entry[2]-1] + traceback_string_B
            elif entry[0] == 'Ix':
                traceback_string_A = seq_a[entry[1]-1] + traceback_string_A
                traceback_string_B = '_' + traceback_string_B
            else: # start is in Iy
                traceback_string_A = '_' + traceback_string_A
                traceback_string_B = seq_b[entry[2]-1] + traceback_string_B
        alignment_pairs.append( (traceback_string_A, traceback_string_B) )
        
    # When the alignment score is very low and the optimal alignments don't start in M, 
    # remove alignments that start or end with gaps.
    new_alignment_pairs = []
    for i in alignment_pairs:
        if not i[0] == '':
            if i[0][0] == '_' or i[1][0] == '_' or i[0][-1] == '_' or i[1][-1] == '_':
                continue
            else:
                new_alignment_pairs.append(i)
    alignment_pairs = new_alignment_pairs   

    # write to the file
    f=open(output_file, 'w')
    f.write(str(high_score)+'\n')
    for i in alignment_pairs:
        f.write('\n'+i[0]+'\n')
        f.write(i[1]+'\n')

    f.close()

def main():

    # check that the file is being properly used
    if (len(sys.argv) !=3):
        print("Please specify an input file and an output file as args.")
        return
        
    # input variables
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # read in the alignment data    
    align_params = AlignmentParameters()
    align_params.load_params_from_file(input_file)

    # populate the score matrices based on the input parameters
    align_object = AlignmentAlgorithm()
    align_object.populate_score_matrices(align_params)

    # perform a traceback and write the output to an output file; also
    # provide variables that can be printed (high_score, traceback_list, traceback_string)
    high_score, traceback_starts = find_traceback_start(align_object, align_params)
    traceback_list, traceback_string = traceback(align_object, align_params, traceback_starts)
    write_alignment(output_file, align_params, traceback_list, high_score)


if __name__=="__main__":
    main()
