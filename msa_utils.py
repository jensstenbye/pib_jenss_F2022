from itertools import combinations
from Bio import SeqIO
import numpy as np
################################################################################
# Pairwise alignment functions
################################################################################
def pairalign_score_matrix(seq_1: str, seq_2: str, score_dict: dict):
    ''' Returns scoring matrix of a pairwise global alignment of two sequences
    '''

    # Create score matrix and fill zero'th column and row with gap scores
    score_matrix = np.zeros((len(seq_1) + 1, len(seq_2) + 1))
    score_matrix[:, 0] = np.arange(0, len(seq_1) + 1) * score_dict['gap']
    score_matrix[0, :] = np.arange(0, len(seq_2) + 1) * score_dict['gap']

    # Fill in rest of score matrix
    for i, base_1 in enumerate(seq_1, 1):
        for j, base_2 in enumerate(seq_2, 1):
            from_left = score_matrix[i, (j-1)] + score_dict['gap']
            from_above = score_matrix[(i-1), j] + score_dict['gap']
            from_across = (score_matrix[(i-1), (j-1)]
                           + score_dict[base_1][base_2])

            score_matrix[i, j] = min(from_left, from_above, from_across)

    return score_matrix


def pairalign_backtracking(seq_1: str, seq_2: str, score_dict: dict):
    ''' Returns aligned sequences from pairwise alignment'''
    score_matrix = pairalign_score_matrix(seq_1, seq_2, score_dict)

    # Start backtracking at bottom right in the scoring matrix
    i, j = score_matrix.shape
    i, j = i - 1, j - 1

    # Build aligned sequences from behind
    alignment_1 = []
    alignment_2 = []
    while (i != 0) and (j != 0):
        if (score_matrix[i, j] == (score_matrix[i-1, j-1]
                                   + score_dict[seq_1[i-1]][seq_2[j-1]])):
            # Coming from across
            i -= 1
            j -= 1
            alignment_1.append(seq_1[i])
            alignment_2.append(seq_2[j])

        elif (score_matrix[i, j] == score_matrix[i-1, j] + score_dict['gap']):
            # Coming from right
            i -= 1
            alignment_1.append(seq_1[i])
            alignment_2.append('-')

        elif (score_matrix[i, j] == score_matrix[i, j-1] + score_dict['gap']):
            # Coming from above
            j -= 1
            alignment_1.append('-')
            alignment_2.append(seq_2[j])

        else:
            raise Exception(f'No backtrack path, exciting at i:{i}, j:{j}')

    # If edge reached, fill with gaps
    while i > 0:
        i -= 1
        alignment_1.append(seq_1[i])
        alignment_2.append('-')
    while j > 0:
        j -= 1
        alignment_1.append('-')
        alignment_2.append(seq_2[j])

    alignment_1 = ''.join(alignment_1[::-1])
    alignment_2 = ''.join(alignment_2[::-1])
    return alignment_1, alignment_2


def pairwise_distance_matrix(sequence_list: list, score_dict: dict):
    ''' Creates matrix of pairwise alignment scores between sequences in list
    '''
    pairwise_score_matrix = np.empty((len(sequence_list), len(sequence_list)))
    pairwise_score_matrix.fill(np.nan)

    pairwise_iterator = combinations(enumerate(sequence_list), 2)
    for (idx_seq1, seq_1), (idx_seq2, seq_2) in pairwise_iterator:
        pairalign_score = pairalign_score_matrix(seq_1, seq_2, score_dict)[-1, -1]

        pairwise_score_matrix[idx_seq1, idx_seq2] = pairalign_score
        pairwise_score_matrix[idx_seq2, idx_seq1] = pairalign_score

    return pairwise_score_matrix


################################################################################
# MSA construction functions
################################################################################
def extend_multiple_alignment(multiple_alignment, msa_seq_idx,
                              align_msa_seq, align_new_seq):
    ''' Extends a multiple alignment with an alignment between
        a string in msa and a new string to add
    '''
    # Track position in string of current multiple alignment
    # and of the same string in the new pairwise alignment
    pointer_multiple = 0
    pointer_pairwise = 0
    extended_sequence = []

    # Extend multiple alignment handling substitution, 
    # insertion and deletion cases until at the end of either string
    while ((pointer_pairwise <= (len(align_msa_seq)-1)) and
           (pointer_multiple <= (len(multiple_alignment[msa_seq_idx])-1))):

        if align_msa_seq[pointer_pairwise] == '-':
            # Insertion (i.e. gap in msa-sequence)
            extended_sequence.append(align_new_seq[pointer_pairwise])

            # If gap not already present in msa, add it to all sequences
            if multiple_alignment[msa_seq_idx][pointer_multiple] != '-':
                for sequence in multiple_alignment:
                    sequence.insert(pointer_multiple, '-')

            pointer_multiple += 1
            pointer_pairwise += 1

        elif multiple_alignment[msa_seq_idx][pointer_multiple] != '-':
            # Substitution or deletion
            extended_sequence.append(align_new_seq[pointer_pairwise])
            pointer_multiple += 1
            pointer_pairwise += 1

        else:
            # Gaps in msa not in alignment (another insertion case)
            extended_sequence.append('-')
            pointer_multiple += 1

    while pointer_multiple < (len(multiple_alignment[msa_seq_idx])):
        # Extend trailing gaps in new sequence
        extended_sequence.append('-')
        pointer_multiple += 1

    while (pointer_pairwise <= (len(align_msa_seq)-1)):
        # Extend trailing gaps in msa and extend sequence with remaining bases
        extended_sequence.append(align_new_seq[pointer_pairwise])
        for sequence in multiple_alignment:
            sequence.append('-')   
        pointer_pairwise += 1     

    multiple_alignment.append(extended_sequence)
    return


def msa_from_guide_tree(guide_tree: list, seq_values: list, score_dict: dict):
    ''' Constructs a multiple sequence alignment from a guide tree'''
    msa = []

    # msa_seq_map maps from index in seq_values to the corresponding sequence
    # in msa used to pair the sequences with their names
    msa_seq_map = [None] * len(seq_values)

    # Add the root and first node to msa
    root_seq, first_node = pairalign_backtracking(seq_values[guide_tree[0][0]],
                                                  seq_values[guide_tree[0][1]],
                                                  score_dict)
    msa.extend((list(root_seq), list(first_node)))
    msa_seq_map[guide_tree[0][0]], msa_seq_map[guide_tree[0][1]] = 0, 1

    # Iteratively add the remaining nodes
    for ext_idx in range(1, len(guide_tree)):
        align_tree, align_node = pairalign_backtracking(seq_values[guide_tree[ext_idx][0]],
                                                        seq_values[guide_tree[ext_idx][1]],
                                                        score_dict)        
        extend_multiple_alignment(msa, msa_seq_map[guide_tree[ext_idx][0]],
                                  align_tree, align_node)
        msa_seq_map[guide_tree[ext_idx][1]] = len(msa) - 1         

    return msa, msa_seq_map


################################################################################
# Score functions
################################################################################
def get_pairalign_score(sequence_1, sequence_2, score_dict):
    ''' Calcualte the score of a pairwise alignment'''
    score = 0

    for base_1, base_2 in zip(sequence_1, sequence_2):
        if (base_1 != '-') and (base_2 != '-'):
            # Base in both sequences
            score += score_dict[base_1][base_2]
        elif (base_1 == '-') and (base_2 == '-'):
            # Gaps in both sequences
            continue
        else:
            # Gap in one sequence
            score += score_dict['gap']  
    return score      


def get_sum_of_pairs_score(multiple_alignment_dict: dict, score_dict: dict):
    ''' Calculates the sum of pairs score for a multiple alignment, by summing
        the scores of all induced alignments'''
    score = 0
    pairwise_iterator = combinations(multiple_alignment_dict.values(), 2)
    for (seq_1, seq_2) in pairwise_iterator:
        score += get_pairalign_score(seq_1, seq_2, score_dict)

    return score


################################################################################
# datastructure functions
################################################################################
def parse_score_matrix(path):
    ''' Parses file containing score matrix into a score dictionary'''

    score_dict = {}
    with open(path, 'r') as score_file:
        # Parse gap score and build datastructure
        base_list = next(score_file).split()
        score_dict['gap'] = int(base_list.pop(0))
        for base in base_list:
            score_dict[base] = {}

        # Fill datastructure with scores
        for line in score_file:
            score_list = line.split()
            base = score_list.pop(0)
            for matched_base, score in zip(base_list, score_list):
                score_dict[base][matched_base] = int(score)

    return score_dict

def fasta_to_dict(path):
    ''' Parses fasta file into dict with name/sequence as key/value pair'''
    sequence_dict = {}
    for record in SeqIO.parse(path, "fasta"):
        sequence_dict[record.name] = str(record.seq)

    return sequence_dict

def dict_to_fasta(path, sequence_dict):
    ''' Writes sequence_dict to a fasta output file'''
    with open(path, 'w') as outfile:
        for name in sequence_dict:
            outfile.write(f'>{name}\n')
            outfile.write(f'{sequence_dict[name]}\n\n')
    return


def construct_msa_dict(msa: list, msa_seq_map: list, seq_names:list):
    ''' Constructs a dictionary with the aligned sequences from msa matrix,
        sequence names and a mapping between the two'''
    msa = [''.join(sequence) for sequence in msa]

    msa_dict = {}
    for name_idx, msa_idx in enumerate(msa_seq_map):
        msa_dict[seq_names[name_idx]] = msa[msa_idx]
    return msa_dict
