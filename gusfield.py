from msa_utils import (pairwise_distance_matrix, msa_from_guide_tree,
                       construct_msa_dict)
import numpy as np


def find_center_seq_idx(distance_matrix: list):
    ''' Finds the key of the sequence with the lowest sum of pairs score'''

    sp_score = np.nansum(distance_matrix, axis=0)
    center_string_idx = np.argmin(sp_score)

    return center_string_idx

def gusfield_guide_tree(distance_matrix:list, extension_order=False):
    ''' Constructs the star tree according to gusfields algortihm. 
        Extention order can be supplied to change the order of extentions'''
    if not type(extension_order) is np.ndarray:
        extension_order = list(range(len(distance_matrix[0])))
    else:
        extension_order = list(extension_order)

    center_string_idx = find_center_seq_idx(distance_matrix)
    extension_order.remove(center_string_idx)
    guide_tree = []
    for string_idx in extension_order:
        guide_tree.append((center_string_idx, string_idx))
    
    return guide_tree


def gusfield_msa(sequence_dict: dict, score_dict: dict, extension_order=False,
                 output_tree=False):
    ''' Performs MSA on a dictionary of sequences using star tree created using
        Gusfields algorithm'''
    (seq_names, seq_values) = zip(*sequence_dict.items())

    adjacency_matrix = pairwise_distance_matrix(seq_values, score_dict)
    star_tree = gusfield_guide_tree(adjacency_matrix, extension_order)
    msa, msa_seq_map = msa_from_guide_tree(star_tree, seq_values, score_dict)

    msa_dict = construct_msa_dict(msa, msa_seq_map, seq_names)
    
    return (msa_dict, star_tree) if output_tree else msa_dict


def gusfield_custom_extension(distance_matrix, center_string_idx, descending=False):
    ''' Generates a list used to alter the Gusfield extension order, to extend
        the MSA by either the closest or furthest sequences first, sorted
        by pairwise distance to the root'''
    best_extension = np.argsort(distance_matrix[center_string_idx])
    worst_extension = np.flip(best_extension)
    if descending:
        return worst_extension
    return best_extension