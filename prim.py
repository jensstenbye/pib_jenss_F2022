import numpy as np
from msa_utils import (pairwise_distance_matrix, msa_from_guide_tree,
                       construct_msa_dict)


def prim_guide_tree(distance_matrix, root_idx=0):
    ''' Constructs the MST using a matrix containing the pairwise distance
        between sequences. Returns a list containing the order of tree 
        extensions'''
    guide_tree = []
    k = len(distance_matrix[0])
    unlinked_node_idx = list(range(k))

    # Prims A holds the minimum distance between nodes in the MST and unadded
    # nodes. Together with their respective indices. First the root is added
    prims_A = [(root_idx, np.inf)] * k
    unlinked_node_idx.remove(root_idx)
    for idx in unlinked_node_idx:
        prims_A[idx] = (root_idx, distance_matrix[root_idx, idx])

    while unlinked_node_idx:  #(O(k))
        # Find the node closest to the tree (O(k))
        node_idx, (tree_idx, score) = min(enumerate(prims_A),
                                          key=lambda x: x[1][1])
        guide_tree.append((tree_idx, node_idx))
        
        # Update prims A (O(k))
        prims_A[node_idx] = (node_idx, np.inf)

        unlinked_node_idx.remove(node_idx)
        for idx in unlinked_node_idx:
            if distance_matrix[node_idx, idx] < prims_A[idx][1]:
                prims_A[idx] = ((node_idx, distance_matrix[node_idx, idx]))

    return guide_tree


def prim_msa(sequence_dict: dict, score_dict: dict, root_idx=0, 
             output_tree=False):
    ''' Performs MSA on a dictionary of sequences using MST created using
        Prims algorithm'''
    (seq_names, seq_values) = zip(*sequence_dict.items())

    adjacency_matrix = pairwise_distance_matrix(seq_values, score_dict)
    mst = prim_guide_tree(adjacency_matrix, root_idx)
    msa, msa_seq_map = msa_from_guide_tree(mst, seq_values, score_dict)

    msa_dict = construct_msa_dict(msa, msa_seq_map, seq_names)

    return (msa_dict, mst) if output_tree else msa_dict


def prim_custom_root(distance_matrix, biggest_pair_dist=False):
    ''' Finds the idx of the root with the largest or shortest sum of 
        pairwise distances to supply to Prim guide tree'''
    sp_score = np.nansum(distance_matrix, axis=0)

    worst_seq_idx = np.argmax(sp_score)
    best_seq_idx = np.argmin(sp_score)
    if biggest_pair_dist:
        return worst_seq_idx
    return best_seq_idx
