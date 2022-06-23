from itertools import combinations
import random
import numpy as np

def shuffle_dict(dictionary):
    ''' Shuffles a dictionary'''
    l = list(dictionary.items())
    random.shuffle(l)
    return dict(l)

def simulate_sequences(length, number, distance, path):
    ''' Simulates set of sequences originating from the same root by some
        distance d.'''

    base_mapping = {0:"A", 1:"T", 2:"G", 3:"C"}
    mutation_dict = {"A":"TGC", "T":"AGC", "G":"ATC", "C":"ATG"}
    sequence_dict = {}

    # Generate root sequence
    root = ''
    for i in range(length):
        root += base_mapping[random.randint(0, 3)]
        sequence_dict['>seq0'] = root

    # Generate leaf sequences
    for seq_idx in range(1, number):
        leaf = ""
        for i in range(length):
            # Generate mutation based on phylogenetic distance
            if random.randint(0, 100) <= int(distance*100):
                leaf += mutation_dict[root[i]][random.randint(0, 2)]
            else:
                leaf += root[i]
        sequence_dict[f'>seq{seq_idx}'] = leaf
    # Shuffle dict so closely related sequence dosnt necessarily follow eachother
    sequence_dict = shuffle_dict(sequence_dict)
    with open(path, 'w') as outfile:
        for name, sequence in sequence_dict.items():
            outfile.write(f'{name}\n')
            outfile.write(f'{sequence}\n\n')
            
    return sequence_dict

def simulate_clusters(length, number_per_cluster, clusters, within_distance, 
                      between_distance, path):
    ''' Simulates set of clusters of sequences closely related within, distantly
        related between'''

    base_mapping = {0:"A", 1:"T", 2:"G", 3:"C"}
    mutation_dict = {"A":"TGC", "T":"AGC", "G":"ATC", "C":"ATG"}
    root_dict = {}

    # Generate root sequence
    root = ''
    for i in range(length):
        root += base_mapping[random.randint(0, 3)]
        root_dict['>seq0'] = root

    # Generate cluster roots
    for seq_idx in range(1, clusters):
        cluster_root = ""
        for i in range(length):
            # Generate mutation based on phylogenetic distance
            if random.randint(0, 100) <= int(between_distance*100):
                cluster_root += mutation_dict[root[i]][random.randint(0, 2)]
            else:
                cluster_root += root[i]
        root_dict[f'>seq{seq_idx}'] = cluster_root
    
    # Generate sequences within each cluster
    cluster_dict = {}
    for name, root in root_dict.items():
        cluster_dict[name] = root
        for seq_idx in range(1, number_per_cluster):
            cluster_leaf = ""
            for i in range(length):
                # Generate mutation based on phylogenetic distance
                if random.randint(0, 100) <= int(within_distance*100):
                    cluster_leaf += mutation_dict[root[i]][random.randint(0, 2)]
                else:
                    cluster_leaf += root[i]   
            cluster_dict[f'{name}.{seq_idx}'] = cluster_leaf    
    # Shuffle dict so closely related sequence dosnt necessarily follow eachother
    cluster_dict = shuffle_dict(cluster_dict)
    with open(path, 'w') as outfile:
        for name, sequence in cluster_dict.items():
            outfile.write(f'{name}\n')
            outfile.write(f'{sequence}\n\n')
    return cluster_dict


def simulate_adj_matrix(k, dist_range=20):
    ''' Simulates a adjacency matrix given a number of sequences. dist_range,
        determines the maximum distance between two seqs'''
    
    sim_pairwise_score_matrix = np.empty((k, k))
    sim_pairwise_score_matrix.fill(np.nan)

    pairwise_iterator = combinations(range(k), 2)
    for (idx_seq1, idx_seq2) in pairwise_iterator:
        
        sim_pairalign_score = random.randint(0, 20)

        sim_pairwise_score_matrix[idx_seq1, idx_seq2] = sim_pairalign_score
        sim_pairwise_score_matrix[idx_seq2, idx_seq1] = sim_pairalign_score

    return sim_pairwise_score_matrix
