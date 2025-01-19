"""
Gene Rearrangement Simulation

Author: Yamila Timmer
Date: 19-01-2025
Version: 0.1

Description:
This script processes the DNA sequence of a somatic cell to simulate gene rearrangement and production of
immunoglobuline of the variable part of the heavy chain of a B-cell. Within this DNA sequence, different clusters
are identified and manipulated with simulated processes such as P and N nucleotide addition.

Key functionalities:
1. Reads a sequence from the file `sequence.txt`.
2. Identifies RSS patterns (as indexes) and uses these to isolate V-, D-, and J-clusters.
3. Simulates processes that add extra variation to the clusters, including: P-, and N nucleotide additions and
exonuclease trimming.
4. Makes merged clusters, each having a randomly chosen V-, D-, and J-cluster. And makes sure that the merged clusters
start with a start codon (ATG).
5. Translates the merged clusters into amino acid sequences.
6. Validates if the resulting sequence represents a functional heavy chain.


Usage:
- Ensure that a valid DNA sequence is present in a text file named `sequentie.txt`.
- Run the script, and it will keep attempting to generate sequences until 10 functional sequences are generated.

Run the script with:
    `python gene_rearrangement.py`

"""

import random


def read_sequence():
    """
    Reads the DNA sequence from the file sequence.txt.

    :return: dna_seq: (str) The DNA sequence.
    """

    with open("../Gene Rearrangement Simulation/sequence.txt", "r") as seq:
        dna_seq = ""
        for lines in seq:
            dna_seq += lines

        return dna_seq


def find_rss(dna_seq, rss):
    """
    Finds indexes of RSS (heptamer and nonamer)
    :param dna_seq: (str) The sequence to search.
    :param rss: (str) The sequence of the RSS (heptamer/nonamer) to search for.
    :return: rss_indexes: (list) The indexes of the RSS (heptamer/nonamer).
    """

    rss_indexes = []
    rss_index = ""
    pos = 0

    # Will keep going as long as .find does not return '-1', .find will return '-1' once no more rss sequences are found
    while rss_index != -1:
        rss_index = dna_seq.find(rss, pos)

        if rss_index!= -1:
            pos = rss_index + 1
            rss_indexes.append(rss_index)

    return rss_indexes

def find_v_cluster(dna_seq, heptamer_indexes, nonamer_indexes):
    """
        Finds sequences of v clusters.
        :param dna_seq: (str) The sequence to search.
        :param heptamer_indexes: (list) The indexes of all found heptameres.
        :param nonamer_indexes: (list) The indexes of all found nonameres.
        :return: cluster_sequences (list) The sequences of v clusters.
        """

    cluster_sequences = []

    # Adds the first V-cluster (from ATG to the first heptamer
    v_cluster = dna_seq[0:heptamer_indexes[0]]
    cluster_sequences.append(v_cluster)

    #
    for nonamer_index in nonamer_indexes:
        for heptamer_index in heptamer_indexes:
            if heptamer_index > nonamer_index + 9:
                v_cluster = dna_seq[nonamer_index + 9:heptamer_index]
                cluster_sequences.append(v_cluster)
                break

    return cluster_sequences


def find_d_cluster(dna_seq, heptamer_indexes):
    """
    Finds sequences of d clusters
    :param dna_seq: (str) The sequence to search.
    :param heptamer_indexes: (list) The indexes of the heptameres.
    :return: cluster_sequences (list) The sequences of D clusters.
    """

    cluster_sequences = []

    # Iterates over the list of heptamer_indexes, '-1' is needed because the final heptamer cannot be iterated with [i+1]
    for i in range(len(heptamer_indexes) - 1):

        # D clusters are found between the end of RSS (heptamer + 7) and the start of a new RSS (heptamer -1)
        d_cluster = dna_seq[heptamer_indexes[i] + 7:heptamer_indexes[i + 1] -1]
        cluster_sequences.append(d_cluster)

    return cluster_sequences


def find_j_cluster(dna_seq, heptamer_indexes, nonamer_indexes):
    """
    Finds sequences of j clusters
    :param dna_seq: (str) The sequence to search.
    :param heptamer_indexes: (list) The indexes of the heptameres.
    :param nonamer_indexes: (list) The indexes of the nonameres.
    :return:
        cluster_sequences (list) The sequences of J clusters.
        heavy_index (int) The index of the start of the heavy chain.
    """

    cluster_sequences = []
    count = 0

    # J-clusters are found between the end of a heptamer [heptamer + 1] and the beginning of a nonamer [nonamer - 1]
    for nonamer_index in nonamer_indexes:
        for heptamer_index in heptamer_indexes:
            if nonamer_index -1 > heptamer_index + 1:
                j_cluster = sequence[heptamer_index + 1:nonamer_index -1]
                cluster_sequences.append(j_cluster)
                count += 1

                break

    # Adding final J-cluster, which is in between the last heptamer and the start of the constant part of the heavy chain
    heavy_constant = "CCTCCACCAAGG"
    heavy_index = dna_seq.find(heavy_constant)

    j_cluster = dna_seq[heptamer_indexes[-1] + 1: heavy_index]
    cluster_sequences.append(j_cluster)

    return cluster_sequences, heavy_index


def exonuclease_trimming(cluster_sequences):
    """
    Trims a random amount of nucleotides on both sides of the cluster.

    :param cluster_sequences (list) The sequences of the clusters (D, V, J) to trim.
    :return: cluster_sequences (list) The sequences of the clusters (D, V, J) after exonuclease trimming.
    """

    for i in range(len(cluster_sequences)):
        remove_nucleotides_int = random.randint(1, 10)
        cluster_sequences[i] = cluster_sequences[i][remove_nucleotides_int:]

    for i in range(len(cluster_sequences)):
        remove_nucleotides_int = random.randint(1, 10)
        cluster_sequences[i] = cluster_sequences[i][:remove_nucleotides_int]

    return cluster_sequences

def n_addition_seq():
    """
    Generates a random sequence of nucleotides (1-10), to be added to the cluster.

    :return: seq (str) The random sequence of nucleotides.
    """

    nucleotides = ["A", "C", "T", "G"]
    addition_number = random.randint(1, 10)
    seq = ""

    for _ in range(addition_number):
        n_nucleotides = random.choice(nucleotides)
        seq += n_nucleotides

    return seq


def n_addition(cluster_sequences):
    """
    Adds randomly generated nucleotides to the cluster.

    :param cluster_sequences (list) The sequences of the clusters (D, V, J) after exonuclease trimming.
    :return: cluster_sequences (list) The sequences of the clusters (D, V, J) after n_addition.
    """

    for i in range(len(cluster_sequences)):
        n_seq = n_addition_seq()  # Generate random n-nucleotides
        cluster_sequences[i] = n_seq + cluster_sequences[i] # Add random sequence in front of cluster


    for i in range(len(cluster_sequences)):
        n_seq = n_addition_seq()
        cluster_sequences[i] += n_seq # Add random sequence behind cluster

    return cluster_sequences


def p_addition(cluster_sequences):
    """

    :param cluster_sequences: (list) The sequences of the clusters (D, V, J) after n_addition.
    :return: cluster_sequences: (list) The sequences of the clusters (D, V, J) after p_addition.
    """
    nucleotide_dict = {"A": "T", "C": "G", "T": "A", "G": "C"}

    for i in range(len(cluster_sequences)):

        # Make "palindrome" and invert it, so that it can properly be appended
        palindrome_nucleotides = (nucleotide_dict[cluster_sequences[i][-2]] + nucleotide_dict[cluster_sequences[i][-1]])[::-1]
        cluster_sequences[i] += palindrome_nucleotides

    for i in range(len(cluster_sequences)):
        # Make "palindrome" and invert it, so that it can properly be appended
        palindrome_nucleotides = (nucleotide_dict[cluster_sequences[i][-2]] + nucleotide_dict[cluster_sequences[i][-1]])[::-1]
        cluster_sequences[i] = palindrome_nucleotides + cluster_sequences[i]

    return cluster_sequences


def merge_clusters(j_cluster, d_cluster, v_cluster, dna_seq, heavy_index):
    """
    Merges the clusters together, randomly choosing a V-, D-, and J-cluster.

    :param j_cluster: (list) The sequences of the found J-clusters, after trimming and p- and n-addition.
    :param d_cluster: (list) The sequences of the found D-clusters, after trimming and p- and n-addition.
    :param v_cluster: (list) The sequences of the found V-clusters, after trimming and p- and n-addition.
    :param dna_seq: (str) The DNA sequence.
    :param heavy_index (int) The index of the start of the heavy chain.
    :return: sequence_merged: (str) The merged sequence of a random D-, V- and J-cluster.
    """

    sequence_merged = ""
    sequence_merged += (random.choice(d_cluster) + random.choice(j_cluster) + random.choice(v_cluster) +
                        dna_seq[heavy_index:])

    return sequence_merged


def start_sequence(sequence_merged):
    """
    Cuts sequence_new to start with ATG (starting codon).

    :param sequence_merged: (str) The merged sequence of a random D-, V- and J-cluster.
    :return: sequence_start: (str) cropped merged sequence, starting from ATG (start codon).
    """

    start_index = sequence_merged.find("ATG")
    sequence_start = sequence_merged[start_index:]

    return sequence_start


def convert_to_aminoacid(sequence_start):
    """
    Converts the DNA sequence to an amino acid sequence.

    :param sequence_start: (str) cropped merged sequence, starting from ATG (start codon).
    :return: amino_seq: (str) The amino acid sequence.
    """
    dna_codons = {
        # 'M' - START, '_' - STOP
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A", "TGT": "C", "TGC": "C", "GAT": "D", "GAC": "D", "GAA": "E",
        "GAG": "E", "TTT": "F", "TTC": "F", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G", "CAT": "H", "CAC": "H",
        "ATA": "I", "ATT": "I", "ATC": "I", "AAA": "K", "AAG": "K", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
        "CTA": "L", "CTG": "L", "ATG": "M", "AAT": "N", "AAC": "N", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R", "TCT": "S",
        "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TGG": "W", "TAT": "Y", "TAC": "Y", "TAA": "_", "TAG": "_",
        "TGA": "_"
    }

    start = True
    amino_seq = ""

    # Convert codons to amino acids
    for i in range(0, len(sequence_start), 3):
            codon = sequence_start[i:i + 3]

            # Stops when the iteration has reached the end of the sequence and there are not enough nucleotides to
            # form one complete codon (len(codon) < 3)
            if len(codon) != 3:
                continue

            amino_acid = dna_codons[codon]

            if amino_acid == 'M':
                start = True
            elif amino_acid == "_":
                start = False

            # Will keep translating to amino acids, until stop codon is reached
            if start:
                amino_seq += amino_acid

    return amino_seq



def functional_check(sequence_start, amino_seq):
    """
    Checks whether the amino acid sequence is functional, based on the constant part of the heavy chain
    (which is lacking a stop codon).

    :param sequence_start: (str) cropped merged sequence, starting from ATG (start codon).
    :param amino_seq: (str) The amino acid sequence.
    :return:
        functional_sequences: (dict) containing the functional amino acid sequences and the corresponding DNA sequence.
        non_functional_sequences: (dict) containing the non-functional amino acid sequences and the corresponding DNA sequence.

    """

    functional_heavy_chain = ("STKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNV"
                              "NHKPSNTKVDKKVGERPAQGGRVSAGSQAQRSCLDASRLCSPSPGQQGRPRLPLHPEASARPTHAQGEG")
    functional_sequences = {}
    non_functional_sequences = {}

    # Adds amino sequence and corresponding dna sequence to functional_sequences dict if constant heavy chain sequence is found
    if functional_heavy_chain in amino_seq:
        functional_sequences[amino_seq] = sequence_start

    # Adds amino sequence and corresponding dna sequence to non_functional_sequences dict when constant heavy chain
    # sequence is not found
    else:
        non_functional_sequences[amino_seq] = sequence_start

    return functional_sequences, non_functional_sequences


def print_result(func_seq, non_func_seq):
    """
    Neatly prints result of the pipeline, printing both the amino and dna sequence of the functional and non-functional
    sequences.

    :param func_seq: (dict) containing the functional amino acid sequences and the corresponding DNA sequence.
    :param non_func_seq: (dict) containing the non-functional amino acid sequences and the corresponding DNA sequence.
    :return:
    """

    count = 1

    print("Functional sequences:")
    for amino_seq, sequence_start in func_seq.items():

        print(count, "Amino Sequence:", amino_seq, "DNA Sequence:", sequence_start)
        count +=1

    count = 1
    print("\nNon-functional sequences:")
    for amino_seq, sequence_start in non_func_seq.items():
        print(count, "Amino Sequence:", amino_seq, "DNA Sequence:", sequence_start)
        count += 1

    print("\nFunctional-sequences found:", len(func_seq), "\nNon-functional-sequences found:", len(non_func_seq))


if __name__ == "__main__":
    # Calling all pipeline functions in the correct order, with the correct arguments
    sequence = read_sequence()

    heptamer_indexes = find_rss(sequence, rss = "CACAGTG")
    nonamer_indexes = find_rss(sequence, rss = "ACAAAAACC")

    v_cluster_sequences = find_v_cluster(sequence, heptamer_indexes, nonamer_indexes)
    d_cluster_sequences = find_d_cluster(sequence, heptamer_indexes)
    j_cluster_sequences, heavy_index = find_j_cluster(sequence, heptamer_indexes, nonamer_indexes)

    functional_sequences = {}
    non_functional_sequences = {}

    max_functional_sequences = 10

    # Pipeline keeps going until 10 functional sequences have been found
    while len(functional_sequences) < max_functional_sequences:
        n_seq = n_addition_seq()

        trimmed_j_clusters = exonuclease_trimming(j_cluster_sequences)
        trimmed_d_clusters = exonuclease_trimming(d_cluster_sequences)
        trimmed_v_clusters = exonuclease_trimming(v_cluster_sequences)

        n_added_j_clusters = n_addition(trimmed_j_clusters)
        n_added_d_clusters = n_addition(trimmed_d_clusters)
        n_added_v_clusters = n_addition(trimmed_v_clusters)

        p_added_j_clusters = p_addition(j_cluster_sequences)
        p_added_d_clusters = p_addition(d_cluster_sequences)
        p_added_v_clusters = p_addition(v_cluster_sequences)

        sequence_merged = merge_clusters(p_added_j_clusters, p_added_d_clusters, p_added_v_clusters, sequence,
                                         heavy_index)
        sequence_start = start_sequence(sequence_merged)

        amino_seq = convert_to_aminoacid(sequence_start)

        func_seq, non_func_seq = functional_check(sequence_start, amino_seq)

        functional_sequences.update(func_seq)
        non_functional_sequences.update(non_func_seq)

    print_result(functional_sequences, non_functional_sequences)
