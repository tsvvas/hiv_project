from patients_and_translation import translate_dna
from Bio import Align
import pandas as pd

# Dealing with hydrophobicity.csv

hydro_df = pd.read_csv('utils/hydrophobicity.csv', index_col = None)
hydro_fact = dict(zip(list(hydro_df['aa']),list(hydro_df['Kyte-Doolittle'])))
hydro_fact['X'] = 0 # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11  added
hydro_fact['*'] = 0
hydro_fact['-'] = 0


def calculating_snp(seq1, seq2):
    """
    Function to calculate SNP in two sequences using their alignment score

    Args:
        seq1: string
        seq2: string

    Return:
        snp: int, num of SNPs
    """
    aligner = Align.PairwiseAligner()
    score = aligner.score(seq1, seq2)

    snp = len(seq1) - score

    '''n = len(seq1)
    m = len(seq2)
    if m > n:
        n, m = m, n
        seq1, seq2 = seq2, seq1
    seq2 += ' '*(n-m)

    snp = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            snp += 1'''

    return snp


def calculate_metric(seq):
    """
    Function to calculate hydrophobicity for DNA sequence
    Args:
        seq: string, sequnce to translate
    Return:
        res: float, metric result
    """
    prot = translate_dna(seq, gap='-')

    res = 0

    for aa in prot:
        res += hydro_fact[aa]

    return res