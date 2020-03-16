from Bio import Align
from Bio import Seq
from Bio.Alphabet import IUPAC, Gapped
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

    # creating aligner
    aligner = Align.PairwiseAligner()

    # calculating score
    score = aligner.score(seq1, seq2)

    # finding out num of SNPs
    snp = len(seq1) - score

    return snp


def calculate_metric(seq):
    """
    Function to calculate hydrophobicity for DNA sequence
    Args:
        seq: string, sequnce to translate
    Return:
        float, metric result
    """

    # translating protein
    prot = Seq(seq, Gapped(IUPAC.unambiguous_dna)).ungap().translate()

    # calculating metrics
    res = 0
    for aa in prot:

        # adding to res hydrophobicity of exact amino acid acording to hydro_fact
        res += hydro_fact[aa]

    return res