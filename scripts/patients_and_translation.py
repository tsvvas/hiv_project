import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import Align
import matplotlib.pyplot as plt
import itertools

# HIV regions and patients lists
patients = ["p{}".format(i) for i in range(1, 12)]
hiv_regions = ["V3", "PR", "psi", "vpr", "vpu", "p1", "p2", "p6", "p7", "p15", "p17", "RRE"]


def get_patients():
    return patients


def get_regions():
    return hiv_regions


def preparing_data(haplo_seq_dict, days):
    """
    Function to help to construct tree
    Args:
        haplo_seq_dict:
        days:

    Returns:
        two lists
    """
    # making re pattern
    patt = r'_[\d]*_'

    # making dict in which we will save all our {name of seq: sequence} pairs
    seq_dict = {}

    # making dict in which we will save all our {day: [names of all seqs for this day]} pairs
    seq_name_days_dict = {day: [] for day in days}

    # fulling created dicts
    # note: type(obj) is SeqIO class
    for obj in haplo_seq_dict:

        # adding name:seq into seq_dict
        seq_dict[obj['desc']] = obj['seq']

        # using re to find out particular day for exact sequence
        res = re.search(patt, obj['desc'])

        # if day was found adding sequence to seq_name_days_dict
        if res != None:
            # note grouping here in order to make days look more understandable
            seq_name_days_dict[re.search(patt, obj['desc']).group(0).replace('_', '')].append(obj['desc'])
        # if can not find out what day, it means seq is reference
        else:
            pass

    return seq_dict, seq_name_days_dict


def translating_seqs(seq_dict):
    """
    Func to translate all seqs from seq_dict to proteins.

    Args:
        seq_dict: dict, dictionary of all sequences

    Returns:
        prot_dict: dict, dictionary of all {name of seq: protein} pairs
    """

    # initializing prot_dict
    prot_dict = {}

    # for-loop to translate all seqs from seq_dict
    for key in seq_dict:
        # making {name of seq: protein} pair in prot_dict
        prot_dict[key] = translate_dna(seq_dict[key], '-')

    return prot_dict



def translate_dna(seq, gap):
    """
    Translating sequence using codon_table (global) without any problems
    Args:
        seq: str, sequence to translate
        gap: str, gap that should be missed

    Returns:
        str, translated seq
    """

    # TODO check for stop codons
    # codon_table for translating
    codon_table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
        '---': '-',
    }

    # our protein
    prot = []
    res = ''

    # we will check every symbol
    for x in seq:

        # if for gap
        if x == gap:
            continue
        res += x

        # only when len(res) == 3 we will translate
        if len(res) == 3:
            if res in codon_table:
                prot.append(codon_table[res])
                res = ''

            # Be careful here
            else:
                prot.append('X')
                res = ''

    return ''.join(prot)
