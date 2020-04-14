import pandas
import os

path, _ = os.path.split(__package__)
kidera_factors = pandas.read_csv(os.path.join(path, 'data', 'kidera.csv'), header=None, index_col=0)
aa_properties = pandas.read_csv(os.path.join(path, 'data', 'aa_property_table.csv'), header=0, index_col=0)
aa_hydrophobicity = pandas.read_csv(os.path.join(path, 'data', 'hydrophobicity.csv'))

symbol_lookup = { 'ALA': 'A', 'ARG': 'R',
                  'ASN': 'N', 'ASP': 'D',
                  'CYS': 'C', 'GLN': 'Q',
                  'GLU': 'E', 'GLY': 'G',
                  'HIS': 'H', 'ILE': 'I',
                  'LEU': 'L', 'LYS': 'K',
                  'MET': 'M', 'PHE': 'F',
                  'PRO': 'P', 'SER': 'S',
                  'THR': 'T', 'TRP': 'W',
                  'TYR': 'Y', 'VAL': 'V' }

kidera_factors.index = kidera_factors.index.map(lambda x: symbol_lookup[x])
aa_properties.index = aa_properties.index.map(lambda x: symbol_lookup[x])
aa_hydrophobicity = aa_hydrophobicity.set_index('aa', drop=True)


def score_sequence(sequence, norm=False):
    """
    Scores Kidera factors for given sequence

    Args:
        sequence: your sequence
        norm: if False returns sum, else mean of Kidera factors

    Returns:
        pd.DataFrame: score
    """
    if '*' in sequence:
        sequence._data = sequence._data.replace('*', '')
    if norm:
        return kidera_factors.loc[list(sequence)].sum() / len(sequence)
    else:
        return kidera_factors.loc[list(sequence)].sum()


def aaprop_sequence(sequence, norm=False):
    """
    Scores most popular amino acid properties for given sequence including
    Kidera factors

    Args:
        sequence: your sequence
        norm: if False returns sum, else mean of Kidera factors

    Returns:
        pd.DataFrame: score
    """
    if '*' in sequence:
        sequence._data = sequence._data.replace('*', '')
    if norm:
        return aa_properties.loc[list(sequence), :].sum() / len(sequence)
    else:
        return aa_properties.loc[list(sequence), :].sum()


def score_hydrophobicity_sequence(sequence, scale="Kyte-Doolittle", norm=False):
    """
    Scores hydrophobicity of amino acids with different scales from
    http://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/650/Hydrophobicity_scales.html

    Args:
        sequence: your sequence
        scale: One of the following scales:
            Kyte-Doolittle  Hopp-Woods  Cornette  Eisenberg   Rose  Janin   Engelman GES
        norm: if False returns sum, else mean of the hydrophobicity scores

    Returns:
        pd.DataFrame: score
    """
    if '*' in sequence:
        sequence._data = sequence._data.replace('*', '')
    if norm:
        return aa_hydrophobicity.loc[list(sequence), scale].sum() / len(sequence)
    else:
        return aa_hydrophobicity.loc[list(sequence), scale].sum()


def build_epitopes(seq, length):
    """
    Function returns possible epitopes of length. If remaining length of the seq is to short it will form one last
    epitope

    Args:
        seq (tuple): translated seq
        length (int): maximum length of one epitope

    Returns:
        list: all epitopes made from seq
    """
    epitopes = []
    i = 0

    # while we can slice the seq
    while i <= len(seq):

        # making new epitope and removing stop symbols
        ep = ''.join(seq[i:i + length]).replace('*', '')
        epitopes.append(ep)
        i += length
    return epitopes

