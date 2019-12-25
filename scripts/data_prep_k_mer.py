import numpy as np
from sklearn import preprocessing
import warnings

warnings.filterwarnings("ignore")
from Bio import SeqIO
import pandas as pd
import tqdm
import os
import itertools


def occurrences_count(string, sub):
    """
        Counting all ocurrances of substring in string using find() method

    Args:
        string: str, string where to find
        sub: str, string to find

    Returns:
        number of ocurrances
    """
    count = start = 0
    while True:
        start = string.find(sub, start) + 1
        if start > 0:
            count += 1
        else:
            return count


amino_string = 'ARNDCQEGHILKMFPSTWYV'


def making_aa_k_mers(k):
    """
    Making all subsequences with length k using aminoacids
    Args:
        k: int, length of k-mer

    Returns:
        list of all possible aminoacid k-mer sequences
    """
    global amino_string

    subseq_iter = itertools.product(amino_string, repeat=k)
    aa_k_mer_list = list(subseq_iter)
    del subseq_iter
    for i in range(len(aa_k_mer_list)):
        tup = aa_k_mer_list[i]
        aa_k_mer_list[i] = ''.join(tup)
    return aa_k_mer_list


def read_fasta_file_prep(path):
    """
    Reading fasta file

    Args:
        path: path to file to read

    Returns:
        record_list: list, list (whole organism) of lists (for each protein) of SeqRecords
        org_name: str, name of particular organism
    """
    fasta_test_file = SeqIO.parse(path, 'fasta')
    record_list = list(fasta_test_file)
    org_name = path.replace('../data/proteomes/', '').replace('.fasta', '')

    return record_list, org_name


def seqio_data(seq_record):
    """
    Working with SeqRecord class

    Args:
        seq_record: SeqRecord class from Biopython

    Returns:
        protein name: str, name of protein
        sequence: str, sequence of protein
    """
    protein = seq_record.name
    seq = str(seq_record.seq)

    return protein, seq


def finding_freq_single_protein(seq, aa_k_mer_list):
    """
    Finding frequnces for subsequences in single protein
    and scaling it with SKlearn StandardScaler()

    Args:
        seq: str, sequence of amino acids in protein
        aa_k_mer_list: lst, all possible k-mers for aminoacids

    Returns:
        vector_freq: list, frequency of all k-mers from aa_k_mer_list,  vector is normalized using sklearn
    """
    n = len(seq)
    k = len(aa_k_mer_list[0])

    vector_freq = []

    for x in aa_k_mer_list:
        vector_freq.append(float(occurrences_count(seq, x)) / n)

    vector_freq = np.array(vector_freq)
    vector_freq = vector_freq.reshape((-1, 1))
    scaler = preprocessing.StandardScaler()
    vector_freq_scaled = scaler.fit_transform(vector_freq)
    del vector_freq

    return list(vector_freq_scaled)


def main_analyzes(path, k_mer_num):
    """
    Construct "organism_name".csv with 2-mer analyzes

    Args:
        path: path to file in fasta format with represantative proteome used in analyzes
        k_mer_num: k-mer length

    Returns:
        None
    """
    # creating dir
    if not os.path.isdir('csv_data'):
        os.mkdir('csv_data')

    # initializing aa_subseqs and DataFrame

    aa_k_mer_list = making_aa_k_mers(k_mer_num)
    table_columns = ['Organism', 'Protein'] + aa_k_mer_list
    proteins_data = pd.DataFrame(columns=table_columns)

    # reading
    prot_records, organism_name = read_fasta_file_prep(path)

    # dealing with human
    if organism_name == 'human_proteome':
        human_list = []
        prot_records_split = np.array_split(prot_records, 100)
        for x in prot_records_split:
            human_list.append(x)
        for j in tqdm.tqdm_notebook(range(0, 100)):

            # as human proteome is very thicc we need split it into pieces to make analyzes faster
            # Creating pd.df
            Proteins_data = pd.DataFrame(columns=table_columns)
            index = 0
            # print(type(j))
            for i in range(len(human_list[j])):
                # calculating metrics and fullfying table
                SeqRecord = human_list[j][i]
                prot_name, seq = seqio_data(SeqRecord)
                freq_vector = finding_freq_single_protein(seq, aa_k_mer_list)

                # making row for table

                adding_row = []
                adding_row.append(organism_name)
                adding_row.append(prot_name)
                adding_row += freq_vector
                Proteins_data.loc[index] = adding_row
                index += 1

            # Writing file
            writing_path = 'csv_data/' + organism_name + '_' + str(j) + '.csv'
            Proteins_data.to_csv(writing_path)

        return None

    index = 0

    # prot_records stuff
    for i in tqdm.tqdm_notebook(range(len(prot_records))):
        # reading
        seq_record = prot_records[i]
        prot_name, seq = seqio_data(seq_record)

        # calculating metrics
        freq_vector = finding_freq_single_protein(seq, aa_k_mer_list)

        # making row for pandas
        adding_row = []
        adding_row.append(organism_name)
        adding_row.append(prot_name)
        adding_row += freq_vector
        proteins_data.loc[index] = adding_row
        index += 1

    # Writing file
    writing_path = 'csv_data/' + organism_name + '.csv'
    proteins_data.to_csv(writing_path)

    del prot_records

    return None


if __name__ == 'main':
    data_files = [f for f in os.listdir('../data/proteomes') if os.path.isfile(os.path.join('../data/proteomes', f))]
    files_path = []
    for i in range(len(data_files)):
        files_path.append('../data/proteomes/' + data_files[i])

    for x in files_path:
        tmp = main_analyzes(x, 2)
