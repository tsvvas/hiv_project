import numpy as np
from sklearn import preprocessing
import warnings
warnings.filterwarnings("ignore")
from Bio import SeqIO
import pandas as pd
import tqdm
import os
import itertools
import work_with_files


def occurrences_count(string, sub):
    """
        Counting all ocurrances of substring in string using find() method

    Args:
        string: str, string where to find
        sub: str, string to find

    Returns:
        number of ocurrances: int
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
    Making all possible subsequences with length k using aminoacids (order is important)
    Args:
        k: int, length of k-mer

    Returns:
        list of all possible aminoacid k-mer sequences
    """
    global amino_string

    # making all possible substrings

    subseq_iter = itertools.product(amino_string, repeat=k)
    aa_k_mer_list = list(subseq_iter)
    del subseq_iter

    # one "for" to deal with tuples which we get from itertools stuff

    for i in range(len(aa_k_mer_list)):
        tup = aa_k_mer_list[i]
        aa_k_mer_list[i] = ''.join(tup)

    return aa_k_mer_list

def seqio_data(seq_record):
    """
    Working with SeqRecord class

    Args:
        seq_record: SeqRecord class from Biopython

    Returns:
        protein name: str, name of protein
        sequence: str, sequence of protein
    """

    # getting protein name
    protein = seq_record.name

    # getting sequence
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
    !!!Warning!!! can take much time, so be prepared and sure that all parameters are good
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
    prot_records, organism_name = work_with_files.read_fasta(path)

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

    return None


def main(k):
    # helpful words
    if not os.path.isdir('../data/proteomes'):
        raise FileExistsError('Have you download all needed represantative proteomes (or added yours maybe) and stored'
                              + 'them in "data/proteomes" folder?')

    # making list with file names
    data_files = [f for f in os.listdir('../data/proteomes') if os.path.isfile(os.path.join('../data/proteomes', f))]
    files_path = []
    for i in range(len(data_files)):
        files_path.append('../data/proteomes/' + data_files[i])

    # using main_analyzes for all files
    for x in files_path:
        tmp = main_analyzes(x, k)

if __name__ == 'main':
    main(2)
