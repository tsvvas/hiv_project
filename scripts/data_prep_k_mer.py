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
        int, number of occurrances
    """
    # Starting index and count num
    count = start = 0

    # Counting
    while True:
        # If found count = count + 1
        # start = {where was found} +1

        start = string.find(sub, start) + 1
        if start > 0:
            count += 1
        else:
            return count


def making_aa_k_mers(k):
    """
    Making all possible subsequences with length k using aminoacids (order is important)
    Args:
        k: int, length of k-mer

    Returns:
        list of all possible aminoacid k-mer sequences
    """
    amino_string = 'ARNDCQEGHILKMFPSTWYV'

    # making all possible substrings
    subseq_iter = itertools.product(amino_string, repeat=k)
    aa_k_mer_list = list(subseq_iter)

    # deleting to free space
    del subseq_iter

    # one "for" to deal with tuples which we get from itertools stuff
    for i in range(len(aa_k_mer_list)):

        # tuples are in list
        tup = aa_k_mer_list[i]

        # rewriting list elements
        aa_k_mer_list[i] = ''.join(tup)

    return aa_k_mer_list

def seqio_data(seq_record):
    """
    Working with SeqRecord class

    Args:
        seq_record: SeqRecord class from Biopython

    Returns:
        protein_name: str, name of protein
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
        list, frequency of all k-mers from aa_k_mer_list,  vector is normalized using sklearn
    """

    # Getting initial sizes
    n = len(seq)
    k = len(aa_k_mer_list[0])

    # Initializing list where all frequencies will be saved
    vector_freq = []

    # Counting frequencies
    for x in aa_k_mer_list:
        vector_freq.append(float(occurrences_count(seq, x)) / n)

    # Making some prep with array
    vector_freq = np.array(vector_freq)
    vector_freq = vector_freq.reshape((-1, 1))

    # Standardizing our frequencies
    scaler = preprocessing.StandardScaler()
    vector_freq_scaled = scaler.fit_transform(vector_freq)

    return list(vector_freq_scaled)


def main_analyzes(path, k_mer_num):
    """
    Construct "organism_name".csv with k-mer analyzes
    !!!Warning!!! can take much time, so be prepared and sure that all parameters are good
    Args:
        path: str, path to file in fasta format with represantative proteome used in analyzes
        k_mer_num: int, k-mer length
    """
    # creating dir to store CSVs produced by function
    if not os.path.isdir('csv_data'):
        os.mkdir('csv_data')

    # initializing aa_subseqs
    aa_k_mer_list = making_aa_k_mers(k_mer_num)

    # initializing DataFrame
    table_columns = ['Organism', 'Protein'] + aa_k_mer_list
    proteins_data = pd.DataFrame(columns=table_columns)

    # reading
    prot_records, organism_name = work_with_files.read_fasta(path)

    # dealing with human, because it needs to be analyzed separately
    if organism_name == 'human_proteome':

        # initializing list
        human_list = []

        # appending all human proteins to list and splitting it into 100 parts
        prot_records_split = np.array_split(prot_records, 100)
        for prot_data_part in prot_records_split:
            human_list.append(prot_data_part)

        # We will split analyze of human, because human proteom is too big to handle
        for j in tqdm.tqdm_notebook(range(0, 100)):

            # Creating pd.df
            proteins_data = pd.DataFrame(columns=table_columns)
            index = 0


            for i in range(len(human_list[j])):

                # taking exact protein and calculating metrics (frequencies)
                SeqRecord = human_list[j][i]
                prot_name, seq = seqio_data(SeqRecord)
                freq_vector = finding_freq_single_protein(seq, aa_k_mer_list)

                # making row for table
                adding_row = []
                adding_row.append(organism_name)
                adding_row.append(prot_name)
                adding_row += freq_vector

                # adding row to the DataFrame
                proteins_data.loc[index] = adding_row
                index += 1

            # Writing file for every part of data, we will combine them later
            writing_path = 'csv_data/' + organism_name + '_' + str(j) + '.csv'
            proteins_data.to_csv(writing_path)

        return

    # Rewriting index
    index = 0

    # working with NOT human proteomes
    for i in tqdm.tqdm_notebook(range(len(prot_records))):

        # reading protein to calculate metrics on it
        seq_record = prot_records[i]
        prot_name, seq = seqio_data(seq_record)

        # calculating metrics (frequencies)
        freq_vector = finding_freq_single_protein(seq, aa_k_mer_list)

        # making row for pandas
        adding_row = []
        adding_row.append(organism_name)
        adding_row.append(prot_name)
        adding_row += freq_vector

        # adding row to the DataFrame
        proteins_data.loc[index] = adding_row
        index += 1

    # Writing file
    writing_path = 'csv_data/' + organism_name + '.csv'
    proteins_data.to_csv(writing_path)

    return


def main(k):
    """
    Main functions to work with. Applying full analyze on all
    represantative proteomes you have.
    Args:
        k: int, k-mer length
    """
    # helpful words
    if not os.path.isdir('../data/proteomes'):
        raise FileExistsError('Have you download all needed represantative proteomes (or added yours maybe) and '
                              'stored them in "data/proteomes" folder?')

    # making list with file names in 'proteomes' dir
    data_files = [f for f in os.listdir('../data/proteomes') if os.path.isfile(os.path.join('../data/proteomes', f))]

    # changing a bit paths from data_files in order to make everything work fine
    files_path = []
    for i in range(len(data_files)):
        files_path.append('../data/proteomes/{data_files[i]}')

    # using main_analyzes for all files
    for file in files_path:
        main_analyzes(file, k)

if __name__ == 'main':
    main(2)
