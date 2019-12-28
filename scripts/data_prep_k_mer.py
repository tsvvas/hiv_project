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
    # TODO: correct types in docsting.
    #  Return data type should be specified, example:
    #  type: description
    #  No names for return args are needed unless there are many of them (more than one).
    #  Please check out https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html
    #  cmd + F: function_with_types_in_docstring
    """
        Counting all ocurrances of substring in string using find() method

    Args:
        string: str, string where to find
        sub: str, string to find

    Returns:
        number of ocurrances
    """
    # TODO: comment
    count = start = 0

    # TODO: comment
    while True:
        start = string.find(sub, start) + 1

        # TODO: comment
        if start > 0:
            count += 1

        # TODO: comment
        else:
            return count


# TODO: what is this? why is it in the middle of code?
amino_string = 'ARNDCQEGHILKMFPSTWYV'


def making_aa_k_mers(k):
    # TODO: docstring
    """
    Making all subsequences with length k using aminoacids
    Args:
        k: int, length of k-mer

    Returns:
        list of all possible aminoacid k-mer sequences
    """

    # TODO: NO GLOBAL VARIABLES
    global amino_string

    # TODO: comment
    subseq_iter = itertools.product(amino_string, repeat=k)

    # TODO: comment. Why list()?
    aa_k_mer_list = list(subseq_iter)

    # TODO: why del?
    del subseq_iter

    # TODO: comment
    for i in range(len(aa_k_mer_list)):

        # TODO: comment
        tup = aa_k_mer_list[i]

        # TODO: comment
        aa_k_mer_list[i] = ''.join(tup)
    return aa_k_mer_list


def read_fasta_file_prep(path):
    # TODO: we do not have fasta yet. Why do you need this?
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

    # TODO: comment
    org_name = path.replace('../data/proteomes/', '').replace('.fasta', '')

    return record_list, org_name


def seqio_data(seq_record):
    # TODO: docstring.
    #  Why is this func neccessary?
    """
    Working with SeqRecord class

    Args:
        seq_record: SeqRecord class from Biopython

    Returns:
        protein name: str, name of protein
        sequence: str, sequence of protein
    """
    protein = seq_record.name

    # TODO: why do you need str()?
    seq = str(seq_record.seq)

    return protein, seq


def finding_freq_single_protein(seq, aa_k_mer_list):
    # TODO: docstring
    """
    Finding frequnces for subsequences in single protein
    and scaling it with SKlearn StandardScaler()

    Args:
        seq: str, sequence of amino acids in protein
        aa_k_mer_list: lst, all possible k-mers for aminoacids

    Returns:
        vector_freq: list, frequency of all k-mers from aa_k_mer_list,  vector is normalized using sklearn
    """
    # TODO: comment
    n = len(seq)

    # TODO: this variable is not used
    k = len(aa_k_mer_list[0])

    # TODO: comment
    vector_freq = []

    # TODO: comment
    for x in aa_k_mer_list:
        # TODO: comment, why to use float()?
        vector_freq.append(float(occurrences_count(seq, x)) / n)

    # TODO: comment
    vector_freq = np.array(vector_freq)

    # TODO: comment
    vector_freq = vector_freq.reshape((-1, 1))

    # TODO: comment
    scaler = preprocessing.StandardScaler()
    vector_freq_scaled = scaler.fit_transform(vector_freq)

    # TODO: comment. Why del?
    del vector_freq

    # TODO: why list()?
    return list(vector_freq_scaled)


def main_analyzes(path, k_mer_num):
    # TODO: docstring. If a func does not return anything Returns can be skipped
    """
    Construct "organism_name".csv with 2-mer analyzes

    Args:
        path: path to file in fasta format with represantative proteome used in analyzes
        k_mer_num: k-mer length

    Returns:
        None
    """
    # creating dir
    # TODO: more info in comment. What dir? What for? Why?
    #  Check mkdir args, it can tolerate exsisting dirs
    if not os.path.isdir('csv_data'):
        os.mkdir('csv_data')

    # TODO: no floating comments
    # initializing aa_subseqs and DataFrame

    # TODO: comment
    aa_k_mer_list = making_aa_k_mers(k_mer_num)

    # TODO: comment
    table_columns = ['Organism', 'Protein'] + aa_k_mer_list

    # TODO: comment
    proteins_data = pd.DataFrame(columns=table_columns)

    # reading
    # TODO: more info. Reading what? What for?
    prot_records, organism_name = read_fasta_file_prep(path)

    # dealing with human
    # TODO: comment. Why only human? Can't u drop non humans in read func?
    if organism_name == 'human_proteome':

        # TODO: comment
        human_list = []

        # TODO: comment
        prot_records_split = np.array_split(prot_records, 100)

        # TODO: comment, x is already taken, make understandable variable names!
        for x in prot_records_split:
            human_list.append(x)

        # TODO: comment
        for j in tqdm.tqdm_notebook(range(0, 100)):

            # TODO: to what this comment is related to?
            # as human proteome is very thicc we need to split it into pieces to make analysis faster

            # Creating pd.df
            # TODO: more info, correct naming
            Proteins_data = pd.DataFrame(columns=table_columns)

            # TODO: comment
            index = 0

            # TODO: why debug code is here?
            # print(type(j))

            # TODO: comment
            #  probably it is more convenient to add data to df using dict. Please check this variant and try to
            #  implement it.
            for i in range(len(human_list[j])):

                # TODO: comment, more info (what metrics, why).
                #  Filling maybe?
                # calculating metrics and fullfying table
                SeqRecord = human_list[j][i]

                # TODO: comment
                prot_name, seq = seqio_data(SeqRecord)

                # TODO: comment
                freq_vector = finding_freq_single_protein(seq, aa_k_mer_list)

                # TODO: no floating comment
                # making row for table

                # TODO: try to change this (see above the loop)
                adding_row = []
                adding_row.append(organism_name)
                adding_row.append(prot_name)
                adding_row += freq_vector
                Proteins_data.loc[index] = adding_row
                index += 1

            # Writing file
            # TODO: more info in comment. What file? What for?
            writing_path = 'csv_data/' + organism_name + '_' + str(j) + '.csv'
            Proteins_data.to_csv(writing_path)

        # TODO: if you need to finish a func just write 'return', it will return None by default
        return None

    # TODO: comment
    index = 0

    # prot_records stuff
    # TODO: this is not an informative comment! Do the same thing as above for humans: change this iterations to
    #  filling df from dict
    for i in tqdm.tqdm_notebook(range(len(prot_records))):

        # reading
        # TODO: more info in comment. What are u reading? Why?
        seq_record = prot_records[i]

        # TODO: comment
        prot_name, seq = seqio_data(seq_record)

        # calculating metrics
        # TODO: more info
        freq_vector = finding_freq_single_protein(seq, aa_k_mer_list)

        # making row for pandas
        # TODO: try to change
        adding_row = []
        adding_row.append(organism_name)
        adding_row.append(prot_name)
        adding_row += freq_vector
        proteins_data.loc[index] = adding_row
        index += 1

    # Writing file
    # TODO: more info
    writing_path = 'csv_data/' + organism_name + '.csv'
    proteins_data.to_csv(writing_path)

    # TODO: why do u even need this?
    del prot_records

    return None


#
if __name__ == 'main':
    if not os.path.isdir('../data/proteomes'):
        # TODO: you can write strings in multiple lines, no need for + or whatever. Change this one
        raise FileExistsError('Have you downloaded all needed represantative proteomes (or added yours maybe) and stored'
                              + 'them in "data/proteomes" folder?')

    # TODO: comment
    data_files = [f for f in os.listdir('../data/proteomes') if os.path.isfile(os.path.join('../data/proteomes', f))]

    # TODO: comment
    files_path = []
    for i in range(len(data_files)):

        # TODO: comment, f-strings
        files_path.append('../data/proteomes/' + data_files[i])

    # TODO: comment, understandable names
    for x in files_path:

        # TODO: comment. main_analyzes function returns None - what is it doing here?
        tmp = main_analyzes(x, 2)

