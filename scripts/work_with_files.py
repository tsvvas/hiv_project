# coding=utf-8
import shutil
import os
import urllib3
import certifi
from Bio import SeqIO
import Bio
import json
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re



# HIV regions and patients lists
patients = ["p{}".format(i) for i in range(1, 12)]
hiv_regions = ["V3", "PR", "psi", "vpr", "vpu", "p1", "p2", "p6", "p7", "p15", "p17", "RRE"]



def read_fasta(path):
    """
    Reading fasta file using Biopython

    Args:
        path: path to file to read

    Returns:
        record_list: list, list (whole organism) of lists (for each protein) of SeqRecords
        org_name: str, name of particular organism
    """
    # reading proteins
    fasta_test_file = SeqIO.parse(path, 'fasta')
    record_list = list(fasta_test_file)

    # making org_name
    org_name = re.search(r'\w*.fasta' , path)[0].replace('.fasta', '')

    return record_list, org_name



def extracting_region_from_reference(region, reference_path, folder):
    """
    Function to extract specified region from exact reference
    Create fasta file with sequence for needed region
    Args:
        region: str, name of the region
        reference_path: str, path to reference file
        folder: str, path to folder to store

    Returns:
        None
    """

    # checking region
    if region not in hiv_regions:
        raise Exception('Wrong region, you can see all regions at https://hiv.biozentrum.unibas.ch')

    # making dir for exact region we will store all region references for all possible patients here
    if not os.path.isdir(folder + region):
        os.mkdir(folder + region)

    # opening reference and finding location
    with open(reference_path) as f:
        reference_info = json.load(f)
    for reg in reference_info['features']:
        if reg['name'] == region:
            loc = reg['location'][0]
            break
        else:
            continue

    # finding out what patient we are working with
    patient = re.search(r"p[\d]*", reference_path)[0]

    # finally finding out needed sequence
    sequence = reference_info['seq'][loc[0]:loc[1]]

    # res stuff to create name for our file
    res = r'/' + re.search(r'[\w]*\.fasta', reference_path)[0].replace(patient, "_".join((patient, region)))

    # writing sequence and info into fasta file
    with open(folder + region + res, 'w') as fasta_file:
        line_1 = reference_info['description'].replace('genomewide', 'region=' + region).lstrip()
        line_2 = sequence
        seq = SeqRecord(Seq(line_2), id='', description=line_1)
        SeqIO.write(seq, fasta_file, 'fasta')



def json2fasta(folder, path_json):
    """
    Сonverting json to fasta

    Args:
        folder: str, path to folder to save
        path_json: str, path to json file

    Returns:
        None
    """

    # making folder for fasta files
    if not os.path.isdir(folder + 'fasta'):
        os.mkdir(folder + 'fasta')

    # loading file
    with open(path_json) as f:
        json_file = json.load(f)

    # using Biopython to write loaded stuff from json to fasta
    path = path_json.replace('data/', 'data/fasta/')
    with open(path, "w") as fasta_file:  # Iterating for objects in json
        for obj in json_file:
            line_1 = obj['name'].lstrip()
            line_2 = obj['sequence']  # Using SeqIO to do stuff
            seq = SeqRecord(Seq(line_2), id='', description=line_1)
            SeqIO.write(seq, fasta_file, 'fasta')


def add_ref_reg2fasta(path_fasta, path_ref_region):
    """

    Args:
        path_fasta: path to fasta file what will be changed
        path_ref_region: path to fasta file with needed region

    Returns:
        Adding exact region from reference file to fasta with haplotypes
    """

    # writing reference of region in fasta with haplotypes
    with open(path_fasta, 'a') as fasta_file, open(path_ref_region) as ref_fasta:
        record = next(SeqIO.parse(ref_fasta, 'fasta'))
        SeqIO.write(record, fasta_file, 'fasta')


def read_fasta_haplo(path):
    """
    Read fasta
    Args:
        path: str, path to fasta file

    Returns:
        list of dicts and list of day:
            dict in list of dicts have 'desc' and 'seq' keys
            list of days contain all days of HIV seq we know for exact patient
    """

    # creating lists
    haplo_seq_dict = []
    days = []

    # patt for finding out days using re
    patt = r'_[\d]*_'  # for days

    # reading fasta
    records, org = read_fasta(path)

    # analyzing every record and getting all needed info from it
    haplo_seq_dict = []
    for record in records:
        haplo_seq_dict.append({})

        # if for reference because reference desc have no days pattern in it
        if 'reference' in record.description:
            haplo_seq_dict[-1]['seq'] = str(record.seq)
            haplo_seq_dict[-1]['desc'] = 'reference'

        if 'reference' not in record.description:  # Making days, lst
            days.append((re.search(patt, record.description).group(0).replace('_', '')))
            haplo_seq_dict[-1]['seq'] = str(record.seq)
            haplo_seq_dict[-1]['desc'] = record.description.lstrip()

    # completing days list
    days = list(set(days))
    days = sort_lst(days)

    return haplo_seq_dict, days


def sort_lst(lst):
    """
    Making list, which contain numbers in str format, sorted, but numbers will still be in str

    Args:
        lst: list, list to sort

    Returns:
        list

    """
    return [str(x) for x in sorted([int(x) for x in lst])]


'''if __name__ = '__main__':
    # folder = "data/"
    # for patient in patients:
    # for region in hiv_regions:
    # download_hivevo_haplotype(patient, region, folder)
    # folder = "data/references/"
    # download_hivevo_references(folder)'''
