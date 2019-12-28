# TODO: this script needs to be replaced by the one from master. All references in your code should be changed to Loader

import shutil
import os
import urllib3
import certifi
from Bio import SeqIO
import json
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re

# HIV properties

patients = ["p{}".format(i) for i in range(1, 12)]
hiv_regions = ["V3", "PR", "psi", "vpr", "vpu", "p1", "p2", "p6", "p7", "p15", "p17", "RRE"]


def download_hivevo_haplotype(patient, hiv_region, folder):
    """
    Making folder for data and downloading from https://hiv.biozentrum.unibas.ch alignment for patient for exact region

    Args:
        patient: string (name of patient)
        hiv_region: string (name of region)
        folder: path (name of folder to which download)

    """

    if patient not in patients:
        raise Exception('Wrong name of patient, you can see patients at https://hiv.biozentrum.unibas.ch')

    if hiv_region not in hiv_regions:
        raise Exception('Wrong region, you can see all regions at https://hiv.biozentrum.unibas.ch')

    http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED',
                               ca_certs=certifi.where())

    api = "https://hiv.biozentrum.unibas.ch/api/data/haplotypes/"

    url = "/".join((api, patient, hiv_region))

    if not os.path.isdir(folder):
        os.mkdir(folder)

    file_name = folder + "_".join(("hivevo", patient, hiv_region)) + ".fasta"

    with http.request('GET', url, preload_content=False) as res, open(file_name, 'wb') as out_file:
        shutil.copyfileobj(res, out_file)


def download_hivevo_references(folder):
    """
    Downloading references for all patients

    Args:
        folder: str, path to folder to download
    """
    global patients

    http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED',
                               ca_certs=certifi.where())

    if not os.path.isdir(folder):
        os.mkdir(folder)

    for patient in patients:
        api = "https://hiv.biozentrum.unibas.ch/api/data/referenceSequence"
        url = "/".join((api, patient))
        file_name = folder + "_".join(("hivevo", "reference", patient)) + ".fasta"

        with http.request('GET', url, preload_content=False) as res, open(file_name, 'wb') as out_file:
            shutil.copyfileobj(res, out_file)


def extracting_region_from_reference(region, reference_path, folder):
    """

    Args:
        region: str, name of the region
        reference_path: str, path to reference file
        folder: str, path to folder to store

    Returns:

    """

    if region not in hiv_regions:
        raise Exception('Wrong region, you can see all regions at https://hiv.biozentrum.unibas.ch')

    if not os.path.isdir(folder + region):
        os.mkdir(folder + region)

    with open(reference_path) as f:
        reference_info = json.load(f)
    for reg in reference_info['features']:
        if reg['name'] == region:
            loc = reg['location'][0]
            break
        else:
            continue
    patient = re.search(r"p[\d]*", reference_path)[0]
    # print(patient)
    sequence = reference_info['seq'][loc[0]:loc[1]]
    # print(sequence)
    res = r'/' + re.search(r'[\w]*\.fasta', reference_path)[0].replace(patient, "_".join((patient, region)))

    with open(folder + region + res, 'w') as fasta_file:
        line_1 = reference_info['description'].replace('genomewide', 'region=' + region).lstrip()
        # print(line_1)
        line_2 = sequence
        seq = SeqRecord(Seq(line_2), id='', description=line_1)
        SeqIO.write(seq, fasta_file, 'fasta')


def json2fasta(folder, path_json):
    """
    Сonverting json to fasta

    Args:
        folder: str, path to folder to save
        path_json: str, path to json file

    """
    if not os.path.isdir(folder + 'fasta'):
        os.mkdir(folder + 'fasta')

    with open(path_json) as f:
        json_file = json.load(f)

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

    with open(path_fasta, 'a') as fasta_file, open(path_ref_region) as ref_fasta:
        record = next(SeqIO.parse(ref_fasta, 'fasta'))
        SeqIO.write(record, fasta_file, 'fasta')


def read_fasta(path):
    """
    Read fasta
    Args:
        path: str, path to fasta file

    Returns:
        list of dicts and list of day:
            dict in list of dicts have 'desc' and 'seq' keys
            list of days contain all days of HIV seq we know for exact patient
    """
    haplo_seq_dict = []
    days = []
    patt = r'_[\d]*_'  # for days
    records = []
    with open(path) as fasta_file:  # Reading
        for record in SeqIO.parse(fasta_file, 'fasta'):
            records.append(record)

    haplo_seq_dict = []
    for record in records:
        haplo_seq_dict.append({})
        if 'reference' in record.description:
            haplo_seq_dict[-1]['seq'] = str(record.seq)
            haplo_seq_dict[-1]['desc'] = 'reference'
        # print('1')
        if 'reference' not in record.description:  # Making days, lst
            days.append((re.search(patt, record.description).group(0).replace('_', '')))
            haplo_seq_dict[-1]['seq'] = str(record.seq)
            haplo_seq_dict[-1]['desc'] = record.description.lstrip()

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
