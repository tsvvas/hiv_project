import shutil
import urllib3
import certifi
import os
import pandas as pd

class Loader:
    """ Класс для загрузки данных """

    def __init__(self):
        self.http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
        self.patients = [f'p{i}' for i in range(1, 12)]
        self.regions = ["V3", "PR", "psi", "vpr", "vpu", "p1", "p2", "p6", "p7", "p15", "p17", "RRE"]

    def load_haplotypes(self, folder='data/hivevo'):
        """
        Загрузка гаплотипов

        Args:
            folder (str): путь, куда будет произведена загрузка
        """
        
        # если папок не существует, то создаём их
        if not os.path.isdir(folder[:folder.find('/')]):
            os.mkdir(folder[:folder.find('/')])
        
        if not os.path.isdir(folder):
            os.mkdir(folder)

        # адрес запроса
        api = "https://hiv.biozentrum.unibas.ch/api/data/haplotypes/{}/{}"

        # для каждого пациента по регионам запрашиваем данные
        for patient in self.patients:
            for region in self.regions:
                url = api.format(patient, region)

                filename = f'{folder}/hivevo_{patient}_{region}.json'

                # сохраняем их
                with self.http.request('GET', url, preload_content=False) as res, open(filename, 'wb') as out_file:
                    shutil.copyfileobj(res, out_file)
                    
    def load_pcr_stats(self, folder='data/pcr_stats'):
        """
        Downloading statistics about PCR for each patients in .tsv format.

        Args:
            folder (str): path dowload to. Default is 'data/pcr_stats'
        """
        if not os.path.isdir(folder[:folder.find('/')]):
            os.mkdir(folder)

        if not os.path.isdir(folder):
            os.mkdir(folder)

        # addres for request
        api = "https://hiv.biozentrum.unibas.ch/download/depth_{}.tsv"

        # для каждого пациента по регионам запрашиваем данные
        for patient in self.patients:
            # making url
            url = api.format(patient)

            filename = f'{folder}/pcr_stats_{patient}.tsv'

            # saving
            with self.http.request('GET', url, preload_content=False) as res, open(filename, 'wb') as out_file:
                shutil.copyfileobj(res, out_file)         
                    
    def load_bam(self, patient, folder='data/bam_data', folder_stats='data/pcr_stats'):
        """
        Downloading bam files for a single patient.
          
        Note 1: PCR statistics can be download via load_pcr_stats method.
        Note 2: PCR statistics files can't be renamed. Be careful!

        Args:
            patient (str): patient from self.patients.
            folder (str): Default is 'data/raw_data'. In this directory folder for each patient will be created.
            folder_stats (str): Default is 'data/raw_data'. Path where statistics about PCR could be found.
        """
        if not os.path.isdir(folder[:folder.find('/')]):
            os.mkdir(folder)

        if not os.path.isdir(folder):
            os.mkdir(folder)
        
        # lisf of PCR regions
        F_list = [f'F{i}' for i in range(1, 7)]

        # addres for request for bam files
        api = "https://hiv.biozentrum.unibas.ch/download/reads_{}_{}_{}.bam"
        
        # stats template to read
        stats_read = "{}/pcr_stats_{}.tsv"

        # reading tsv for patient
        stats_df = pd.read_csv(stats_read.format(folder_stats, patient), delimiter='\t', encoding='utf-8')

        # finding out sequencing days for patients
        days = list(stats_df['days since infection'])
        
        # folder for patient
        patient_folder = f'{folder}/{patient}'
        
        if not os.path.isdir(patient_folder):
            os.mkdir(patient_folder)

        # for-loops for downloading
        for num, day in enumerate(days):
            for F in F_list:
                # making url. +1 should be added to num, because we need (1, 2, ...) not (0, 1, ...)
                url = api.format(patient, num+1, F)

                # making filename
                filename = f'{patient_folder}/{patient}_{day}_{F}.bam'

                # saving
                with self.http.request('GET', url, preload_content=False) as res, open(filename, 'wb') as out_file:
                    shutil.copyfileobj(res, out_file)

    def load_references(self, folder='data/references'):
        """
        Загрузка референсов

        Args:
            folder (str): путь, куда будет произведена загрузка
        """
        # если папки не существует, то создаём её
        if not os.path.isdir(folder):
            os.mkdir(folder)

        # адрес запроса
        api = "https://hiv.biozentrum.unibas.ch/api/data/referenceSequence/{}"

        # для каждого пациента запрашиваем данные и сохраняем
        for patient in self.patients:
            url = api.format(patient)
            filename = f'{folder}/hivevo_reference_{patient}.json'
            with self.http.request('GET', url, preload_content=False) as res, open(filename, 'wb') as out_file:
                shutil.copyfileobj(res, out_file)
    
    def load_all(self, load_bam=False):
        """
        Download for each patient: haplotypes, pcr statistics, bam files, references using all methods of the class
        with default parameters. Will last about 1-2 hours (mostly because of bam_files). Require 3 GB of empty space.
        This will make 'data' folder with such structure:
        data
            |
            |-> hivevo (folder for haplotypes)
            |-> references (folder for references)
            |-> pcr_stats (folder for PCR statistics)
            |-> bam_files (folder for bam files)
            
        Args:
            load_bam (bool): default=False, if False script wont download bam files, if True it will.
        
        Note: bam files are huge, so be careful.
        """
        self.load_haplotypes()
        self.load_references()
        self.load_pcr_stats()
        
        if load_bam:
            for patient in self.patients:
                self.load_bam(patient)