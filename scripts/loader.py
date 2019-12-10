import shutil
import urllib3
import certifi
import os


class Loader:

    def __init__(self):
        self.http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
        self.patients = [f'p{i}' for i in range(1, 12)]
        self.regions = ["V3", "PR", "psi", "vpr", "vpu", "p1", "p2", "p6", "p7", "p15", "p17", "RRE"]

    def load_haplotypes(self, folder='data/hivevo'):

        if not os.path.isdir(folder):
            os.mkdir(folder)
        api = "https://hiv.biozentrum.unibas.ch/api/data/haplotypes/{}/{}"
        for patient in self.patients:
            for region in self.regions:
                url = api.format(patient, region)

                filename = f'{folder}/hivevo_{patient}_{region}.fasta'

                with self.http.request('GET', url, preload_content=False) as res, open(filename, 'wb') as out_file:
                    shutil.copyfileobj(res, out_file)

    def load_reference(self, folder='data/hivevo'):
        if not os.path.isdir(folder):
            os.mkdir(folder)

        api = "https://hiv.biozentrum.unibas.ch/api/data/referenceSequence/{}"
        for patient in self.patients:
            url = api.format(patient)
            filename = f'{folder}/hivevo_reference_{patient}.fasta'
            with self.http.request('GET', url, preload_content=False) as res, open(filename, 'wb') as out_file:
                shutil.copyfileobj(res, out_file)

