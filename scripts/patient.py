import pandas as pd
import json
import os


class Patient:

    def __init__(self, pid, data_path='data/hivevo'):
        self.id = pid
        self.reference = Reference(data_path, patient=pid)
        self.status = None
        self.data_path = data_path
        self.regions = None
        self.init_regions()

    def init_regions(self):
        regions = pd.DataFrame()

        for name in ["V3", "PR", "psi", "vpr", "vpu", "p1", "p2", "p6", "p7", "p15", "p17", "RRE"]:
            filename = f'hivevo_{self.id}_{name}.fasta'

            with open(os.path.join(self.data_path, filename)) as f:
                region_file = json.load(f)

            if type(region_file) == dict:
                print(f'No haplotype for patient {self.id} for region {name}')
                continue

            region = pd.DataFrame(data=region_file)
            region.columns = ["days", "frequency", "sequence", "name", "description"]
            region['nreads'] = region.description.str.split(': ').str[-1]
            region['name'] = name
            region.drop(['description'], axis=1, inplace=True)
            regions = pd.concat([regions, region], ignore_index=True)

        self.regions = Region(regions)


class Reference:
    """ reference data """

    def __init__(self, reference_path, patient=None):

        self.reference_df = pd.DataFrame()

        # all files
        reference_list = [x for x in os.listdir(reference_path) if 'reference' in x]

        # if 1 patient is needed
        if patient:
            reference_list = [x for x in reference_list if f'_{patient}.' in x]

        # compiling from jsons
        for name in reference_list:
            with open(os.path.join(reference_path, name)) as f:
                json_file = json.load(f)
                t = pd.DataFrame(data=json_file).drop(['name', 'description'], axis=1)
                t = pd.concat([t.drop(['features'], axis=1), t.features.apply(pd.Series)], axis=1)
                t.location = t.location.apply(pd.Series)
                t['region_seq'] = t.apply(lambda x: x.seq[x.location[0]:x.location[1]], axis=1)
                t.rename(mapper={'region_seq': 'sequence', 'seq': 'full_reference'}, axis=1, inplace=True)
                self.reference_df = pd.concat([self.reference_df, t], ignore_index=True)

        # if 1 patient, assigning a region class to property
        if patient:
            self.region = Region(self.reference_df)

    def patient(self, patient, region=None):
        """ access to random keys """

        if region:
            return self.reference_df.loc[(self.reference_df.id == f'reference_{patient}') &
                                         (self.reference_df.name == region)]
        else:
            return self.reference_df.loc[self.reference_df.id == f'reference_{patient}']


class Region:

    def __init__(self, data):
        self.region = data
        self.keys = data.name.unique()
        self.score()

    def __getitem__(self, key):
        # dict like access
        return self.region.loc[self.region.name == key]

    def __len__(self):
        return len(self.keys)

    def __repr__(self):
        # df
        return repr(self.region)

    def score(self):
        pass

    def major(self, region=None):
        if 'days' in self.region.columns:
            if region:
                return self.region.loc[self.region.name == region].groupby(['days']).max()
            return self.region.groupby(['name', 'days']).max()
        return None

