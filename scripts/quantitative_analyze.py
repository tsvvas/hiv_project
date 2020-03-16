import numpy as np
import os
import sys
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import Align
from Bio.Alphabet import IUPAC, Gapped
import matplotlib.pyplot as plt
import itertools
from joblib import load
from igraph import Graph as igraph_graph, plot, rescale

import graph
import patients_data
import tree_building
import data_prep_k_mer

class Quantitative:
    def __init__(self, patients_list):
        self.patients_list = patients_list
        self.patients_evolution = None
        self.classificator = None

    def loading_classificator(self, path):
        # Our classificator
        self.classificator = load(path)

    def clf_metric_2_mer_path(self, path_tree, prot_dict, aa_k_mer_list):
        """
        Converting sequence from path_tree into object that can be used in our classificator

        Args:
            path_tree: list, list with names of haplotypes in exact path in tree
            prot_dict: dict,
            aa_k_mer_list: list, all posible k-mers from aminoacids

        Returns:
            list of all possible aminoacid k-mer sequences
        """
        metric_path = []

        # for-loop for names in path_tree
        for id_ in path_tree:
            metric_path.append([])

            # finding prot for exact name
            prot = prot_dict[id_]

            # finding frequency metrics
            tmp = data_prep_k_mer.finding_freq_single_protein(prot, aa_k_mer_list)

            # KOCTblJLb (our vector from finding_freq_bla_bla is not feeling well (my bad :))
            for x in tmp:
                metric_path[-1].append(float(x[0]))

        return metric_path

    def quantitative_analyzes(self):
        """
        Function to analyze all paths in created Phylo_tree for each patient.
        Here probability to be human's protein is stored in Y, which then will be plotted.
        Args:
            patients: list, list of all patients
        Return:
            patients_evolution: dict, dict of dicts -> {patient:{'X': X, 'Y': [Ys for all paths]}}
        """

        # preparing dict to return
        self.patients_evolution = {}

        # making k-mers
        aa_k_mer_list = data_prep_k_mer.making_aa_k_mers(2)

        # preparing references for all patients
        ref = patients_data.Reference('data/hivevo')

        # for-loop for patients
        for patient in self.patients_list:

            # We will not use patient#3 and patient#10 because their HIV wasn't cool at all
            # joke, additional info can be found here (https://elifesciences.org/articles/11282)

            if patient != 'p3' and patient != 'p10':

                # creating dataset for patient
                pat_class = patients_data.Patient(patient)
                pat_data = pat_class.regions['V3']

                # extracting reference
                ref_data = ref.get_patient(patient, region='V3')

                # adding reference to dataset -> now we are ready to construct tree
                pat_data = pd.concat([ref_data, pat_data], ignore_index=True).sort_values(by=['days'])
                # print(pat_data)

                # Constructing tree
                tree = tree_building.Tree(pat_data)
                tree.build()

                # Seqs data converting
                seq_data = tree.mapping

                prot_dict = {}  # making protein dictionary

                for day_seq in list(seq_data.keys()):
                    id_ = seq_data[day_seq]
                    prot_dict[id_] = Seq(day_seq[1], Gapped(IUPAC.unambiguous_dna)).ungap().translate()

                # Dealing with graph
                vertices = [i for i in range(len(tree.mapping))]
                edges = tree.graph

                g = graph.Graph()
                g.add_vertices(vertices)

                # setting correct weights
                g.set_edge_weights(edges)

                # getting all paths
                phylo_paths = g.all_paths()

                # Creating unique days
                days = set()
                for day, _ in list(tree.mapping.keys()):
                    days.add(day)
                days = sorted(list(days))

                # adding patient
                self.patients_evolution[patient] = {}
                self.patients_evolution[patient]['X'] = None
                self.patients_evolution[patient]['Y'] = []

                # Making X
                self.patients_evolution[patient]['X'] = days

                # Using classificator to find out probability to be human's gene

                for path in phylo_paths:
                    met = self.clf_metric_2_mer_path(path, prot_dict, aa_k_mer_list)
                    Y = self.classificator.predict_proba(met)[:, 1]
                    self.patients_evolution[patient]['Y'].append(Y)

    def plot_paths(self):

        # Making font BIGGER!
        plt.rcParams['font.size'] = 20

        # Plotting
        fig, axs = plt.subplots(3, 3, figsize=(35, 20), sharex='none', sharey='none')

        i = 0
        patients_used = list(self.patients_evolution.keys())
        for x in range(3):
            for y in range(3):
                pat = patients_used[i]
                pat_X = self.patients_evolution[pat]['X']
                pat_Y = self.patients_evolution[pat]['Y']
                for Y in pat_Y:
                    axs[x, y].plot(pat_X, Y)
                    axs[x, y].set_title('patient {}'.format(patients_used[i]))
                axs[x, y].grid(alpha=1)
                i += 1

        # Sadly it should be here because of my newbie skills in pyplot
        axs[0, 0].set(ylabel='Prob to be human prot, %')
        axs[1, 0].set(ylabel='Prob to be human prot, %')
        axs[2, 0].set(ylabel='Prob to be human prot, %')

        axs[2, 0].set(xlabel='time, days')
        axs[2, 1].set(xlabel='time, days')
        axs[2, 2].set(xlabel='time, days')

        fig.suptitle('Evolution of HIV of V3 genome region', fontsize=40)
        plt.show()