import re
import os
import numpy as np
import pandas as pd
import networkx as nx
import patients_data
import matplotlib.pyplot as plt

from scipy.stats import mode
from io import StringIO
from Bio.Seq import Seq
from Bio import Phylo, SeqIO, AlignIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.Applications import PhymlCommandline
from Bio.Align.Applications import MuscleCommandline
from Bio.Phylo.TreeConstruction import *

if 'trees' not in next(os.walk('./'))[1]:
    os.mkdir('trees')


class Phylo_tree:
    def __init__(self, patient, region, data='DNA', alignment=None):
        '''
        Class constructor
        
        Parameters
        ----------
        patient : str
            Patient name e.g. "p1".
        
        region : str, 
            Region name e.g. "V3".
        
        data : {"DNA", "protein"}, default: "DNA"
            Type of required data.
        
        alignment : str, optional
            Name of file this aligned haplotype sequenses.
            Note: sequenses names must be in format 
                "patient{patient name}_region{region name}_day{day}_id{id}"
                and all id must be unique!
        '''
        fmt = 'phylip-relaxed'
        self.plot_properties = -1, -1
        self.dir = './trees/'
        self.file_name = f'patient{patient}_region{region}_{data}'
        self.align_file_name = f'{self.file_name}_align.{fmt}'
        self.data_type = data
        if self.align_file_name not in next(os.walk(self.dir))[2] and not alignment:
            pat = patients_data.Patient(patient).regions.region
            ref = patients_data.Reference('data/references/').get_patient(patient)
            pat = pd.concat([pat, ref], ignore_index=True).sort_values(by=['days'])
            pat = pat[pat.name == region]
            if data == 'DNA':
                alphabet = IUPAC.unambiguous_dna
                column = 'sequence'
            elif data == 'protein':
                alphabet = IUPAC.protein
                column = 'translated'
            else:
                raise ValueError('Incorrect format! Format must be "DNA" or "protein".')

            records = []
            for i in pat.index:
                haplo = pat.loc[i]
                name = f'patient{patient}_region{region}_day{haplo.days}_id{i}'
                seq = Seq(haplo[column], alphabet=alphabet) if data == 'DNA' else haplo[column]
                records.append(SeqRecord(seq, id=name, description=data))
            lens = np.vectorize(len)(pat[column])
            p_mode = mode(lens)[0][0]
            if data == 'protein' or not np.all(lens == p_mode):
                if data == 'DNA':
                    print(f'Warning! Sequences are not all the same length!\nmode len = {p_mode}\ndeviations:')
                    for r in np.array(records)[lens != p_mode]:
                        print(f'{r.id}  len = {len(r)}')
                    print('Alignment occurs...')
                handle = StringIO()
                SeqIO.write(records, handle, 'fasta')
                muscle = MuscleCommandline(quiet=True)
                stdout, stderr = muscle(stdin=handle.getvalue())
                if stderr:
                    return stderr
                align = StringIO(stdout)
                SeqIO.convert(align, 'fasta', f'{self.dir}{self.align_file_name}', fmt)
            else:
                SeqIO.write(records, f'{self.dir}{self.align_file_name}', fmt)
        elif alignment:
            try:
                AlignIO.convert(*alignment, f'{self.dir}{self.align_file_name}', fmt)
            except ValueError:
                print(f'Warning! Sequences are not all the same length!\nAlignment occurs...')
                with open(alignment[0], 'r') as handle:
                    muscle = MuscleCommandline(quiet=True)
                    stdout, stderr = muscle(stdin=handle.read())
                    if stderr:
                        return stderr
                    align = StringIO(stdout)
                    SeqIO.convert(align, 'fasta', f'{self.dir}{self.align_file_name}', fmt)

    class Haplo():
        '''
        Haplotype storing class
        '''
        def __init__(self, name, tree):
            '''
            Haplotype storing class constructor
            
            Parameters
            ----------
            name : str
                SeqRecord.id

            tree : Bio.Phylo.Newick.Tree 
                Phylogenetic tree
            '''
            self.name = name
            self.patient = re.search(r'patient(.*?\d+)', name).groups()[0]
            self.haplo_id = int(re.search(r'id(.*?\d+)', name).groups()[0])
            self.day = int(re.search(r'day(.*?\d+)', name).groups()[0])
            self.path = tree.get_path(name)                             # путь по дереву до элемента
            self.rand_for_plot = [0]                                    # идентификатор для сортировки
            self.virtual = 0                                            # виртуальная вершина? (для длинных связей)
            self.x = 0                                                  # координата

    def build_tree(self, p, tree_type='ML', quiet=False):
        '''
        Tree building function.
        
        Построение происходит по следующему принципу:
        - Строим вспомогательное ML или MP дерево с помощью сторонних тулов.
          Листья этого дерева - гаплотипы, которые должны стать вершинами нового
          дерева.
        - Составляем список из листьев вспомогательного дерева, представив их как
          объекты специального класса Haplo
        - Вытаскиваем из списка гаплотипы начиная с самого старшего, для каждого 
          извлеченного гаплотипа среди оставшихся в списке ищем наиболее вероятного
          родителя, минимизирующего специальную метрику.
          
        Метрика для поиска родителя.
        У каждого гаплотипа есть предок, присутствовавший в предыдущий день
        секвенирования. Но вероятность обнаружить гаплотип при однократном 
        секвенировании равна р. То есть вероятность найти предка в предыдущий день -
        р, двумя днями ранее - (1-р)р и т.д., то есть имеет геометрическое 
        распределение. Но из длины пути по вспомогательному дереву нам известно
        "эволюционное расстояние" между всемии гаплотипами. Таким образом вероятность 
        одному гаплотипу быть предком другого в данной модели пропорциональна отношению
        вероятности найти предка через данное количество заборов крови к "эволюционному
        расстоянию" между данными гаплотипами.
          
        Parameters
        ----------
        p : float
            Probability of detecting an existing haplotype during one blood sampling.
        
        tree_type : {"ML", "MP"}, default: "ML"
            Type of phylogenetic tree: maximum likehood or maximum parsimony.
            
        quiet : optional, default: False
            If True, additional information about tree building is not displayed.
        '''
        
        # Строим таблицу по всем гаплотипам, из нее можно будет дергать пути
        self.DF = pd.DataFrame()
        for record in SeqIO.parse(f'{self.dir}{self.align_file_name}', 'phylip-relaxed'):
            rec_name = record.id
            rec_id = int(re.search(r'id(.*?\d+)', rec_name).groups()[0])
            rec_day = int(re.search(r'day(.*?\d+)', rec_name).groups()[0])
            rec_dict = {rec_id : {'day':rec_day, 'seq':record.seq}}
            self.DF = self.DF.append(pd.DataFrame.from_dict(rec_dict, orient='index'))
            
        self.prob = p
        
        # Проверяем наличие построенного ML/MP дерева. Если такогого нет - строим
        if f'{self.file_name}_{tree_type}tree.newick' not in next(os.walk(self.dir))[2]:
            if tree_type == 'ML':
                dt = 'nt' if self.data_type == 'DNA' else 'aa'
                command_tree = PhymlCommandline(input=f'{self.dir}{self.align_file_name}', search='BEST', datatype=dt)
                stdout, stderr = command_tree()
                if stderr:
                    return stderr
                os.rename(f'{self.dir}{self.align_file_name}_phyml_tree.txt', f'{self.dir}{self.file_name}_MLtree.newick')
                with open(f'{self.dir}{self.align_file_name}_phyml_stats.txt', 'r') as stats:
                    self.tree_stats = stats.read()
                if not quiet:
                    print(self.tree_stats)
                os.remove(f'{self.dir}{self.align_file_name}_phyml_stats.txt')
                self.tree = Phylo.read(f'{self.dir}{self.file_name}_MLtree.newick', 'newick')
            elif tree_type == 'MP':
                align = AlignIO.read(f'{self.dir}{self.align_file_name}', 'phylip-relaxed')
                model = 'blastn' if self.data_type == 'DNA' else 'blosum62'
                constructor = DistanceTreeConstructor(DistanceCalculator(model), 'nj')
                nj_tree = constructor.build_tree(align)
                searcher = NNITreeSearcher(ParsimonyScorer())
                constructor = ParsimonyTreeConstructor(searcher, nj_tree)
                self.tree = constructor.build_tree(align)
                Phylo.write(self.tree, f'{self.dir}{self.file_name}_MPtree.newick', 'newick')
            else:
                raise ValueError('Incorrect tree type! Tree must be "ML" or "MP".')
        else:
            self.tree = Phylo.read(f'{self.dir}{self.file_name}_{tree_type}tree.newick', 'newick')
            
            
        self.G = nx.DiGraph()    # граф для нашего дерева
        nodes = []               # список будущих вершин дерева
        
        # Добавляем в nodes все листья референсного ML/MP дерева в виде объектов класса Haplo
        for clade in self.tree.get_terminals():
            nodes.append(self.Haplo(str(clade), self.tree))
        nodes.sort(key=lambda node: node.day)                  # сортируем по дню
        self.days = sorted(set(map(lambda x: x.day, nodes)))   # список дней
        while len(nodes) > 1:
            node = nodes.pop()      # изымаем из списка самый старший гаплотип 
            parents = {}            # его потенциальные родители
            
            # Среди оставшихся гаплотипов:
            for nd in nodes:
                if nd.day != node.day:                    # если предшествует рассматриваемому гаплотипу по времени
                    path = set(node.path) ^ set(nd.path)  # путь по вспомогательному дереву
                    prob = p * (1 - p) ** (self.days.index(node.day) - self.days.index(nd.day) - 1) # вероятность
                    path_len = sum(map(lambda item: item.branch_length, path)) # длина пути по вспомогательному дереву
                    if path_len:
                        parents[nd] = prob / path_len     # сопоставляем значение метрики
                    else:
                        parents[nd] = np.inf
            par, weight = max(parents.items(), key=lambda item: item[1])  # выбираем наиболее вероятного предка
            self.G.add_edge(par, node, weight=weight)                     # добавляем ребро в граф

    def plot(self, fig_size=(10,7), days_scaling=True, virt_nodes=True, iters=1000, topo=0):
        '''
        Tree plotting function
        
        Parameters
        ----------
        fig_size : (int, int), optional, default: (10,7)
            Image size.
        
        days_scaling : bool, optional, default: True
            If False distance between levels becomes proportional to the 
            difference in days, else levels are equidistant.
            
        virt_nodes : bool, optional, default: True
            If True, when connecting vertices at different levels, arrows 
            are drawn between every two levels, else it is drowing only one
            arrow for every two connected nodes.
            
        iters : int, optional, default: 1000
            Number of iterations for constructing an optimally compact tree.
            
        topo : int, optional, default: 0
            Tree configuration number. Used for plot fine tuning.
        '''
                
        def branch_shufle(seed):
            np.random.seed(seed)

            for day in self.days[1:]:
                n = len(self.haplos_by_day[day])
                rand = np.random.choice(np.arange(n), n, replace=0)
                self.haplos_by_day[day].sort(key=lambda x: x.haplo_id)
                for i in range(n):
                    haplo = self.haplos_by_day[day][i]
                    haplo.rand_for_plot = self.edges[self.edges[:, 1] == haplo][:, 0][0].rand_for_plot + [haplo.virtual, rand[i]]
                    
        def descent():
            self.haplos_by_day[self.days[0]][0].x = 0
            for t in range(len(self.days)):
                self.haplos_by_day[self.days[t]].sort(key=lambda x: x.rand_for_plot)
                n = len(self.haplos_by_day[self.days[t]])
                pos_last = -np.inf
                par_shift = 0
                last_virtual = 0
                first = True
                for j in range(n):
                    par = self.haplos_by_day[self.days[t]][j]
                    par.x += par_shift
                    children = sorted(self.G_virt.successors(par), key=lambda x: x.rand_for_plot)
                    m = len(children)
                    if m:
                        for i in range(j):
                            self.haplos_by_day[self.days[t]][i].x += par_shift
                        dist = 1 if children[0].virtual or last_virtual else 1 # тут можно поставить разные значения
                        last_virtual = children[-1].virtual
                        if par.x - (m - 1)/2 < pos_last + dist:
                            shift = (pos_last + dist) - (par.x - (m - 1)/2)
                            par.x += shift
                            par_shift += shift
                            if first:
                                first = False
                                for i in range(j):
                                    self.haplos_by_day[self.days[t]][i].x += par_shift
                        pos_0 = par.x - (m - 1)/2
                        pos_last = pos_0 + (m - 1)
                        for i in range(m):
                            children[i].x = pos_0 + i

        def climb():
            for day in self.days[:-1][::-1]:
                shift = 0
                first = True
                for i in range(len(self.haplos_by_day[day])):
                    haplo = self.haplos_by_day[day][i]
                    haplo.x += shift
                    children = sorted(self.G_virt.successors(haplo), key=lambda x: x.rand_for_plot)
                    m = len(children)
                    if m:
                        x = (children[0].x + children[-1].x)/2
                        shift += x - haplo.x
                        haplo.x = x
                        if first:
                            first = False
                            for j in range(i):
                                self.haplos_by_day[day][j].x += shift

        def minimax():
            x_min, x_max = np.inf, -np.inf
            for node in self.G_virt.nodes:
                x_min = min(x_min, node.x)
                x_max = max(x_max, node.x)
            return x_min, x_max
        
        if (self.prob, iters) != self.plot_properties:
            self.plot_properties = self.prob, iters
            self.G_virt = self.G.copy()
            self.haplos_by_day = {}
            for node in self.G_virt.nodes:
                self.haplos_by_day[node.day] = self.haplos_by_day.get(node.day, []) + [node]

            mult_lev_edges = []
            virt_edges = []
            for edge in self.G_virt.edges:
                n = self.days.index(edge[1].day) - self.days.index(edge[0].day)
                if n > 1:
                    mult_lev_edges.append(edge)
                    par = edge[0]
                    for i in range(n-1):
                        child = self.Haplo(par.name, self.tree)
                        child.virtual = 1
                        child.day = self.days[self.days.index(par.day) + 1]
                        child.rand_for_plot = par.rand_for_plot + [0]
                        virt_edges.append((par, child))
                        par = child
                    virt_edges.append((par, edge[1]))

            self.G_virt.add_edges_from(virt_edges)
            self.G_virt.remove_edges_from(mult_lev_edges)

            self.haplos_by_day = {}
            for node in self.G_virt.nodes:
                self.haplos_by_day[node.day] = self.haplos_by_day.get(node.day, []) + [node]

            self.edges = np.array(self.G_virt.edges)
            self.topologys = []
            for s in range(iters):
                branch_shufle(s)
                descent()
                climb()
                x_min, x_max = minimax()
                self.topologys.append((x_max - x_min, s))
            self.topologys.sort()

            self.plot_properties = self.prob, iters
        
        s = self.topologys[topo][1]
        branch_shufle(s)
        descent()
        climb()

        self.nodelist = []
        self.pos = {}
        for node in self.G_virt.nodes:
            if not node.virtual:
                self.nodelist.append(node)
            self.pos[node] = node.x, -self.days.index(node.day) if days_scaling else -node.day
            
        if virt_nodes:
            G_for_plot = self.G_virt
        else:
            G_for_plot = self.G

        x_min, x_max = minimax()

        size = fig_size[0] * 0.8
        shift = (x_max - x_min) * 0.06
        fig_scale = min(*fig_size)

        p = plt.figure(figsize=fig_size)
        plot = nx.draw_networkx(G_for_plot, pos=self.pos, with_labels=0, nodelist=self.nodelist, labels=None, 
                                node_size=fig_scale*12, linewidths=0, width=fig_scale*0.05, arrowsize=fig_scale*0.7)
        plt.text(x_max + shift, 0, 'reference', fontsize=size, horizontalalignment='center')
        for day in self.days[1:]:
            plt.text(x_max + shift, -self.days.index(day) if days_scaling else -day, f'day {day}', 
                     fontsize=size, horizontalalignment='center')
        plt.show()
        
    def paths(self):
        '''
        Function obtaining all the paths in the tree from the 
        first to the last day. 
        
        Returns
        -------
        list of haplo_id lists.
        '''
        ref = self.haplos_by_day[self.days[0]][0]
        paths = []
        for term in self.haplos_by_day[self.days[-1]]:
            paths.append([haplo.haplo_id for haplo in nx.shortest_path(self.G, ref, term)])
        return paths
    
    def get_path(self, path):
        '''
        Function returning a DataFrame of evolutionary path.
        
        Parameters
        ----------
        path : list of int
            List of path's haplo_id.
        
        Returns
        -------
        pandas.DataFrame.
        '''
        return self.DF.loc[path]