from Bio import Align

class Node:
    def __init__(self, parent=None, kids=[], name=None, seq=None, level=None):
        """

        """
        self.parent = parent  # node type
        self.kids = kids
        self.name = name
        self.level = level
        self.seq = seq

    def __repr__(self):
        return 'Node for haplo: ' + self.name

    def find_parent(self, levels_lst, levels_nodes):
        '''
        Finding parent from lower level for node
        '''
        score_max = -float('inf')

        low_level = levels_lst[levels_lst.index(self.level) - 1]
        low_level_nodes = levels_nodes[low_level]

        aligner = Align.PairwiseAligner()  # Aligner to use
        for node in low_level_nodes:  # Finding max score for pairs nodes from low_level (parent nodes) with exact node
            score = aligner.score(self.seq, node.seq)
            if score > score_max:
                self.parent = node
                score_max = score
        self.parent.kids.append(self)

    def make_path(self):
        '''
        Making path from final level to reference
        '''

        if self.level == 'root':
            return self.name
        else:
            return self.name + ' ' + self.parent._make_path()

    def path(self):
        """
        Making _make_path more useful

        Return:
            list, list of all possible paths in this tree
        """
        return self.make_path().split()[::-1]


class Phylo_tree:
    def __init__(self, root=None):
        self.root = root
        self.levels_nodes = None
        self.levels_lst = None
        self.levels_days = None

    def construct_levels(self, days):
        '''
        Fixing construction of our tree
        '''
        self.levels_days = {'level_' + str(i): days[i] for i in range(len(days))}
        self.levels_lst = list(self.levels_days.keys())
        self.levels_nodes = {x: [] for x in self.levels_lst}
        self.levels_nodes['root'] = None

    def constructing_tree(self, seq_dict, seq_name_days_dict):
        '''
        Constructing tree using prepared data
        '''
        root_node = Node(name='reference', seq=seq_dict['reference'], level='root', kids=[])
        self.levels_nodes['root'] = root_node
        for level in self.levels_lst:  # Iterating levels and nodes to find structure of the tree
            if level == 'level_0':  # level_0 should be analyzed differently
                for haplotype in seq_name_days_dict[self.levels_days[level]]:
                    node = Node(name=haplotype, seq=seq_dict[haplotype], level=level, parent=root_node)
                    root_node.kids.append(node)
                    self.levels_nodes[level].append(node)
            else:
                for haplotype in seq_name_days_dict[self.levels_days[level]]:
                    node = Node(name=haplotype, seq=seq_dict[haplotype], level=level)
                    node.find_parent(self.levels_lst, self.levels_nodes)
                    self.levels_nodes[level].append(node)

    def create_all_path(self):
        '''
        Making all paths from final nodes to reference
        '''
        all_paths = []
        for node in self.levels_nodes[self.levels_lst[-1]]:
            all_paths.append(node.path())
        return all_paths