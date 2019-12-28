
from Bio import Align


class Node:
    """
    TODO: docstring
    """

    def __init__(self, name, seq, level, parent=None, kids=None):
        """
        TODO: WRITE A COMPLETE DOCSTRING WITH CLEAR DESCRIPTIONS OF ITS ARGUMENTS
        TODO: ARGS SHOULD LOOK LIKE:
        TODO: arg1 (its type): description
        TODO: user classes as types should be referenced as:
        TODO: arg2 (:class:`~path_from_project_dir_to_script_separated_by_dots.Class_name`)
        TODO: example: arg3 (:class: `~scripts.phylo_tree.Node`)

        TODO: set a correct docstring format in your PyCharm (see instruction)
        TODO: PyCharm can form a docsting on its own, all you need to do is to start typing double quotes
        TODO: (", not ') three times. Than just press Enter and the pointer will go to the next line and an automatic
        TODO: docstring will appear below.
        TODO: write a description of logic of a function on the first line after quotes and fill everything for the
        TODO: arguments.

        TODO: classes also need to have a docsting.

        TODO: set inspection rules (will send later)

        Args:
            name:
            seq:
            level:
            parent:
            kids:
        """
        # TODO: remove trailing comments
        #  all comments should be placed above the line reffered in the comment
        self.parent = parent  # node type

        # TODO: why not to initialize kids with empty list in the first place rather than passing it as an argument
        #  later in code? See line 158
        self.kids = kids
        self.name = name
        self.level = level
        self.seq = seq

    def __repr__(self):
        """
        # TODO: docstring
        Returns:

        """

        # TODO: f-strings
        return 'Node for haplo: ' + self.name

    def find_parent(self, levels_lst, levels_nodes):
        """
        # TODO: docstring
        Finding parent from lower level for node
        Args:
            levels_lst:
            levels_nodes:

        Returns:

        """

        # TODO: comment. Also why do you need to specify inf type? There is no need in that
        score_max = -float('inf')

        # TODO: comment. Explain what the index is, why - 1?
        low_level = levels_lst[levels_lst.index(self.level) - 1]

        # TODO: comment. What is this?
        low_level_nodes = levels_nodes[low_level]

        # TODO: comments should be placed above
        aligner = Align.PairwiseAligner()  # Aligner to use

        # TODO: comments should be placed above
        for node in low_level_nodes:  # Finding max score for pairs nodes from low_level (parent nodes) with exact node

            # TODO: comment.
            score = aligner.score(self.seq, node.seq)

            # TODO: comment
            if score > score_max:
                self.parent = node
                score_max = score

            # TODO: this part is redundant: if the condition is not met than nothing happens anyway
            else:
                continue

        # TODO: comment
        self.parent.kids.append(self)

    def _make_path(self):
        # TODO: change quotes
        '''
        Making path from final level to reference
        '''

        if self.level == 'root':
            return self.name
        else:
            # TODO: f-strings
            #  why protected method is used? If its not ment to be 'private', read naming conventions, please.
            return self.name + ' ' + self.parent._make_path()

    def path(self):
        # TODO: quotes
        '''
        Making _make_path more useful
        '''

        # TODO: comment. Why this is so complicated? Maybe you can preprocess this in make_path?
        return self._make_path().split()[::-1]


class Phylo_tree:
    # TODO: docstring

    def __init__(self, root=None):
        # TODO: docstring
        #  can there be a tree without a root? Does 'root' arg need to be a keyword parameter?
        self.root = root
        self.levels_nodes = None
        self.levels_lst = None
        self.levels_days = None

    def construct_levels(self, days):
        # TODO: quotes
        '''
        Fixing construction of our tree
        '''

        # TODO: f-strinngs, comment
        self.levels_days = {'level_' + str(i): days[i] for i in range(len(days))}

        # TODO: comment
        self.levels_lst = list(self.levels_days.keys())

        # TODO: comment
        self.levels_nodes = {x: [] for x in self.levels_lst}

        # TODO: comment
        self.levels_nodes['root'] = None

    def constructing_tree(self, seq_dict, seq_name_days_dict):
        # TODO: docsting, quotes
        '''
        Constructing tree using prepared data
        '''

        # TODO: comment
        root_node = Node(name='reference', seq=seq_dict['reference'], level='root', kids=[])
        self.levels_nodes['root'] = root_node

        # TODO: split code into logical blocks with new lines
        #  comments above
        for level in self.levels_lst:  # Iterating levels and nodes to find structure of the tree

            # TODO: comments above
            #  r u sure level 0 needs to be processed seperately?
            if level == 'level_0':  # level_0 should be analyzed differently

                # TODO: comment. Do not forget to explain the key
                for haplotype in seq_name_days_dict[self.levels_days[level]]:

                    # TODO: comment
                    node = Node(name=haplotype, seq=seq_dict[haplotype], level=level, parent=root_node)
                    root_node.kids.append(node)
                    self.levels_nodes[level].append(node)
            else:

                # TODO: comment
                for haplotype in seq_name_days_dict[self.levels_days[level]]:
                    node = Node(name=haplotype, seq=seq_dict[haplotype], level=level)

                    # TODO: if u separate loop by level, explain why you use find_parent here and can not use the same
                    #  func for level 0?
                    node.find_parent(self.levels_lst, self.levels_nodes)
                    self.levels_nodes[level].append(node)

    def create_all_path(self):
        # TODO: quotes, docsting
        '''
        Making all paths from final nodes to reference
        '''
        all_paths = []

        # TODO: comment
        for node in self.levels_nodes[self.levels_lst[-1]]:
            all_paths.append(node.path())
        return all_paths
