from Bio.Align import PairwiseAligner


class Tree:

    def __init__(self, data):
        self.start_day = data.days.min()
        self.end_day = data.days.max()
        self.data = data
        self.region = data.name.unique()[0]
        self.root = None
        self.mapping = {tuple(x): y for (x, y) in
                        zip(data.loc[:, ['days', 'sequence']].sort_values(by='days').itertuples(index=False),
                            range(len(data)))}
        self.graph = {}

    def build(self):
        days = sorted(self.data.days.unique(), reverse=True)

        raw_tree = dict.fromkeys(days, {})

        for i, day in enumerate(days[:-1]):

            current_df = self.data.loc[self.data.days == day]

            previous_day = days[i + 1]
            previous_df = self.data.loc[self.data.days == previous_day]

            for row in current_df.itertuples(index=False):
                if row.sequence not in raw_tree[day]:
                    raw_tree[day][row.sequence] = Node(row.sequence, row.frequency, row.translated)

                # here it is
                parent = raw_tree[day][row.sequence].align(previous_df.loc[:, ['sequence', 'frequency']])

                if parent not in raw_tree[previous_day]:
                    parent_node = previous_df.loc[previous_df.sequence == parent]
                    raw_tree[previous_day][parent] = Node(parent, parent_node.frequency.values[0],
                                                          parent_node.translated.values[0])

                raw_tree[previous_day][parent].children.append(raw_tree[day][row.sequence])
                raw_tree[day][row.sequence].parent = raw_tree[previous_day][parent]

                self.graph[(self.mapping[(previous_day, parent)],
                            self.mapping[(day, raw_tree[day][row.sequence].seq)])] = raw_tree[day][row.sequence].score

        self.root = raw_tree[previous_day][parent]

class Node:

    def __init__(self, seq, freq, translated):
        self.parent = None
        self.children = []
        self.alignment_score = None
        self.seq = seq
        self.frequency = freq
        self.parent_scores = {}
        self.score = None
        self.translation = translated
        self.features = None

    def pairwise(self, potential_parent):
        aligner = PairwiseAligner()
        aligner.mode = 'local'
        seq = self.seq.replace('-', '_')
        potential_parent_undersores = potential_parent.replace('-', '_')
        score = aligner.score(seq, potential_parent_undersores)
        self.parent_scores[potential_parent] = score

    def align(self, parents):
        self.parent_scores = {}

        # parents.sequence and parents.frequency
        for row in parents.itertuples():
            # print(row.frequency, self.frequency)
            # if row.frequency >= self.frequency:
            self.pairwise(row.sequence)
        return self.find_parent()

    def find_parent(self):
        sorted_seqs = sorted(self.parent_scores, key=self.parent_scores.get, reverse=True)
        if sorted_seqs:
            self.score = self.parent_scores[sorted_seqs[0]]
            return sorted_seqs[0]
        print('Cannot find a parent')

    def __repr__(self):
        return f'Node {self.seq}, freq {self.frequency}, has {len(self.children)} children'

