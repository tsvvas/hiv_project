from Bio.Align import PairwiseAligner
import pandas as pd


class Tree:
    """
    Класс для дерева по дням. Реализует поиск связей между гаплотипами.
    """

    def __init__(self, data):
        """
        Args:
            data (pd.DataFrame): данные о регионе
        """
        self.start_day = data.days.min()
        self.end_day = data.days.max()
        self.data = data
        self.region = data.name.unique()[0]
        self.root = None

        # переименование вершин в номера, т.к. последовательности и дни не влезают
        # также могут быть дублирования в самих последовательностях, т.к. может быть обратная мутация
        # поэтому часть ключа - день
        self.mapping = {tuple(x): y for (x, y) in
                        zip(data.loc[:, ['days', 'sequence']].sort_values(by='days').itertuples(index=False),
                            range(len(data)))}

        # дополнительный маппинг для вывода данных на печать
        self._mapping = {tuple(x): y for (x, y) in zip(data.loc[:, ['days', 'sequence', 'frequency']].sort_values(
            by='days').itertuples(index=False), range(len(data)))}
        self.graph = {}

    def build(self):
        """
        Функция построения дерева
        """
        # получаем последовательные дни
        days = sorted(self.data.days.unique(), reverse=True)

        # инициализируем словарь днями
        raw_tree = dict.fromkeys(days, {})

        # для каждого дня
        for i, day in enumerate(days[:-1]):

            # будем рассматривать данные по текущему дню
            current_df = self.data.loc[self.data.days == day]

            # и данные по предыдущему дню
            previous_day = days[i + 1]
            previous_df = self.data.loc[self.data.days == previous_day]

            # инициализируем узлы
            for row in current_df.itertuples(index=False):
                if row.sequence not in raw_tree[day]:
                    raw_tree[day][row.sequence] = Node(row.sequence, row.frequency, row.translated)

                # найдём родителя текущего листа
                parent = raw_tree[day][row.sequence].align(previous_df.loc[:, ['sequence', 'frequency']])

                # если такого родительского узла ещё нет, то инициализируем его
                if parent not in raw_tree[previous_day]:
                    parent_node = previous_df.loc[previous_df.sequence == parent]
                    raw_tree[previous_day][parent] = Node(parent, parent_node.frequency.values[0],
                                                          parent_node.translated.values[0])

                # обновляем детей узла
                raw_tree[previous_day][parent].children.append(raw_tree[day][row.sequence])

                # делаем ссылку от ребёнка к родителям
                raw_tree[day][row.sequence].parent = raw_tree[previous_day][parent]

                # сохраняем вес выравнивания как вес для графа
                self.graph[(self.mapping[(previous_day, parent)],
                            self.mapping[(day, raw_tree[day][row.sequence].seq)])] = raw_tree[day][row.sequence].score

        # сохраняем корень дерева
        self.root = raw_tree[previous_day][parent]

    def path_info(self, path):
        """
        Возвращает датафрейм с подробностями о пути
        Args:
            path (list): список идентификаторов вершин

        Returns:
            pd.DataFrame: отсортирован по дню, с последовательностями и их частотами
        """
        # получаем словарь соответсвий идентификаторов вершин и данных о записи
        vertex_id = {key: val for (val, key) in self._mapping.items()}
        id_data = []

        # получаем список записей по идентификатам вершин
        for x in path:
            id_data.append(vertex_id[x])
        return pd.DataFrame(id_data, columns=['day', 'seq', 'freq']).sort_values('day')


class Node:
    """
    Класс листа дерева
    """
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
        """
        Парное выравнивание последовательности листа на потенциальных родителей и сохранение скора

        Args:
            potential_parent (str): последовательность потенциального родителя
        """
        aligner = PairwiseAligner()

        # используем локальное выравнивание
        aligner.mode = 'local'

        # заменим дефолтные символы пропуска на другие
        seq = self.seq.replace('-', '')
        potential_parent_undersores = potential_parent.replace('-', '')

        # выравнивание
        score = aligner.score(seq, potential_parent_undersores)

        # сохраняем скор
        self.parent_scores[potential_parent] = score

    def align(self, parents):
        """
        Внешняя функция для выравнивания
        Args:
            parents (pd.DataFrame): последовательности с возможными родителями и их частотами

        Returns:
            str: родитель листа
        """
        self.parent_scores = {}

        for row in parents.itertuples():
            self.pairwise(row.sequence)
        return self.__parent

    @property
    def __parent(self):
        """
        Выбор родителя

        Returns:
            str: последовательность родителя
        """
        # сортируем словарь по скору
        sorted_seqs = sorted(self.parent_scores, key=self.parent_scores.get, reverse=True)

        # если он не пустой
        if sorted_seqs:

            # сохраянем макисмальный скор
            self.score = self.parent_scores[sorted_seqs[0]]

            # возвращаем последовательность родителя
            return sorted_seqs[0]
        print('Cannot find a parent')

    def __repr__(self):
        return f'Node {self.seq}, freq {self.frequency}, has {len(self.children)} children'

