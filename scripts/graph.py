from igraph import Graph as igraph_graph


class Graph(igraph_graph):

    def set_edge_weights(self, edges):
        for edge in edges:
            self.add_edge(edge[0], edge[1], weight=edges[edge])

    def find_longest_path(self, tree):
        """ Finding the longest path in a graph

        Args:
            tree (tree): patient region tree
        Return:
            max_weighted_lenght (list): vertices of longest path
            longest_path (float): weight of longest path """
        paths = self.get_shortest_paths(0, weights='weight', mode='ALL')
        longest_path = None
        longest_path_id = -1
        max_weighted_lenght = -1

        for path_id, path in enumerate(paths):

            if path:
                first = 0
                path_weight = 0
                for second in range(1, len(path)):
                    first_el, second_el = path[first], path[second]
                    path_weight += tree.graph[(first_el, second_el)]
                    first = second
                if max_weighted_lenght < path_weight:
                    max_weighted_lenght = path_weight
                    longest_path = path
                    longest_path_id = path_id

        max_edges = self.get_shortest_paths(0, weights='weight', mode='ALL', output="epath")[longest_path_id]
        return max_weighted_lenght, longest_path, max_edges

    def colour_path(self, edges):
        """ colouring the edges of the longest path """
        for edge_id in edges:
            self.es[edge_id]['color'] = 'red'
