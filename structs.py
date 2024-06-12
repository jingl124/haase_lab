class GeneNode:
    def __init__(self, g, e=None):
        '''
        gene (string): name of gene
        edges (list): list of Edges that this gene regulates
        '''
        self.gene = g
        self.edges = e if e is not None else []
    
    def add_edge(self, t, a):
        edge = Edge(t, a)
        # # if we have a global edge tracker in grn_finder.py, then we don't need this for loop
        # for e in self.edges:
        #     if e.target == t and e.act == a:
        #         return
        self.edges.append(edge)

class Edge:
    def __init__(self, t, a):
        '''
        target (GeneNode): the gene the Edge is pointing to/regulating
        act (boolean): True if Edge is activating, False if Edge is inhibiting
        '''
        self.target = t
        self.act = a
