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
        self.edges.append(edge)

class Edge:
    def __init__(self, t, a):
        '''
        target (GeneNode): the gene the Edge is pointing to/regulating
        act (int): 1 if Edge is activating, 0 if Edge is inhibiting
        '''
        self.target = t
        self.act = a