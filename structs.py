class GeneNode:
    def __init__(self, g, e=None):
        '''
        gene (string): name of gene
        edges (list): list of Edges that this gene regulates
        '''
        self.gene = g
        self.edges = e if e is not None else []
    
    def add_edge(self, t, a):
        '''
        Add edge if it doesn't already exist.

        Parameters:
        t (GeneNode): target of Edge
        a (boolean): True if Edge is activating, False if Edge is inhibiting
        '''
        if self.get_edge(t, a) is not None:
            return
        edge = Edge(t, a)
        self.edges.append(edge)

    def get_edge(self, t, a):
        '''
        Return the Edge containing a given target with a given type of regulation.

        Parameters: 
        t (GeneNode): target of desired Edge
        a (boolean): True if Edge is activating, False if Edge is inhibiting

        Returns: 
        edge (Edge): Edge containing t and proper regualtion type
        '''
        for edge in self.edges:
            if edge.target == t:
                if edge.act == a:
                    return edge
                return None
        return None

class Edge:
    def __init__(self, t, a):
        '''
        target (GeneNode): the gene the Edge is pointing to/regulating
        act (boolean): True if Edge is activating, False if Edge is inhibiting
        '''
        self.target = t
        self.act = a
