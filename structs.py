class GeneNode:
    def __init__(self, g, e=None):
        '''
        gene (string): name of gene
        edges (list): list of Edges that this gene regulates
        '''
        self.gene = g
        self.edges = e if e is not None else []
    
    def add_edge(self, target, act):
        '''
        Add edge if it doesn't already exist.

        Parameters:
        target (GeneNode): target of Edge
        act (boolean): True if Edge is activating, False if Edge is inhibiting
        '''
        if self.get_edge(target, act) is not None:
            return
        edge = Edge(target, act)
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
    def __init__(self, target, act, time):
        '''
        target (GeneNode): the gene the Edge is pointing to/regulating
        act (boolean): True if Edge is activating, False if Edge is inhibiting
        time (float): the time of activation/inhibition
        '''
        self.target = target
        self.act = act
        self.time = time
