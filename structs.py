class GeneNode:
    def __init__(self, g, e=None):
        '''
        gene (string): name of gene
        edges (list): list of Edges that this gene regulates
        '''
        self.gene = g
        self.edges = e if e is not None else []
    
    def add_edge(self, target, act, time=None):
        '''
        Add edge if it doesn't already exist.

        Parameters:
        target (GeneNode): target of Edge
        act (boolean): True if Edge is activating, False if Edge is inhibiting
        time (float): time of activation/inhibition, None if not applicable
        '''
        if self.get_edge(target, act) is not None:
            return
        edge = Edge(target, act, time)
        self.edges.append(edge)

    def get_edge(self, t, a=None):
        '''
        Return the Edge containing a given target with a given type of regulation.

        Parameters: 
        t (GeneNode): target of desired Edge
        a (boolean): True if Edge is activating, False if Edge is inhibiting. 
            None if type of regulation doesn't matter.

        Returns: 
        edge (Edge): Edge containing t and proper regualtion type
        '''
        for edge in self.edges:
            if edge.target == t:
                if a == None or edge.act == a:
                    return edge
                return None
        return None
    
    def remove_edge(self, t, a=None):
        '''
        Remove the Edge of a GeneNode given the edge's target.

        Parameters:
        t (GeneNode): 
        a (boolean): True if Edge is activating, False if Edge is inhibiting. 
            None if type of regulation doesn't matter. 
        '''
        e = self.get_edge(t, a)
        if e == None:
            return
        for edge in self.edges:
            if edge.target == t and (a == None or edge.act == a):
                self.edges.remove(edge)
                return

class Edge:
    def __init__(self, target, act, time):
        '''
        target (GeneNode or GroupNode): the gene the Edge is pointing to/regulating
        act (boolean): True if Edge is activating, False if Edge is inhibiting
        time (float): the time of activation/inhibition
        '''
        self.target = target
        self.act = act
        self.time = time

class GroupNode:
    def __init__(self, nodes, type=None):
        '''
        nodes (list of GeneNodes or GroupNodes): elements of the ortholog group or complex
        type (string): "ortholog" or "complex"
        '''
        self.nodes = nodes
        self.type = type

    def get_nodes(self):
        return self.nodes  

