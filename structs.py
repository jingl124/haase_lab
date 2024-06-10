import graphviz

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

def visualize_gene_network(gene_nodes):
    dot = graphviz.Digraph(comment='Gene Regulatory Network')
    
    # Add nodes
    for gene in gene_nodes:
        dot.node(gene.gene)
    
    # Add edges
    for gene in gene_nodes:
        for edge in gene.edges:
            if edge.act:
                color = 'green'
                arrowhead = 'normal'
            else:
                color = 'red'
                arrowhead = 'tee'
            dot.edge(gene.gene, edge.target.gene, color=color, arrowhead = arrowhead)
    
    # Render the graph
    dot.render('gene_network', view=True)

# Example usage
gene_a = GeneNode('GeneA')
gene_b = GeneNode('GeneB')
gene_c = GeneNode('GeneC')

gene_a.add_edge(gene_b, True)  # GeneA activates GeneB
gene_b.add_edge(gene_c, False) # GeneB inhibits GeneC

gene_nodes = [gene_a, gene_b, gene_c]

visualize_gene_network(gene_nodes)