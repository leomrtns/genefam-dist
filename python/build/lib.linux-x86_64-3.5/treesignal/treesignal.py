import _treesignalc, dendropy, numpy

class TreeSignal(object):
    """
    This basic class contains generates a common coordinate system for gene families based on reference species trees.

    Its main parts are the __init__ (for initial setup) and __call__ (for generation of feature matrix), and does not
    need to store the gene tree info, only the reference (a.k.a. anchoring) species trees.

    -- Currently this class only works to test the algos in development
    """
    sp_trees = None
    sp_string = None
    spectrum_fromtrees = _treesignalc.fromtrees # default is to use original distances
    #TODO: work with shortened taxa names

    def __init__(self, sp_trees = None, n_species = 16, chain_size = 100, n_spr = 1, pvalue = False):
        if pvalue: # instead of original distances, use p-value scores as distances 
            self.spectrum_fromtrees = _treesignalc.fromtrees_pvalue
        if sp_trees is None:
            if n_species < 4:
                n_species = 4
            if chain_size < 2:
                chain_size = 2
            self.sp_string = _treesignalc.generate_spr_trees(n_species, chain_size, n_spr)
            self.sp_trees = dendropy.TreeList.get( data=self.sp_string, schema="newick")
        else:
            self.update_spstring_from_trees(sp_trees = sp_trees)

    def __call__(self, genetree = None):
        if genetree is None:
            feat_matrix = []
            for gtree in self.sp_trees:
                feat_matrix.append(self.spectrum_fromtrees(gtree.as_string(schema="newick"), self.sp_string))
            return numpy.array(feat_matrix)
        else:
            return numpy.array(self.spectrum_fromtrees(genetree.as_string(schema="newick"), self.sp_string))

    def update_sptrees_from_string(self, sp_string = None):
        self.sp_string = sp_string
        self.sp_trees = dendropy.TreeList.get( data=self.sp_string, schema="newick")
    def update_spstring_from_trees(self, sp_trees = None):
        self.sp_trees = dendropy.TreeList(sp_trees)
        self.sp_string = self.sp_trees[0].as_string("newick").rstrip()
        for i in range (1,len(self.sp_trees)):
            self.sp_string += self.sp_trees[i].as_string("newick").rstrip()
