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
    replicates = 1
    spectrum_fromtrees = _treesignalc.fromtrees_rescale # default is to use rescaled distances
    #TODO: work with shortened taxa names

    def __init__(self, sp_trees = None, replicates = 1):
        """ if replicates is zero, then we use non-scaled distances; if replicate is one, then we use theoretical upper
        bounds to normalise distances; otherwise we use p-value normalisation
        """
        if (replicates > 1) and (replicates < 4):
            replicates = 4 # too few samples to estimate bounds
        self.replicates = replicates 
        if self.replicates == 0:
            self.spectrum_fromtrees = _treesignalc.fromtrees # default is to use rescaled distances
        if self.replicates > 1:
            self.spectrum_fromtrees = _treesignalc.fromtrees_pvalue # default is to use rescaled distances
        if sp_trees is not None:
            self.update_spstring_from_trees(sp_trees = sp_trees)

    def __call__(self, genetree = None, string = False):
        if genetree is None:
            feat_matrix = []
            if (self.replicates > 1): # p-value needs an extra term
                for gtree in self.sp_trees: # dendropy adds quotes whenever there's underscores
                    gt_string = gtree.as_string(schema="newick",suppress_edge_lengths=True).rstrip().replace("'","")
                    feat_matrix.append(self.spectrum_fromtrees(gt_string, self.sp_string, self.replicates))
            else: # unscaled or divided by theoretical maximum
                for gtree in self.sp_trees: # dendropy adds quotes whenever there's underscores
                    gt_string = gtree.as_string(schema="newick",suppress_edge_lengths=True).rstrip().replace("'","")
                    feat_matrix.append(self.spectrum_fromtrees(gt_string, self.sp_string))
            return numpy.array(feat_matrix)
        else:
            if not string: # dendropy Tree object
                genetree = genetree.as_string(schema="newick",suppress_edge_lengths=True).rstrip().replace("'", "")
            if (self.replicates > 1): # p-value needs an extra term
                return numpy.array(self.spectrum_fromtrees(genetree, self.sp_string, self.replicates))
            else:
                return numpy.array(self.spectrum_fromtrees(genetree, self.sp_string))

    def update_sptrees_from_string(self, sp_string = None):
        self.sp_string = sp_string
        self.sp_trees = dendropy.TreeList.get( data=self.sp_string.replace("'", ""), schema="newick", preserve_underscores=True)

    def update_spstring_from_trees(self, sp_trees = None):
        self.sp_trees = dendropy.TreeList(sp_trees)
        self.sp_string = self.sp_trees[0].as_string("newick",suppress_edge_lengths=True).rstrip().replace("'", "")
        for i in range (1,len(self.sp_trees)):
            self.sp_string += self.sp_trees[i].as_string("newick",suppress_edge_lengths=True).rstrip().replace("'", "")

    def generate_sptrees(self, n_species = 16, chain_size = 100, n_spr = 1):
        if n_species < 4:
            n_species = 4
        if chain_size < 2:
            chain_size = 2
        self.sp_string = _treesignalc.generate_spr_trees(n_species, chain_size, n_spr)
        #self.sp_trees = dendropy.TreeList.get( data=self.sp_string, schema="newick", preserve_underscores=True, quote_underscores=False)
        self.sp_trees = dendropy.TreeList.get( data=self.sp_string, schema="newick", preserve_underscores=True)
