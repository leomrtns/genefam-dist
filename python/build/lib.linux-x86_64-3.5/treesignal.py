import __treesignalc, dendropy

gtre_string   = dendropy.simulate.treesim.birth_death_tree(birth_rate=1.0, death_rate=0.5, ntax=8).as_string("newick")

def create_sptre_string ():
    sptre_string  = dendropy.simulate.treesim.birth_death_tree(birth_rate=1.0, death_rate=0.5, ntax=8).as_string("newick")
    sptre_string += dendropy.simulate.treesim.birth_death_tree(birth_rate=1.0, death_rate=0.5, ntax=8).as_string("newick")
    sptre_string += dendropy.simulate.treesim.birth_death_tree(birth_rate=1.0, death_rate=0.5, ntax=8).as_string("newick")
    sptre_string += dendropy.simulate.treesim.birth_death_tree(birth_rate=1.0, death_rate=0.5, ntax=8).as_string("newick")
    return sptre_string

#t1.find_node_with_taxon_label("T6").taxon.label="T1" # enough to change the newick format: doesnt change Taxon, only its label 
