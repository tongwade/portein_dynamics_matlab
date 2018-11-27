function [pdb_phobic_matrix] = getHydrophobic_matrix(pdb_a)
pdb_ca = as('name CA',pdb_a);
pdb_phobic_matrix = as('resn ALA ILE LEU PHE VAL PRO GLY',pdb_ca,1)';