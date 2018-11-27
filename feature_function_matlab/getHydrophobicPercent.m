function [hydrophobic_precent] = getHydrophobicPercent(pdb_a)
pdb_ca = as('name CA',pdb_a);

pdb_hydorphobic_matrix = as('resn THR GLY SER HIS GLN ARG LYS ASN GLU ASP',pdb_ca,1)';
hydrophobic_precent = (sum(pdb_hydorphobic_matrix)/length(pdb_ca))*100;