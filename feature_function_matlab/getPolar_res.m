function [inter_matrix] = getPolar_res(pdb_a)
pdb_ca = as('name CA',pdb_a);
pdb_polar = as('resn TYR HIS ASP ASN GLU GLN LYS ATG SER THR CYS',pdb_ca);
inter = [pdb_polar.internalResno];
inter_matrix = zeros(length(pdb_ca),1)
inter_matrix(inter) = 1;