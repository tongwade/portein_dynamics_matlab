function [inter_matrix] = getCharge_res(pdb_a)
pdb_ca = as('name CA',pdb_a);
pdb_charge = as('resn ASP GLU HIS LYS ARG',pdb_ca);
inter = [pdb_charge.internalResno];
inter_matrix = zeros(length(pdb_ca),1);
inter_matrix(inter) = 1;