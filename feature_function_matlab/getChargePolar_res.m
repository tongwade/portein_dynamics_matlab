function [pdb_charge_matrix,pdb_polar_matrix] = getChargePolar_res(pdb_a)
pdb_ca = as('name CA',pdb_a);
pdb_charge_matrix = as('resn ASP GLU HIS LYS ARG',pdb_ca,1)';
pdb_polar_matrix = as('resn TYR HIS ASP ASN GLU GLN LYS ATG SER THR CYS',pdb_ca,1)';