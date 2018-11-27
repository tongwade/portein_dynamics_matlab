function [feature] = getGNMeigenvalue_vector1(pdb_a)
pdb_ca = as('name CA',pdb_a);
[pdbStructure,GNMValue]=GNM(pdb_ca);
[gnmEigVector]=getGNM(pdbStructure,1);
eig_1 = 1/GNMValue(1);
feature = gnmEigVector*eig_1;