function [feature] = getGNMeigenvalue_vector2(pdb_a)
pdb_ca = as('name CA',pdb_a);
[pdbStructure,GNMValue]=GNM(pdb_ca);
[gnmEigVector]=getGNM(pdbStructure,2);
eig_2 = 1/GNMValue(2);
feature = gnmEigVector*eig_2;
