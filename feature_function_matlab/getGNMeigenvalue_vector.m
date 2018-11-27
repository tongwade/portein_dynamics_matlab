function [gnmEigVector_eig_1,gnmEigVector_eig_2] = getGNMeigenvalue_vector(pdbStructure,GNMValues)
[gnmEigVector_1]=getGNM(pdbStructure,1);
eig_1 = 1/GNMValues(1);
gnmEigVector_eig_1 = abs(gnmEigVector_1)*eig_1;
[gnmEigVector_2]=getGNM(pdbStructure,2);
eig_2 = 1/GNMValues(2);
gnmEigVector_eig_2 = abs(gnmEigVector_2)*eig_2;
