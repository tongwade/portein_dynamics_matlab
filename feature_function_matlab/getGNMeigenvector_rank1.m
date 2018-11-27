function [norm_vector_1,nor_Rank] = getGNMeigenvector_rank1(pdbStructure)
[gnmEigVector_1]=getGNM(pdbStructure,1);
norm_vector_1 = abs(gnmEigVector_1);
[~, ~, ranking] = unique(norm_vector_1);
Rank = max(ranking) - ranking  + 1;
nor_Rank = 1-(Rank./length(pdbStructure));



