function [norm_vector_2,nor_Rank] = getGNMeigenvector_rank2(pdbStructure)
[gnmEigVector_2]=getGNM(pdbStructure,2);
norm_vector_2 = abs(gnmEigVector_2);
[~, ~, ranking] = unique(norm_vector_2);
Rank = max(ranking) - ranking  + 1;
nor_Rank = 1-(Rank./length(pdbStructure));