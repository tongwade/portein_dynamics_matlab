function [pca_3_1,pca_2_1,S] = getPCAfeature(pdb_ca)
ca_center = getGeometrycenter(pdb_ca);
res_coord = getCoordfromca(pdb_ca);
Q_matrix = zeros(3,length(pdb_ca));
Q_matrix(1,1:length(pdb_ca)) = res_coord(:,1)' - ca_center(1);
Q_matrix(2,1:length(pdb_ca)) = res_coord(:,2)' - ca_center(2);
Q_matrix(3,1:length(pdb_ca)) = res_coord(:,3)' - ca_center(3);
C = Q_matrix*transpose(Q_matrix);
[U,S] = eig(C);
S=diag(S);
%[u_1,s_1]=getPrincipalAxis(pdb_ca)
pca_3_1 = sqrt(min(S)/max(S));
S = sort(S);
pca_2_1 = sqrt(S(2)/max(S));


