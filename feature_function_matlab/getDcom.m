function [stand_dist,res_dist,res_center] = getDcom(pdb_a)
pdb_a = assignMass(pdb_a);
protein_center = getCenterOfMass(pdb_a);
tbl = tabulate([pdb_a.internalResno]);
count_internal = tbl(:,2);
count_internal(count_internal == 0) = '';
res_center = zeros(1,1);
n = 0;
m = 1;
for i = 1:length(count_internal)
    n = n + count_internal(i);
    res_cen = getCenterOfMass(pdb_a(m:n));
    res_center(i,1) = res_cen(1);
    res_center(i,2) = res_cen(2);
    res_center(i,3) = res_cen(3);
    m = m +count_internal(i);
end
res_dist = zeros(1,1);
for i = 1:length(count_internal)
    X = [protein_center;res_center(i,:)];
    dist = pdist(X);
    res_dist(i,1) = dist;
end
stand_dist = (res_dist - mean(res_dist))./std(res_dist);
