function [gyradius]=getRadiusOfgyration(pdb_a)
tbl = tabulate([pdb_a.internalResno]);
count_internal = tbl(:,2);
count_internal(count_internal == 0) = '';
res_center = zeros(1,1);
n = 0;
m = 1;
for i = 1:length(count_internal)
    n = n + count_internal(i);
    res_cen = getGeometrycenter(pdb_a(m:n));
    res_center(i,1) = res_cen(1);
    res_center(i,2) = res_cen(2);
    res_center(i,3) = res_cen(3);
    m = m +count_internal(i);
end

r_mean = sum(res_center)/length(res_center);
R_sum = 0;
for i = 1:length(res_center)
    R_sum = R_sum+(sum((res_center(i)-r_mean).^2)) ;
end
gyradius = sqrt(R_sum/length(res_center));