function [sum_contact] = getAtomicContact_5(pdb_a)
[a_contact] = getContactMatrix(pdb_a,pdb_a,4.5);
a_contact(eye(size(a_contact))~=0)=0;
a_contact_1 = a_contact;

tbl = tabulate([pdb_a.internalResno]);
count_internal = tbl(:,2);
count_internal(count_internal == 0) = '';

sum_contact = zeros(1,1);
res1_count = count_internal(1)+count_internal(2)+count_internal(3);
res_matrix = a_contact_1(1:end,1:res1_count);
%sum_row = arrayfun(@(x) sum(res_matrix(x,:)),(1:size(res_matrix,1)).');
sum_row=sum(res_matrix,2);
sum_row(sum_row==0)= '';
res1_minus = count_internal(1)+count_internal(2)+count_internal(3)+count_internal(4);
minus_matrix = a_contact_1(1:res1_minus,1:res1_count);
%sum_minus_row = arrayfun(@(x) sum(minus_matrix(x,:)),(1:size(minus_matrix,1)).');
sum_minus_row=sum(minus_matrix,2);
sum_minus_row(sum_minus_row==0)= '';
sum_contact(1,1) = length(sum_row)-length(sum_minus_row);

res2_count = count_internal(1)+count_internal(2)+count_internal(3)+count_internal(4);
res_matrix = a_contact_1(1:end,1:res2_count);
%sum_row = arrayfun(@(x) sum(res_matrix(x,:)),(1:size(res_matrix,1)).');
sum_row=sum(res_matrix,2);
sum_row(sum_row==0)= '';
res2_minus = count_internal(1)+count_internal(2)+count_internal(3)+count_internal(4)+count_internal(5) ;
minus_matrix = a_contact_1(1:res2_minus,1:res2_count);
%sum_minus_row = arrayfun(@(x) sum(minus_matrix(x,:)),(1:size(minus_matrix,1)).');
sum_minus_row=sum(minus_matrix,2);
sum_minus_row(sum_minus_row==0)= '';
sum_contact(2,1) = length(sum_row)-length(sum_minus_row);

res2_count = count_internal(1)+count_internal(2)+count_internal(3)+count_internal(4)+count_internal(5);
res_matrix = a_contact_1(1:end,1:res2_count);
%sum_row = arrayfun(@(x) sum(res_matrix(x,:)),(1:size(res_matrix,1)).');
sum_row=sum(res_matrix,2);
sum_row(sum_row==0)= '';
res2_minus = count_internal(1)+count_internal(2)+count_internal(3)+count_internal(4)+count_internal(5)+count_internal(6) ;
minus_matrix = a_contact_1(1:res2_minus,1:res2_count);
%sum_minus_row = arrayfun(@(x) sum(minus_matrix(x,:)),(1:size(minus_matrix,1)).');
sum_minus_row=sum(minus_matrix,2);
sum_minus_row(sum_minus_row==0)= '';
sum_contact(3,1) = length(sum_row)-length(sum_minus_row);

n = 1;
m = 1;
res_count = count_internal(1);
res_minus = count_internal(1);
for i = 1:length(count_internal)-6
    n = n +count_internal(i);
    res_count = res_count + count_internal(i+1)+count_internal(i+2)+count_internal(i+3)+count_internal(i+4)+count_internal(i+5);
    res_matrix = a_contact_1(1:end,(n:res_count));
    %sum_row = arrayfun(@(x) sum(res_matrix(x,:)),(1:size(res_matrix,1)).');
    sum_row=sum(res_matrix,2);
    sum_row(sum_row==0)= '';
    
    res_minus = res_minus +count_internal(i+1)+count_internal(i+2)+count_internal(i+3)+count_internal(i+4)+count_internal(i+5)+count_internal(i+6);
    minus_matrix = a_contact_1(m:res_minus,n:res_count);
    %sum_minus_row = arrayfun(@(x) sum(minus_matrix(x,:)),(1:size(minus_matrix,1)).');
    sum_minus_row=sum(minus_matrix,2);
    sum_minus_row(sum_minus_row==0)= '';
    
    sum_contact(i+3,1) = length(sum_row)-length(sum_minus_row);
    
    res_count = res_count -(count_internal(i+2)+count_internal(i+3)+count_internal(i+4)+count_internal(i+5));
    res_minus = res_minus -(count_internal(i+2)+count_internal(i+3)+count_internal(i+4)+count_internal(i+5)+count_internal(i+6));
    m = m +count_internal(i);
end

n = n +count_internal(i+1);
res_count = res_count + count_internal(i+2)+count_internal(i+3)+count_internal(i+4)+count_internal(i+5);
res_matrix = a_contact_1(1:end,(n:res_count));
%sum_row = arrayfun(@(x) sum(res_matrix(x,:)),(1:size(res_matrix,1)).');
sum_row=sum(res_matrix,2);
sum_row(sum_row==0)= '';

res_minus = res_minus+count_internal(i+2)+count_internal(i+3)+count_internal(i+4)+count_internal(i+5)+count_internal(i+6);
minus_matrix = a_contact_1(m:res_minus,n:res_count);
%sum_minus_row = arrayfun(@(x) sum(minus_matrix(x,:)),(1:size(minus_matrix,1)).');
sum_minus_row=sum(minus_matrix,2);
sum_minus_row(sum_minus_row==0)= '';

sum_contact(i+4,1) = length(sum_row)-length(sum_minus_row);

n = n +count_internal(i+2);
res_count = res_count ;
res_matrix = a_contact_1(1:end,(n:end));
%sum_row = arrayfun(@(x) sum(res_matrix(x,:)),(1:size(res_matrix,1)).');
sum_row=sum(res_matrix,2);
sum_row(sum_row==0)= '';

m = m +count_internal(i+1);
res_minus = res_minus;
minus_matrix = a_contact_1(m:end,n:end);
%sum_minus_row = arrayfun(@(x) sum(minus_matrix(x,:)),(1:size(minus_matrix,1)).');
sum_minus_row=sum(minus_matrix,2);
sum_minus_row(sum_minus_row==0)= '';

sum_contact(i+5,1) = length(sum_row)-length(sum_minus_row);

n = n +count_internal(i+3);
res_count = res_count ;
res_matrix = a_contact_1(1:end,(n:end));
%sum_row = arrayfun(@(x) sum(res_matrix(x,:)),(1:size(res_matrix,1)).');
sum_row=sum(res_matrix,2);
sum_row(sum_row==0)= '';

m = m +count_internal(i+2);
res_minus = res_minus;
minus_matrix = a_contact_1(m:end,n:end);
%sum_minus_row = arrayfun(@(x) sum(minus_matrix(x,:)),(1:size(minus_matrix,1)).');
sum_minus_row=sum(minus_matrix,2);
sum_minus_row(sum_minus_row==0)= '';

sum_contact(i+6,1) = length(sum_row)-length(sum_minus_row);
