function [sum_contact] = getAtomicContact_1(pdb_a)
[a_contact] = getContactMatrix(pdb_a,pdb_a,4.5);
a_contact(eye(size(a_contact))~=0)=0;
a_contact_1 = a_contact;

tbl = tabulate([pdb_a.internalResno]);
count_internal = tbl(:,2);
count_internal(count_internal == 0) = '';
n = 0;
m = 1;
for i = 1:length(count_internal)-1
    n = n + count_internal(i)+count_internal(i+1);
    a_contact_1(m:n,m:n) = 0;
    n = n - count_internal(i+1);
    m = m +count_internal(i);
end
%remove_contact = sum(a_contact_1,2);
%pdb_ON = as('name O N',pdb_a);
%minus_atom = [pdb_ON.atomno];
%for i = minus_atom(2:end-1)
    %remove_contact(i) = remove_contact(i)-1;
%end
n = 0;
m = 1;
sum_contact = zeros(1,1);
for i = 1:length(count_internal)
    n = n + count_internal(i);
    res_matrix = a_contact_1(1:end,(m:n));
    sum_row=sum(res_matrix,2);
    sum_row(sum_row==0)= '';
    sum_contact(i,1) = length(sum_row);
    m = m +count_internal(i);
end
