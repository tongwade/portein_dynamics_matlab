function [ssbond_feature,protein_hole_ss_bond] = getssbondMatrix(pdb_a)
pdb_CYS = as('resn CYS',pdb_a);

pdb_ca = as('name CA',pdb_a);

pdb_ss = as('name SG',pdb_CYS);
[ss_contact] = getContactMatrix(pdb_ss,pdb_ss,2.4);
ss_contact(eye(size(ss_contact))~=0)=0;

%real_contact = ss_contact & angle_matrix;
real_contact = ss_contact;
%找到距離2.4 角度大於90度、小於270度的contact
%acosd(dot(vector(1,:),vector(6,:))/norm(vector(1,:))/norm(vector(6,:)))

ss_interno = [pdb_ss.internalResno];
close_matrix = abs(bsxfun(@minus,ss_interno',ss_interno));
close_matrix(close_matrix == 1) = 0;
close_matrix(close_matrix >1) = 1;
%find neighbor atom
no_close_matrix = real_contact & close_matrix;
%remove neighbor atom in real_contact matrix

ssDistance = getPairwiseDistance(pdb_ss,pdb_ss);
ssDistance(ssDistance >2.4) = 0;
real_ssDistance = ssDistance.*no_close_matrix;
real_ssDistance(real_ssDistance == 0) = 3;
%把matrix0值改成3 方便取最小值(因為ssbond小於2.4)
remove_double_contact = zeros(length(real_ssDistance),length(real_ssDistance));
for i = 1:length(real_ssDistance)
    one_matrix = real_ssDistance(i,1:length(real_ssDistance));
    close_num = min(one_matrix);
    if close_num == 3
        index = one_matrix == close_num;
        remove_double_contact(i,index) = 0;
    else
        index = find(one_matrix == close_num);
        remove_double_contact(i,index) = close_num;
    end
end
result_contact_matrix = remove_double_contact & real_contact;
contact_index = find(sum(result_contact_matrix) == 1);
contact_internalResno = [pdb_ss(contact_index).internalResno];
%找到contact residue的internalresno
ca_internal_matrix = [pdb_ca.internalResno];
ssbond_feature = zeros(length(pdb_ca),1);
for i = 1:length(contact_index)
    contact_resno = contact_internalResno(i);
    contact_inca_index = find(ca_internal_matrix == contact_resno);
    ssbond_feature(contact_inca_index,1) = 1;
end
count_ssbond = sum(ssbond_feature);
protein_hole_ss_bond = count_ssbond/2;