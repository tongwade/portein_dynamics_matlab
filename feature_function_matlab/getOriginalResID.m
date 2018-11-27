function [ori_resid] = getOriginalResID(PDBID,resid)
load(['D:\ML_project\test_data_muscle\new_resid\' PDBID '_new_resid.mat']);
[~,~,all_ori_resid] = xlsread(['D:\ML_project\test_data_muscle\ori_resid\' PDBID '_ori_resid.csv']);
[~,index] = ismember(resid,new_resid,'R2012a');
ori_resid = arrayfun(@(x) all_ori_resid(x),index);
end
