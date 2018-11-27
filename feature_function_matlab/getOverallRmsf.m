function [rmsf_overall] = getOverallRmsf(pdb_ca,GNMValue)
gnm_sig = 0;
for i = 1:length(pdb_ca)-1
    gnm_sig = gnm_sig + (1/GNMValue(i));
end
rmsf_overall = sqrt(gnm_sig/length(pdb_ca));
    