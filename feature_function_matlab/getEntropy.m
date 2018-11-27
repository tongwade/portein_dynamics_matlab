function [entropy] = getEntropy(GNMValues,pdb_ca)
log_GNMValue = log(GNMValues);
eig_muti = sum(log_GNMValue);
N = length(pdb_ca);
entropy = (0.5*N)+(((N*log(2*pi))-eig_muti)*0.5);