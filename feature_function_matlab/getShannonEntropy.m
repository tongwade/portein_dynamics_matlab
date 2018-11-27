function [shannon_entropy] = getShannonEntropy(pdb_ca,GNMValue)
binrange = linspace(GNMValue(1),GNMValue(end),length(pdb_ca));
bincounts = histc(GNMValue,[binrange]);
%要把最後一個值加到最後一個bin 裡面
bincounts = bincounts(1:end-1);
bincounts(end) = bincounts(end)+1;
bincounts(end+1) = 0;
bincounts = bincounts(bincounts~=0);
len_ca = sum(bincounts);
shannon_entropy = ((-1)*sum((bincounts./len_ca).*log2(bincounts./len_ca))) / log2((length(GNMValue)));
%shannon entropy normalize =divide log2(ca length)


