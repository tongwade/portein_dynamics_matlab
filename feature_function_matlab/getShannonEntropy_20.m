function [shannon_entropy20] = getShannonEntropy_20(GNMValues)
top_20 = round(length(GNMValues)*0.2);
top_20_GNMValues = GNMValues(1:top_20);
binrange = linspace(top_20_GNMValues(1),top_20_GNMValues(end),length(top_20_GNMValues)+1);
bincounts = histc(top_20_GNMValues,[binrange]);
%要把最後一個值加到最後一個bin 裡面
bincounts = bincounts(1:end-1);
bincounts(end) = bincounts(end)+1;
bincounts(end+1) = 0;
bincounts = bincounts(bincounts~=0);
len_ca = sum(bincounts);
shannon_entropy20 = ((-1)*sum((bincounts./len_ca).*log2(bincounts./len_ca))) / log2((length(top_20_GNMValues)));
%shannon entropy normalize =divide log2(ca length)

