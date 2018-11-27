function [polarity] = getPolarity(pdbStructureOrSequence)
order = 'LPMWAVFIGSTCNQYHDEKR';
% res_name = {'LEU','PRO','MET','TRP','ALA','VAL','PHE','ILE','GLY','SER','THR','CYS','ASN','GLN','TYR','HIS','ASP','GLU','LYS','ARG'};
polaritys = [4.9 8 5.7 5.4 8.1 5.9 5.2 5.2 9 9.2 8.6 5.5 11.6 10.5 6.2 10.4 13 12.3 11.3 10.5];

if isstruct(pdbStructureOrSequence)
    allAtomNames = {pdbStructureOrSequence.atmname};
    noh = pdbStructureOrSequence(cellfun(@isempty,regexp(allAtomNames,'H*')));
    chainSeq = getSequence(noh);
    seq = [];
    for i = 1: length(chainSeq)
        seq = [seq chainSeq{i}];
    end
elseif ischar(pdbStructureOrSequence)
    seq = pdbStructureOrSequence;
else
    error('UsefulFunc:Polarity','Input should PDB structure or sequence (char)')
end

[~,index] = ismember(seq,order);
polarity = polaritys(index);