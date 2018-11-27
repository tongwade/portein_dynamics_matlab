function [sideChainLength] = getSideChainLength(pdbStructureOrSequence)
order = 'GAVLIMFWPSTCYNQDEKRH';
numSideChainHeavyAtom = [0,1,3,4,4,4,7,10,3,2,3,2,8,4,5,4,5,5,7,6];

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
    error('UsefulFunc:SideChainLength','Input should PDB structure or sequence (char)')
end

[~,index] = ismember(seq,order);
sideChainLength = numSideChainHeavyAtom(index);

