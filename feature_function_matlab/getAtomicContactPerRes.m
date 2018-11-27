function [groupTotalContact]=getAtomicContactPerRes(pdb,cutoff,groupNeighbor,skipNeighbor)
    numAtom = length(pdb);
    atomNumPerRes = getAtomNumPerRes(pdb);
    numRes = length(atomNumPerRes);
    contactM = getContactMatrix(pdb,pdb,cutoff);
    index = zeros(numRes + 1,1 );
    index(2:end) = cumsum(atomNumPerRes);
    
    resNumPerGroup = zeros(1,numRes);
    groupByGroup = zeros(numAtom,numRes);
    for i = 1:numRes
       startResIndex = i - groupNeighbor;
       endResIndex = i + groupNeighbor;
       if startResIndex <1
           startResIndex = 1;
       end
       if endResIndex > numRes
           endResIndex = numRes;
       end
       resNumPerGroup(i) = endResIndex - startResIndex + 1;
       groupByGroup(:,i) = sum(contactM(:,index(startResIndex)+1:index(endResIndex+1)),2);
    end
    groupByGroup = groupByGroup > 0;
    for i = 1:numRes
        startResIndex = i - groupNeighbor - skipNeighbor; 
        endResIndex = i + groupNeighbor + skipNeighbor;
       if startResIndex < 1
           startResIndex = 1;
       end
       if endResIndex > numRes
           endResIndex = numRes;
       end
       groupByGroup(index(startResIndex)+1:index(endResIndex+1),i) = 0;
    end
    groupTotalContact = sum(groupByGroup,1)./resNumPerGroup;
end