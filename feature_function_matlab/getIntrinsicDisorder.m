function [intrinsicDisorder] =  getIntrinsicDisorder(sequence)
    fastaSeq = sequence2fasta(sequence,'TMP');
    fileName = char(java.util.UUID.randomUUID);
    tmpfile = fopen(fileName,'w');
    try
        fprintf(tmpfile,'%s',fastaSeq);
        fclose(tmpfile);
        [exitCode, data] = system(['RONN ' fileName]);
        if exitCode ~= 0
            error('Useful_func:IntrinsicDisorder','Failed to call RONN.');
        end
        intrinsicDisorder =  sscanf(data(6:end),'%f');
        delete(fileName);
    catch e
        fclose(tmpfile);
        delete(fileName);
        rethrow(e);
    end
end