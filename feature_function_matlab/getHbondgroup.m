function [h_bond_group] = getHbondgroup(pdbname,index,output_path,input_path)
%index is residue id
cd(output_path)
try
    exitCode = system(['hbplus ' input_path lower(pdbname(1:4)) '.pdb cc.new -d 3.5 -A 120 0 0 -e HOH " O " 0 -O']);
    if exitCode ~= 0
            error('Useful_func:Hbondgroup','Failed to call hbplus.');
    end
catch e
    rethrow(e);
end
try
    texttt = textread([output_path '/' lower(pdbname(1:4)) '.hb2'],'%s','delimiter','\n','whitespace','');

    textData1 = regexp(texttt,'\s+','split'); 
    data = textData1(9:end);
    cell_matrix = cell(length(data),2);
    for i = 1:length(data)
        cell_matrix{i,1} = data{i}{1};
        cell_matrix{i,2} = data{i}{3};
        %cell_matrix{i,1:13} = data{i}
    end
    if length(data{1}{2}) > 6
        hbond_matrix = zeros(1,1);

        for i = 1:length(data)
            hbond_matrix(i,1) = str2double(data{i}{1}(2:5));
            hbond_matrix(i,2) = str2double(data{i}{2}(2:5));
        end
    else
        hbond_matrix = zeros(1,1);

        for i = 1:length(data)
            hbond_matrix(i,1) = str2double(data{i}{1}(2:5));
            hbond_matrix(i,2) = str2double(data{i}{3}(2:5));
        end
    end
    for i = 1:length(hbond_matrix)
        if hbond_matrix(i,1) == hbond_matrix(i,2) | hbond_matrix(i,1)-1 == hbond_matrix(i,2) | hbond_matrix(i,1)+1 == hbond_matrix(i,2)
            hbond_matrix(i,1) =0; 
            hbond_matrix(i,2) =0; 
        end
    end
    all_resid = [hbond_matrix(:,1) ; hbond_matrix(:,2)];
    all_resid = all_resid(all_resid~=0);
    tbl = tabulate(all_resid);

    resid = tbl(:,1);
    h_bond_count = tbl(:,2);

    resid(h_bond_count == 0) = '';
    h_bond_count(h_bond_count == 0) = '';
    
    h_bond_group = zeros(length(index),1);
    for i = 1:length(h_bond_count)
        id = resid(i);
        y = find(index == id);
        h_bond_group(y,1) = h_bond_count(i);
    end
    delete([lower(pdbname(1:4)) '.hb2']);
    delete([lower(pdbname(1:4)) '.h']);
    

catch e
    rethrow(e);
end
