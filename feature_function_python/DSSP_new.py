def getDSSP_final(pdbPATH):
    a = subprocess.check_output([dssp_exec,pdbPATH])
    buf = StringIO.StringIO(a)
    lines = buf.readlines()


    name_idex  = lines.index('  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA \n')

    lines = lines[name_idex+1:]

    data = []
    for line in lines:
        line_data = []
        line_data.append(line[6:10].strip())
        line_data.append(line[11:12].strip())
        line_data.append(line[13:14].strip())
        line_data.append(line[15:17].strip())
        line_data.append(line[18:25])
        line_data.append(line[26:29].strip())
        line_data.append(line[30:34].strip())
        line_data.append(line[35:38].strip())
        data.append(line_data)

    Columns = ['resID',' chain','residue','type','structure','BP1','BP2','ACC']

    DSSP_Data_Frame=pd.DataFrame(data,columns= Columns)


    secondary2score = {'H':0.5,'I':0.5,'G':0.5,'T':1,'S':1,'B':0,'E':0}
    assigneds = secondary2score.keys()
    secondaryScore = DSSP_Data_Frame['type'].apply(lambda x: secondary2score[x] if x in assigneds else 1)
    final_dssp_table = pd.DataFrame(zip(DSSP_Data_Frame['resID'],secondaryScore),columns=['resID','dssp_result'])
    res_len = final_dssp_table.shape[0]

    loop_num = sum (final_dssp_table['dssp_result'] == 1)

    loop_percent = float (loop_num)/float (res_len) * float (100)
    final_dssp_table['loop_percent'] = loop_percent
    final_dssp_table['resID'] = final_dssp_table['resID'].replace('',np.nan, regex=True)

    final_dssp_table = final_dssp_table[~final_dssp_table['resID'].isnull()]
    return final_dssp_table