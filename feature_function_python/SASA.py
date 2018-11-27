import os
import pandas as pd
import re
import shutil
import uuid
import subprocess

RSAColumns = ['Header','ResName','chainID','resID','All-atoms-ABS','All-atoms-REL','Total-Side-ABS','Total-Side-REL','Main-Chain-ABS','Main-Chain-REL','Non-polar-ABS','Non-polar-REL','All polar-ABS','All polar-REL']
respattern = re.compile('^RES')
naccessExec ='/opt/Naccess/naccess'
def getSASA(pdbPath):  
    # generate tmpfile name
    tmpName = str(uuid.uuid1());
    tmpFullName = tmpName + '.pdb';
    rsa_file= tmpName + '.rsa'
    
    # create tmpfile
    shutil.copy(pdbPath,os.path.join(os.curdir,tmpFullName)) # create tmpfile
    subprocess.call([naccessExec,tmpFullName])
    with open(rsa_file, 'r') as naccessData:      
        lines = naccessData.readlines()
        data = []
        for line in lines:
            if respattern.search(line):
                subdata = []
                subdata.append(line[0:3].strip()) # add header
                subdata.append(line[4:7].strip()) #add resname
                subdata.append(line[7:9].strip()) #add chian
                subdata.append(line[9:13].strip()) #add resid
                subdata += line[13:].strip().split() #add sasa data
                data.append(subdata)
            
        SASA_Data_Frame=pd.DataFrame(data,columns= RSAColumns)
        SASA_Data_Frame.drop('Header',axis=1,inplace=True)
    
    # remove tmporary files
    for fileName in os.listdir(os.curdir):
        if re.search('^%s'%tmpName,fileName):
            os.remove(fileName) 
    return SASA_Data_Frame                   