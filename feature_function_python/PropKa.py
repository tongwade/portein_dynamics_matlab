from __future__ import division
from __future__ import print_function

import os, sys

import propka.pdb, propka.version, propka.output, propka.conformation_container, propka.group, propka.lib
import propka.lib
import pandas as pd

options,_ = propka.lib.loadOptions()
def getpKa(pdbPath,outType='pd'):
    molecularContainer = Molecular_container_no_output_file(pdbPath,options)
    for name in molecularContainer.conformation_names:
        molecularContainer.conformations[name].calculate_pka(molecularContainer.version,options)
    molecularContainer.find_non_covalently_coupled_groups()
    molecularContainer.average_of_conformations()
    groups = molecularContainer.conformations['AVR'].groups
    propkaData= [(i.residue_type,i.atom.resNumb,i.atom.chainID,i.pka_value) for i in groups]
    if outType in ['pandas','pd']:
        output = pd.DataFrame(propkaData,columns=['ResName','ResID','ChainID','pKa'])
    elif outType in ['np','numpy']:
        output = np.array(propkaData)
    else:
        output = propkaData
    return output

class Molecular_container_no_output_file(propka.molecular_container.Molecular_container):
    def __init__(self, input_file, options=None):
         # printing out header before parsing input
        propka.output.printHeader()

        # set up some values
        self.options = options
        self.input_file = input_file
        self.dir = os.path.split(input_file)[0]
        self.file = os.path.split(input_file)[1]
        self.name = self.file[0:self.file.rfind('.')]
        input_file_extension = input_file[input_file.rfind('.'):]

        # set the version
        if options:
            parameters = propka.parameters.Parameters(self.options.parameters)
        else:
            parameters = propka.parameters.Parameters('propka.cfg')
        try:
            exec('self.version = propka.version.%s(parameters)'%parameters.version)
        except:
            raise Exception('Error: Version %s does not exist'%parameters.version)

        # read the input file
        if input_file_extension[0:4] == '.pdb':
            # input is a pdb file
            # read in atoms and top up containers to make sure that all atoms are present in all conformations
            [self.conformations, self.conformation_names] = propka.pdb.read_pdb(input_file, self.version.parameters,self)
            if len(self.conformations)==0:
                sys.exit(-1)

            self.top_up_conformations()

            # make a structure precheck
            propka.pdb.protein_precheck(self.conformations, self.conformation_names)

            # set up atom bonding and protonation
            self.version.setup_bonding_and_protonation(self)

            # Extract groups
            self.extract_groups()

            # sort atoms
            for name in self.conformation_names:
                self.conformations[name].sort_atoms()

            # find coupled groups
            self.find_covalently_coupled_groups()

            # write out the input file
            filename = self.file.replace(input_file_extension,'.propka_input')
            #propka.pdb.write_input(self, filename)


        elif input_file_extension == '.propka_input':
            #input is a propka_input file
            [self.conformations, self.conformation_names] = propka.pdb.read_input(input_file, self.version.parameters, self)

            # Extract groups - this merely sets up the groups found in the input file
            self.extract_groups()

            # do some additional set up
            self.additional_setup_when_reading_input_file()

        else:
            sys.exit(-1)