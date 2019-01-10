import os
import Bio.PDB
import pandas as pd

dssp_exec ='/opt/dssp/dssp'
def getDSSP(PDBPath):
    idx = pd.IndexSlice
    
    # Use DSSP to assign secondary structure type to each residue.
    pdb_parser = Bio.PDB.PDBParser()
    structure = pdb_parser.get_structure('PDB', PDBPath)
    model = structure[0]
    dssp_results = Bio.PDB.DSSP(model, PDBPath, dssp=dssp_exec)

    chain_ID_resid_list = dssp_results.keys()
    dssp_results_table = pd.DataFrame(list(dssp_results))
    dssp_results_table = dssp_results_table[[0,1,2,3,4,5]]
    dssp_results_table.columns = ['DSSP Index', 'Residue Type', 'SS', 'Acc', 'Phi', 'Psi']
    dssp_results_table.index = pd.MultiIndex.from_tuples([(chain_ID_resid_list[x][0], chain_ID_resid_list[x][1][1]) for x in range(len(chain_ID_resid_list))])
    dssp_results_table = dssp_results_table.reindex()

    # Make sure all helices have negative phi and psi angles, otherwise change them to '-'
    helix_only_rows = (dssp_results_table['SS'] == 'H')
    helix_only_table = dssp_results_table[helix_only_rows]
    pos_phi_or_psi_idx = helix_only_table[(helix_only_table['Phi'] >= 0).values + (helix_only_table['Psi'] >= 0).values].index.tolist()
    dssp_results_table.loc[idx[:,pos_phi_or_psi_idx], 'SS'] = '-'
    return dssp_results_table