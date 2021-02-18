import glob
import os
import subprocess
import sys
from pele_platform.constants import constants as cs
import pele_platform.Utilities.Helpers.helpers as hp
from rdkit import Chem
from rdkit.Chem import rdmolfiles
from rdkit.Chem.AtomPairs import Pairs
from itertools import combinations
from rdkit.Chem import AllChem
import tempfile

OUTPUT = "input.conf"

'''
def growing_sites_symmetry(fragment,user_bond):
    bonds = []
    mol = Chem.MolFromPDBFile(fragment, removeHs=False)
    #mol2 = Chem.MolToSmiles(mol)
    symmetryList = []
    counter = 0
    for atom in mol.GetAtoms():
        fing = AllChem.GetMorganFingerprint(atom,2)
        print(fing)
        previouscounter = counter
        idx = atom.GetIdx()
        neighborList = [neighbor.GetIdx() for neighbor in atom.GetNeighbors()]
        if len(neighborList) < 2:
            continue
        neighborCombinations = list(combinations(neighborList,2))
        for neighborSet in neighborCombinations:
            NewMol = Chem.RWMol(mol) 
            NewMol.RemoveBond(neighborSet[0],idx)
            NewMol.GetAtomWithIdx(neighborSet[0]).SetNoImplicit(True)
            NewMol.RemoveBond(neighborSet[1],idx)
            NewMol.GetAtomWithIdx(neighborSet[1]).SetNoImplicit(True)
            NewMol.GetAtomWithIdx(idx).SetNoImplicit(True)
            Chem.SanitizeMol(NewMol)
            Fragmented_Smiles = str(Chem.MolToSmiles(NewMol))
            FragA, FragB, FragC = sorted(Fragmented_Smiles.split('.'))
            if (FragA == FragB) or (FragA == FragC ) or (FragB == FragC):
                if previouscounter == counter: 
                    symmetryList.append([idx])
                    symmetryList[counter].append(neighborSet)
                    counter += 1
                else:
                    symmetryList[counter - 1].append(neighborSet)
    print(symmetryList)
    sys.exit()
    if mol:
        heavy_atoms = [a for a in mol.GetAtoms() if a.GetSymbol() != "H"]
        for a in heavy_atoms:
            hydrogens = [n for n in a.GetNeighbors() if n.GetSymbol() == "H" and n.GetIdx() in number]
            at_name = a.GetMonomerInfo().GetName().strip()
            for h in hydrogens:
                if h not in matches:
                    h_name = h.GetMonomerInfo().GetName().strip()
                    bonds.append("{} {} {}-{}".format(fragment, user_bond, at_name, h_name))
    return bonds
    
'''
    
def growing_sites(fragment, user_bond):
    """
    Retrieves all possible growing sites (hydrogens) on the fragment. Takes PDB fragment file as input.
    Output - list of strings represeting sites, e.g. "benzene.pdb C6-H6 C1-H2"
    """
    if hp.is_rdkit():
        from rdkit import Chem

    bonds = []
    
    #fragment_pdb= Chem.PDBWriter(fragment.name)
        
    mol = Chem.MolFromPDBFile(fragment.name, removeHs=False)
    
    if mol:
        heavy_atoms = [a for a in mol.GetAtoms() if a.GetSymbol() != "H"]
        for a in heavy_atoms:
            hydrogens = [n for n in a.GetNeighbors() if n.GetSymbol() == "H"]
            at_name = a.GetMonomerInfo().GetName().strip()
            for h in hydrogens:
                h_name = h.GetMonomerInfo().GetName().strip()
                bonds.append("{} {} {}-{}".format(fragment, user_bond, at_name, h_name))
    
    return bonds


def sdf_to_pdb(file_list, path, logger):

    out = []
    fo = tempfile.NamedTemporaryFile()
    if file_list:
        converted_mae = []
        output = []

        # convert all SDF to MAE
        schrodinger_path = os.path.join(cs.SCHRODINGER, "utilities/structconvert")
        command_mae = "{} -isd {} -omae {}"
        command_pdb = "{} -imae {} -opdb {}"
        #with tempfile.TemporaryFile(mode = 'w+') as fp:
        for f in file_list:
            fout = os.path.splitext(os.path.basename(f))[0] + ".mae"
            fout_path = os.path.join(os.path.dirname(f), fout)
            try:
                command_mae = command_mae.format(schrodinger_path, f, fout_path)
                #fp = tempfile.TemporaryFile()
                #subprocess.Popen(command_mae.split(), stdout=fp, stderr=subprocess.STDOUT)
                subprocess.call(command_mae.split())
                converted_mae.append(fout_path)
            except Exception as e:
                logger.info("Error occured while converting SD files to mae.", e)
         
        # convert all MAE to PDB, it will result in a lot of numbered pdb files
        for c in converted_mae:
            
            fout = c.replace(".mae",".pdb")
            try:
                command_pdb = command_pdb.format(schrodinger_path, c, fout)
                subprocess.call(command_pdb.split())
                os.remove(c)
            except Exception as e:
                logger.info("Error occured while converting mae to PDB.", e)
            
        pdb_pattern = os.path.splitext(converted_mae[0])
        converted_pdb = glob.glob(pdb_pattern[0]+"*"+".pdb")        
        # ~~~ If it's stupid but it works (?), it isn't stupid. ~~~
        
        # read in PDB file created by Schrodinger, substitute residue name and add chain ID
        out = []
        fp = tempfile.NamedTemporaryFile(delete=False, suffix='.pdb')
        for c in converted_pdb:
            with open(c, "r") as fin:
                lines = fin.readlines()
                new_lines = []
                for line in lines:
                    if line.startswith("HETATM") or line.startswith("ATOM"):
                        new_lines.append(line)
            
            new_lines = [l.replace("UNK", "GRW") for l in new_lines if "UNK" in l]
            new_lines = [l[:21]+"L"+l[22:] for l in new_lines]

            with open(c, "r+") as fout:
                for line in new_lines:
                    fout.write(line)
                for line in fout.readlines():
                    fp.write(line.encode()) 
            out.append(fp)
            os.remove(c)
        #out[0].seek(0)
        #print(out[0].read())
        #sys.exit()
        #out =  converted_pdb
    return out

def get_library(frag_library):
    directory = os.path.dirname(os.path.abspath(__file__))                                                                                                                              
    path = frag_library if os.path.exists(frag_library) else os.path.join(directory, "Libraries", frag_library.strip())                                                                 
    if not os.path.exists(path):                                                                                                                                                        
        raise OSError(f"File {frag_library} doesn't exist and is not one of our internal libraries. Please check the frag_library flag in input.yaml.")   
    return path


def get_fragment_files(path, logger):

    fragment_files = []                                                                                                                                                                 
    extensions = ['*.pdb', '*.sdf']                                                                                                                                                     
    
    for e in extensions:                                                                                                                                                                
        fragment_files.extend(glob.glob(os.path.join(path, e.upper())))                                                                                                                 
        fragment_files.extend(glob.glob(os.path.join(path, e.lower())))                                                                                                                 
                                                                                                                                                                                        
    # convert SDF to PDB, if necessary                                                                                                                                                  
    sdf_files = [elem for elem in fragment_files if ".sdf" in elem.lower()]                                                                                                             
    pdb_files = [elem for elem in fragment_files if ".pdb" in elem.lower()]                                                                                                             
    all_files = pdb_files + sdf_to_pdb(sdf_files, path, logger)                                                                                                                                                                               
    return all_files


def write_config_file(output_name, bond_list):

    with open(output_name, "w+") as conf_file:
        for line in bond_list:
            conf_file.write(line+"\n")


def main(user_bond, frag_library, logger):

    # find the library and extract fragments
    path = get_library(frag_library)
    all_files = get_fragment_files(path, logger) 
    
    # get all possible growing sites
    bond_list = []
    for file in all_files:
        
        bond_list.extend(growing_sites(file, user_bond))
    
    # write input.conf 
    write_config_file(OUTPUT, bond_list)
    
    return OUTPUT
