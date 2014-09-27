__author__ = 'Patrick Grinaway'

import openeye.oechem as oechem
import openeye.oeiupac as oeiupac
import openeye.oeomega as oeomega
import openeye.oequacpac as oequacpac
import openeye.oegrapheme as oegrapheme
import copy




def _gen_from_iupac(iupac_name):

    iupac_mol = oechem.OEMol()
    oeiupac.OEParseIUPACName(iupac_mol,iupac_name)
    oechem.OEAddExplicitHydrogens(iupac_mol)
    omega = oeomega.OEOmega()
    omega.SetMaxConfs(1)
    omega(iupac_mol)
    return iupac_mol


def _read_mol2(mol2_filename):
    mol = oechem.OEMol()
    istream = oechem.oemolistream()
    istream.open(mol2_filename)
    oechem.OEReadMol2File(istream, mol)
    return mol

def _atom_name_list(list_of_names):
    return "PX = TitratableResidue('p-xylene', %s, pKa=7.0)\n" % str(list_of_names)


def _generate_charge_python(list_of_charges, ref_ene, prot_state=1):
    return "PX.add_state(protcnt=%d, refene=%s, charges=%s)\n" % (prot_state, ref_ene, str(list_of_charges))

def _generate_zero_refenergies(num_refene):
    refene_string = ""
    for i in range(1,num_refene+1):
        str_to_add = """refene%d = _ReferenceEnergy(igb2=0, igb5=0)\n
        refene%d.solvent_energies()\n
        refene%d.dielc2_energies(igb2=0, igb5=0, igb8=0)\n
        refene%d.dielc2.solvent_energies()\n""" % (i, i, i, i)
        refene_string+=str_to_add
    return refene_string

def generate_fluorinated_mols(starting_mol):
    """
    Takes an organic molecule, and creates a list of molecules with each hydrogen substituted with a fluorine

    Parameters
    ----------
    starting_mol : OEMol
        OEMol containing the molecule to fluorine scan. Need not have explicit hydrogens

    Returns
    -------
    fluorinated_list : list of OEMol
        list of OEMol objects where each hydrogen in starting_mol is replaced by a fluorine

    >>>import openeye.oechem as oechem
    >>>import openeye.iupac as oeiupac
    >>>mol = oechem.OEMol()
    >>>oeiupac.OEParseIUPACName(mol, "p-xylene")
    """
    fluorinated_list = []
    #make sure explicit hydrogens are present
    oechem.OEAddExplicitHydrogens(starting_mol)

    #get locations of hydrogens:
    h_list = [atm.GetIdx() for atm in starting_mol.GetAtoms(oechem.OEIsHydrogen())]

    #iterate through list, creating new mol and replacing each H with F:
    for hydrogen in h_list:
        mol_to_modify = copy.deepcopy(starting_mol)
        mol_to_modify.GetAtom(oechem.OEHasAtomIdx(hydrogen)).SetAtomicNum(oechem.OEElemNo_F)
        fluorinated_list.append(mol_to_modify)
    return fluorinated_list

def generate_charges(mol_to_charge, charge_model=oequacpac.OECharges_AM1BCCSym):
    """
    Takes an organic molecule with geometry, and generates charges using specified charge model

    Parameters
    ----------
    mol_to_charge : OEMol
        An OEMol object containing the molecule to charge. Must have a 3D geometry.
    charge_model (optional) : OECharges_x
        An OECharges specifying the charge model to use. Default: OECharges_AM1BCCSym

    Returns
    -------
    charged_mol : OEMol
        an OEMol with charges on atoms
    """
    charged_mol = copy.deepcopy(mol_to_charge)
    oequacpac.OEAssignPartialCharges(charged_mol, charge_model, False, False)
    return charged_mol





if __name__ == "__main__":

    mol = _read_mol2('/Users/grinawap/px.mol2')
    atom_names = [atom.GetName() for atom in mol.GetAtoms()]
    fluorinated_list = generate_fluorinated_mols(mol)
    p_charge_molecules = [generate_charges(molecule) for molecule in fluorinated_list]

    #iterate through the charges and write them out:
    chargelist = []
    for molecule in p_charge_molecules:
        chargelist.append([atom.GetPartialCharge() for atom in molecule.GetAtoms()])
    outstring = ""
    i=1
    for charges in chargelist:
        outstring+=_generate_charge_python(charges, "ref_ene%d"%i)
        i+=1
    name_list = _atom_name_list(atom_names)
    name_list+=outstring
    outf = open('ex_ens.py','w')
    outf.writelines(name_list)
    outf.close()



