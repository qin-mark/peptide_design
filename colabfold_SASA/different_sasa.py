import os
import uuid
from typing import Optional, Dict
from Bio.PDB import PDBParser, SASA
from Bio import PDB


class SurfaceArea:
    def __init__(self, probe_radius: float = 1.4, n_points: int = 100, radii_dict: Optional[Dict] = None):
        if radii_dict is None:
            radii_dict = {'X': 2.0}

        self.parser = PDBParser(QUIET=1)
        self.structure_computer = SASA.ShrakeRupley(probe_radius=probe_radius, n_points=n_points, radii_dict=radii_dict)

    def __call__(self, name, loc) -> float:
        struct = self.parser.get_structure(name, loc)
        # self.structure_computer.compute(struct, level="C")
        # return struct[0]['C'].sasa
        self.structure_computer.compute(struct, level="S")
        return struct.sasa


class ChainSplitter:
    def __init__(self, out_dir=None):
        """ Create parsing and writing objects, specify output directory. """
        self.parser = PDB.PDBParser()
        self.writer = PDB.PDBIO()
        if out_dir is None:
            out_dir = os.path.join(os.getcwd(), "chain_PDBs")
        self.out_dir = out_dir

    def make_pdb(self, pdb_path, chain_letters, overwrite=False, struct=None):

        pdb_id = self.out_dir.split('\\')[-1]
        # Get structure, write new file with only given chains
        if struct is None:
            struct = self.parser.get_structure(pdb_id, self.out_dir)
        self.writer.set_structure(struct)
        self.writer.save(pdb_path, select=SelectChains(chain_letters))

        return pdb_path


class SelectChains(PDB.Select):
    """ Only accept the specified chains when saving. """
    def __init__(self, chain_letters):
        self.chain_letters = chain_letters

    def accept_chain(self, chain):
        return (chain.get_id() in self.chain_letters)


def get_complex_pdb_name(path, complex_name):
    name_list = os.listdir(path)
    pdb_name = ''
    for name in name_list:
        # TODO add file name befor the 'complex_name' use '//' to split
        if name.startswith(complex_name + '_unrelaxed_rank_1') and name[-4:] == '.pdb':
            pdb_name = name
    if pdb_name == '':
        return ValueError
    else:
        return pdb_name


def get_sasa_difference(complex_name):
    # TODO replace the path
    path = os.getcwd() + '\\test_pdb\\'
    complex_pdb_name = get_complex_pdb_name(path, complex_name)
    # print(complex_pdb_name)
    piptide_pdb_name = 'piptide.pdb'
    akt_pdb_name = 'akt.pdb'
    akt_sasa = 28407.02
    splitter = ChainSplitter(path + complex_pdb_name)
    splitter.make_pdb(path + piptide_pdb_name, 'B')
    sasa_fn = SurfaceArea()
    complex_sasa = sasa_fn(uuid.uuid4().hex, path + complex_pdb_name)
    independent_sasa = sasa_fn(uuid.uuid4().hex, path + piptide_pdb_name)
    # akt_sasa = sasa_fn(uuid.uuid4().hex, path + akt_pdb_name)
    # print(akt_sasa, independent_sasa, complex_sasa)
    return akt_sasa + independent_sasa - complex_sasa


if __name__ == '__main__':
    print('%.2f' % get_sasa_difference('test_9216f'))