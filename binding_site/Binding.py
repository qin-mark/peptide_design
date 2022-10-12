from os import path
from optparse import OptionParser
from Bio.PDB import PDBParser


def binding_ratio(cutoff,pdbfile, anchor_before_len=0, anchor_after_len=0, peptide_chain='B', AKT_chain='C',
                  AKT_domain_start=149,
                  AKT_domain_end=409, CA_only=False):
    '''
    pdbfile: pdbfile position for AKT and peptide complex.
    in result pdb file, peptide and anchor_before,anchor_after are concatenate
    peptide_len: the length of original peptide
    anchor_before_len: the length of anchor_before
    anchor_after_len: the length of anchor_after
    peptide_chain indicate the peptide chain_AKT
    AKT_chain indicate the AKT chain
    AKT_domain_start: domain start position on AKT
    AKT_domain_end: domain end position on AKT
    CA_only: whether consider all atoms or CA only, default is False same with current ppt
    '''
    parser = PDBParser()
    structure = parser.get_structure('peptide_AKT', pdbfile)

    model = structure[0]
    chain_peptide = model[peptide_chain]
    chain_AKT = model[AKT_chain]

    peptide_len=len(chain_peptide)-anchor_before_len-anchor_after_len
    bind_num = 0
    cutoff = float(cutoff)
    peptide_bind_sites = {}
    AKT_bind_sites = {}
    AKT_bind_sites_domain = {}
    for residue1 in chain_peptide:
        peptide_pos = residue1.id[1]
        if anchor_before_len != 0:
            if peptide_pos <= anchor_before_len:  # peptide_pos in anchor_before
                continue

        if anchor_after_len != 0:
            if peptide_pos > anchor_before_len + peptide_len:  # peptide_pos in anchor_after
                break

        # print(peptide_pos)
        for residue2 in chain_AKT:
            AKT_pos = residue2.id[1]
            if CA_only:
                # compute distance between CA atoms
                try:
                    distance = residue1['CA'] - residue2['CA']
                except KeyError:
                    ## no CA atom, e.g. for H_NAG
                    print("no ca ")
                    continue

                if distance <= cutoff:
                    peptide_bind_sites[peptide_pos] = distance
                    AKT_bind_sites[AKT_pos] = distance
                    if AKT_pos >= AKT_domain_start and AKT_pos <= AKT_domain_end:
                        AKT_bind_sites_domain[AKT_pos] = distance

            else:
                # compute distance between all atoms
                for atom1 in residue1.get_atoms():
                    for atom2 in residue2.get_atoms():
                        distance = atom1 - atom2
                        if distance <= cutoff:
                            peptide_bind_sites[peptide_pos] = distance
                            AKT_bind_sites[AKT_pos] = distance
                            if AKT_pos >= AKT_domain_start and AKT_pos <= AKT_domain_end:
                                AKT_bind_sites_domain[AKT_pos] = distance

                            break

    num_of_AKT_bind_sites = len(AKT_bind_sites)
    num_of_peptide_bind_sites = len(peptide_bind_sites)
    num_of_AKT_bind_sites_domain = len(AKT_bind_sites_domain)
    print('binding sites of AKT domain: '+str(AKT_bind_sites_domain))
    print('binding sites of AKT: '+str(AKT_bind_sites))
    print('binding sites of AKT peptide: '+str(peptide_bind_sites))
    return num_of_AKT_bind_sites_domain / (AKT_domain_end - AKT_domain_start + 1), num_of_AKT_bind_sites / len(chain_AKT),num_of_peptide_bind_sites/peptide_len



def main():
    """

    """
    # Construct the usage.
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("-f", "--pdbfile", type="string",dest="pdbfile", default=None,
                      help="Input PDB file path.")
    parser.add_option("-b", "--anchor_before_len", type="int", dest="anchor_before_len", default=0,
                      help="Input anchor before length.")
    parser.add_option("-t", "--anchor_after_len", type="int", dest="anchor_after_len", default=0,
                      help="Input anchor after length.")
    parser.add_option("-p", "--peptide_chain", type="string", dest="peptide_chain", default="B",
                      help="Input anchor after length.")
    parser.add_option("-a", "--AKT_chain", type="string", dest="AKT_chain", default="C",
                      help="Input AKT chain.")
    parser.add_option("-s", "--AKT_domain_start", type="int", dest="AKT_domain_start", default=149,
                      help="Input AKT domain start.")
    parser.add_option("-e", "--AKT_domain_end", type="int", dest="AKT_domain_end", default=409,
                      help="Input AKT domain end.")
    parser.add_option("-o", "--CA_only", dest="CA_only", default=False,
                      help="CA only.")
    parser.add_option("-c", "--cutoff", dest="cutoff", default=5,
                      help="Input cutoff.")

    # Parse the options and args input by users.
    (args, options) = parser.parse_args()
    # args = parser.parse_args()

    filePath = path.normpath(args.pdbfile)
    if not path.isfile(filePath) or not filePath.endswith(".pdb"):
        raise Exception("The input file is not exist or a available PDB file.")
    num1, num2,num3 = binding_ratio(args.cutoff,args.pdbfile,  args.anchor_before_len, args.anchor_after_len, args.peptide_chain,
                  args.AKT_chain, args.AKT_domain_start, args.AKT_domain_end, args.CA_only)
    print('binding ratio of AKT domain: ' + str(num1))
    print('binding ratio of AKT: ' + str(num2))
    print('binding ratio of peptide: ' + str(num3))


if __name__ == "__main__":
    main()