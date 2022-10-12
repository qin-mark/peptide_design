from Bio.PDB import PDBParser


def binding_ratio(pdbfile, mcherrystart_site=29, peptide_chain='B', AKT_chain='C', AKT_domain_start=149,
                  AKT_domain_end=409, CA_only=False):
    '''
    pdbfile: pdbfile position for AKT and peptide complex.
    in result pdb file, peptide and mchaeery are concatenate, and mcherrystart_site is the start site for mcherry.
    peptide_chain indicate the peptide chain_AKT
    AKT_chain indicate the AKT chain
    if mcherrystart_site=-1, no mcherry
    AKT_domain_start: domain start position on AKT
    AKT_domain_end: domain end position on AKT
    CA_only: whether consider all atoms or CA only, default is False same with current ppt
    '''
    parser = PDBParser()
    # read structure from file
    structure = parser.get_structure('peptide_AKT', pdbfile)

    model = structure[0]
    chain_peptide = model[peptide_chain]
    chain_AKT = model[AKT_chain]

    bind_num = 0
    cutoff = 5
    peptide_bind_sites = {}
    AKT_bind_sites = {}
    AKT_bind_sites_domain = {}
    for residue1 in chain_peptide:
        peptide_pos = residue1.id[1]
        if mcherrystart_site != -1:
            if peptide_pos >= mcherrystart_site:
                break

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
    # print(num_of_AKT_bind_sites_domain)
    # print(num_of_AKT_bind_sites)
    return num_of_AKT_bind_sites_domain / (AKT_domain_end - AKT_domain_start + 1), num_of_AKT_bind_sites / 480
