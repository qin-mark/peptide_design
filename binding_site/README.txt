Input:
  -f PDBFILE, --pdbfile=PDBFILE
                        Input PDB file path.
Setting:
  -b ANCHOR_BEFORE_LEN, --anchor_before_len=ANCHOR_BEFORE_LEN
                        Input anchor before length.
  -t ANCHOR_AFTER_LEN, --anchor_after_len=ANCHOR_AFTER_LEN
                        Input anchor after length.
  -p PEPTIDE_CHAIN, --peptide_chain=PEPTIDE_CHAIN
                        Input anchor after length.
  -a AKT_CHAIN, --AKT_chain=AKT_CHAIN
                        Input AKT chain.
  -s AKT_DOMAIN_START, --AKT_domain_start=AKT_DOMAIN_START
                        Input AKT domain start.
  -e AKT_DOMAIN_END, --AKT_domain_end=AKT_DOMAIN_END
                        Input AKT domain end.
  -o CA_ONLY, --CA_only=CA_ONLY
                        CA only.
  -c cutoff, --cutoff=cutoff
                        Input cutoff.

Example:
python Binding.py -f ./EHGAIPRLVQLLMRAHQDTQRRTSMASS_unrelaxed_rank_1_model_3.pdb
python Binding.py --pdbfile ./EHGAIPRLVQLLMRAHQDTQRRTSMASS_unrelaxed_rank_1_model_3.pdb

output:
binding ratio of AKT domain: 0.02681992337164751
binding ratio of AKT: 0.014583333333333334
binding ratio of peptide: 0.2857142857142857