import Primer_Diversity.primer_diversity_check as pdc 
import pandas as pd
import numpy as np

def get_specific_primers(amplicon_df,scheme_number, primer_df,primer_number,raw_fasta_file,primer_len_dict,primer_sequence_dict,file_path='single_scheme_diverse_check_results'):
    base_primer_name,primer_names = pdc.get_primer_name(primer_df,primer_number)
    fp_name = base_primer_name + '_' + str(scheme_number)+'_LEFT'
    rp_name = base_primer_name + '_' + str(scheme_number)+'_RIGHT'
    primer_specific_df = amplicon_df[(amplicon_df['f_id']==fp_name) & (amplicon_df['r_id']==rp_name)]

    fp_len=primer_len_dict[fp_name]
    rp_len=primer_len_dict[rp_name]
    forward_primer = primer_sequence_dict[fp_name]
    reverse_primer = primer_sequence_dict[rp_name]

    fp_counts_mat,rp_counts_mat = pdc.get_primer_sequence(primer_specific_df,raw_fasta_file,fp_len,rp_len)
    pdc.plot_matrix(fp_counts_mat,forward_primer,scheme_number,base_primer_name,primer_type='Forward',file_path=file_path)
    pdc.plot_matrix(rp_counts_mat,reverse_primer,scheme_number,base_primer_name,primer_type="Reverse",file_path=file_path)
    return 




