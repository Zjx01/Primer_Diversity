import Primer_Diversity.primer_diversity_check as pdc
import pandas as pd
import numpy as np
import logging 
import os 



def get_all_diverse_primers(amplicon_df, primer_names,raw_fasta_file,primer_len_dict,primer_sequence_dict,file_path='all_primers_diverse_check_results'):
    
    if not os.path.exists(file_path):
        try:
            os.makedirs(file_path)
            os.chdir(file_path)
        except:
            log_message = f'Error in creating directory: {file_path}'
            logging.error(log_message)
            return

    logging.basicConfig(filename='primer_diversity.log', level=logging.INFO, 
                format='%(asctime)s - %(levelname)s - %(message)s')
        
    for primer in primer_names:
        fp_name = primer[0]
        rp_name = primer[1]
        sequence_type = primer[0].split('_')[0]
        cur_primer_num = primer[0].split('_')[1]
        primer_specific_df=amplicon_df[(amplicon_df['f_id']==fp_name) & (amplicon_df['r_id']==rp_name)]
        if primer_specific_df.empty:
            log_message = f'No amplicon found for primer pair: {primer}'
            logging.warning(log_message)  # Log as a warning
            return

        else:
            fp_len=primer_len_dict[fp_name]
            rp_len=primer_len_dict[rp_name]
            forward_primer = primer_sequence_dict[fp_name]
            reverse_primer = primer_sequence_dict[rp_name]
            fp_counts_mat,rp_counts_mat = pdc.get_primer_sequence(primer_specific_df,raw_fasta_file,fp_len,rp_len)

            fp_diverse_check = pdc.get_primer_diversity(fp_counts_mat,threshold=0.3,nucleotide_number=2)
            fp_mismatch_check = pdc.get_primer_mismatch(fp_counts_mat,forward_primer)

            rp_diverse_check = pdc.get_primer_diversity(rp_counts_mat,threshold=0.3,nucleotide_number=2)
            rp_mistmatch_check = pdc.get_primer_mismatch(rp_counts_mat,reverse_primer)
            

            if fp_diverse_check > 0 or fp_diverse_check > 0:
                if fp_mismatch_check > 0:
                    log_message = f'Mismatch in Forward Primer {cur_primer_num}: {forward_primer}'
                    logging.info(log_message)
                if fp_diverse_check > 0:
                    log_message = f'Diverse in Forward Primer {cur_primer_num}: {forward_primer}'
                    logging.info(log_message)
                pdc.plot_matrix(fp_counts_mat, forward_primer, cur_primer_num, sequence_type, primer_type='Forward',sc=False, file_path=file_path)
            else:
                log_message = f'Forward Primer {cur_primer_num} {forward_primer} is conserved'
                logging.info(log_message)

            if rp_diverse_check > 0 or rp_mistmatch_check > 0:
                if rp_mistmatch_check > 0:
                    log_message = f'Mismatch in Reverse Primer {cur_primer_num}: {reverse_primer}'
                    logging.info(log_message)
                if rp_diverse_check > 0:
                    log_message = f'Diverse in Reverse Primer {cur_primer_num}: {reverse_primer}'
                    logging.info(log_message)
                pdc.plot_matrix(rp_counts_mat, reverse_primer, cur_primer_num, sequence_type, primer_type='Reverse', sc=False, file_path=file_path)
            else:
                log_message = f'Reverse Primer {cur_primer_num} {reverse_primer} is conserved'
                logging.info(log_message)

    return    

