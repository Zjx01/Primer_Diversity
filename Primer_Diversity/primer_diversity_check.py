import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio import AlignIO
import logomaker as lm
import os 



def read_amplicon_hits(file):
    df = pd.read_csv(file, sep='\t')
    return df


def read_primer_design(file):
    primer_df = pd.read_csv(file, sep='\t')
    primer_len_dict={}
    primer_sequence_dict={}
    for name in primer_df['name']:
        primer_len_dict[name] = primer_df.loc[primer_df['name']==name]['size'].values[0]
        if name.endswith('LEFT'):
            primer_sequence_dict[name]= primer_df.loc[primer_df['name']==name].seq.values[0]
        else:
            primer_sequence_dict[name]= reverse_complement(primer_df.loc[primer_df['name']==name].seq.values[0])
    return primer_df,primer_len_dict, primer_sequence_dict


def get_primer_number(amplicon_df):
    primer_number = amplicon_df['f_id'].nunique()
    return primer_number


def get_primer_name(primer_df,primer_number):
    base_primer_name = primer_df['name'].str.split('_').str[0][0]
    primer_names = []
    for i in range(1,primer_number+1):
        pnumber = i
        primer_name_left = base_primer_name + '_' + str(pnumber) + '_LEFT'
        primer_name_right = base_primer_name + '_' + str(pnumber) + '_RIGHT'
        primer_names.append((primer_name_left,primer_name_right))
    return base_primer_name,primer_names



def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'R':'Y','Y':'R','N': 'N', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n','y':'r'}
    return "".join(complement[base] for base in reversed(seq)).upper()
            


def get_primer_sequence(primer_specific_df,raw_fasta_file,fp_len,rp_len):
    primer_forward_matrix = []
    primer_reverse_matrix = []
    for record in SeqIO.parse(raw_fasta_file, "fasta"):
        if record.id in primer_specific_df['target_id'].values:
            sequence = str(record.seq)
            primer_sequence_info = primer_specific_df[primer_specific_df['target_id']==record.id]
            primer_f_start,primer_f_end = primer_sequence_info['f_start'].values[0],primer_sequence_info['f_end'].values[0]
            primer_r_start,primer_r_end = primer_sequence_info['r_start'].values[0],primer_sequence_info['r_end'].values[0]
            
            primer_f_sequence = sequence[primer_f_start:primer_f_end+1].upper()
            primer_r_sequence = sequence[primer_r_start:primer_r_end+1].upper()
            
            # here, we only use the primer sequence with the exact length of arctic primer 
            if len(primer_f_sequence) == fp_len:
                primer_forward_matrix.append(primer_f_sequence)
            if len(primer_r_sequence) == rp_len:
                primer_reverse_matrix.append(primer_r_sequence)


    fp_counts_mat = lm.alignment_to_matrix(primer_forward_matrix)
    rp_counts_mat = lm.alignment_to_matrix(primer_reverse_matrix)
            
    return fp_counts_mat,rp_counts_mat


def primer_similarity(primer1,primer2):
    mismatch_score = sum([i!=j for i,j in zip(primer1,primer2)])
    return mismatch_score


def get_primer_mismatch(counts_mat,primer):
    max_freq = ''.join(counts_mat.idxmax(axis=1))
    primer_mismatch = primer_similarity(primer,max_freq)
    return primer_mismatch


def get_primer_diversity(counts_mat,threshold=0.3,nucleotide_number=2):
    counts_mat = counts_mat.div(counts_mat.sum(axis=1), axis=0)
    each_nucleotide_evaluate = lambda x: len([i for i in x if i > threshold]) 
    diverse_number = counts_mat.apply(each_nucleotide_evaluate, axis=1)
    nucleotide_all_evaluate = sum(diverse_number.values>=nucleotide_number)
    return nucleotide_all_evaluate

    
    
    
def plot_matrix(primer_counts_mat,primer,specific_primer_number,sequence_type,file_path,primer_type='forward',sc=True):
    lm.Logo(primer_counts_mat,color_scheme='classic')
    plt.title(sequence_type+' Scheme '+str(specific_primer_number)+' '+primer_type+' Primer Nuclotide Logo Plot')
    positions = range(len(primer))
    plt.xticks(positions, list(primer));

    if file_path is None:
        file_path = os.getcwd()
        
    if os.path.exists(file_path):
        if file_path[-1] != '/':
            file_path = file_path+'/'
    else:
        os.makedirs(file_path)

    scheme_name = sequence_type+'_scheme_'+str(specific_primer_number)+'_'+primer_type+'_primer_logo_plot.png'
    plt.savefig(file_path+scheme_name)
    if sc == True:
        plt.show()
    plt.close()
    return





