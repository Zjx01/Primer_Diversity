import argparse
from ast import arg
from .primer_diversity_check import *
from .multi_scheme_diversity_check import get_all_diverse_primers
from .single_scheme_diversity_check import get_specific_primers


def get_parser():
    parser=argparse.ArgumentParser(description="Output diversity plot for primer scheme")

    ### neccessary arguments
    parser.add_argument("-f", "--fasta", help="Input raw fasta sequence file with primer scheme", required=True)
    parser.add_argument("-p", "--primer_design", help="primer_design file conting primer names and primer sequences", required=True)
    parser.add_argument("-a", "--amplicon_hits", help="Amplicon hits file", required=True)
    parser.add_argument("--singles_scheme", "-sc", help="Output file with single scheme diversity plot", action='store_true')
    parser.add_argument('--multiple_scheme','-mc', help="Output file with multiple scheme diversity plot", action='store_true')
    parser.add_argument("--output", "-o", help="output directory for the primer diversity plot", default="None")

    ### option arguments - single scheme
    parser.add_argument("-sn", "--scheme_number", help="Number of the primer scheme for visualization", default="None")

    ### optional arguments - multiple schemes
    parser.add_argument("-t", "--threshold", help="threshold for nucleotide diversity (<0.5)", default=0.3)
    parser.add_argument("-n", "--nucleotide_number", help="number of nucleotide diversity", default=2)

    return parser


def main():
    parser = get_parser()
    args, unknown = parser.parse_known_args()
    amplicon_df=args.amplicon_hits
    raw_fasta_file=args.fasta
    primer_design= args.primer_design
    file_path=args.output


    amplicon_df= read_amplicon_hits(amplicon_df)
    primer_df,primer_len_dict, primer_sequence_dict = read_primer_design(primer_design)
    primer_number=get_primer_number(amplicon_df)
    base_primer_names,primer_names=get_primer_name(primer_df,primer_number)

    if args.singles_scheme:
        scheme_number = args.scheme_number
        get_specific_primers(amplicon_df,scheme_number, primer_df,primer_number,raw_fasta_file,primer_len_dict,primer_sequence_dict,file_path)

    if args.multiple_scheme:
        get_all_diverse_primers(amplicon_df, primer_names,raw_fasta_file,primer_len_dict,primer_sequence_dict,file_path)

                

if __name__ == '__main__':
    main()

