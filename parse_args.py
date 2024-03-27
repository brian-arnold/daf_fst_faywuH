import sys
import os
import argparse

def parse_args():
    description = "FST."
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--winsize', type=int, required=True, help="Window size in SNPs, must be integer > 0")
    parser.add_argument('--offset', type=int, required=False, default=None, help="The offset between windows in SNPs, must be integer, default is winsize")
    parser.add_argument('--min_AF', type=int, required=False, default=1, help="minimum allele frequency as integer, e.g. 3 means the ALT allele must have been detected on at least 3 chromosomes")
    parser.add_argument('--pop_file', type=str, required=True, help="a file of individuals from all 3 populations, names must match those in VCF; column1 has ind name, column2 has population index")
    parser.add_argument('--focal_pop', type=str, required=True, help="The population to use as the focal population to compute DAF, FST, and Fay-Wu H statistics")
    parser.add_argument('--sister_pop', type=str, required=True, help="The population to use as the sister to the focal population, more closely related than the outgroup population, used to compute DAF, FST")
    parser.add_argument('--outgroup_pop', type=str, required=True, help="The population to use for polarizing mutations and computing DAF")
    
    parser.add_argument('--output_fst', type=str, required=False, default = "./FST.txt", help="Name of the output file for the FST statistic")
    parser.add_argument('--output_faywuh', type=str, required=False, default = "./FayWuH.txt", help="Name of the output file for the FayWuH statistic")
    parser.add_argument('--output_ddaf', type=str, required=False, default = "./DDAF.txt", help="Name of the output file for DDAF statistic")
    
    parser.add_argument('--vcf', type=str, required=True, help="the VCF file to analyze")

    return parser.parse_args()

args = parse_args()
if args.offset is None:
    args.offset = args.winsize



