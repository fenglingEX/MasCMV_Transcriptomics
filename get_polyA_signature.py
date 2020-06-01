from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import pandas as pd 
import sys
import re
from argparse import ArgumentParser

def get_args():
    parser = ArgumentParser()
    # List arguments
    parser.add_argument("-g", "--genome_file", required = True, help = "Genome file in .fasta format.")
    parser.add_argument("-o", "--output_file", required = True, help = "Output .tsv file to store results in.")
    parser.add_argument("-i", "--input_file", required = True, help = ".bed format file containing polyA coordinates. Output by ContextMap.")
    args = vars(parser.parse_args())
    return args

# Reads in .bed file output by ContextMap and converts it into a pandas dataframe.
def read_in_polya_site_bed_as_dataframe(input_file):
    col_names = ["Chromosome", "Start", "End", "Name", "Score", "Strand"]
    bed_df = pd.read_csv(input_file, sep = "\t", usecols = [0,1,2,3,4,5], names = col_names, index_col = 3)
    # Turn numeric columns into ints.
    bed_df.Start = bed_df.Start.astype('int64')
    bed_df.End = bed_df.End.astype('int64')
    bed_df.Score = bed_df.Score.astype('int64')
    return bed_df


def read_in_genome_sequence(genome_file):
    genome_seq = SeqIO.read(genome_file, "fasta").seq
    return str(genome_seq)

# Extracts set lenght sequence upstream of a polyA site (relative to the transcript's direction)
def extract_upstream_sequence(genome_seq, position, bases, strand):
    # If transcript is plus sense, get preceding bases in sequence
    if strand == '+':
        # PolyA signal present upstream of position in genome
        sequence = genome_seq[(position-bases):position]
    # If transcript is minus sense, get subsequent bases, and reverse complement them.
    elif strand == '-':
        # PolyA signal present downstream of position in genome
        sequence = genome_seq[position:(position+bases)]
        sequence = str(Seq(sequence, IUPAC.unambiguous_dna).reverse_complement())
    return sequence


# Iterates over a given sequence and counts how many times
# the given motif occurs. Return count.
def count_motif_in_sequence(sequence, motif):
    count = 0
    motif_length = len(motif)
    # For each substring with length equal to the motif
    # check if it matches the motif.
    for a in range(len(sequence)-motif_length):
        if sequence[a:a+motif_length] == motif:
            count += 1
    return count

# Counts the occurences of two polyA signals upstream of polyA sites
# At 50 and 100 bases distance.
# Returns a dataframe containing polyA signal counts for each site.
def find_poly_a_signals_for_each_site(polya_sites, genome_seq):
    # For each polyA site (row in dataframe)
    for index, site in polya_sites.iterrows():
        # Get strand
        strand = site["Strand"]
        # Strand determines which position is used to extract upstream sequence, end or start.
        if strand == '+':
            position = site["Start"]
        elif strand == '-':
            position = site["End"]
        # 
        seq_50_base = extract_upstream_sequence(genome_seq, position, 50, strand)
        seq_100_base = extract_upstream_sequence(genome_seq, position, 100, strand)

        # For each motif, count the occurrences within 50 and 100 bases of the polyA site.
        # TODO: This could be turned into a nested for-loop that calls count_motif_in_sequence.
        motif_1 = "AATAAA"
        motif_1_50_base = count_motif_in_sequence(seq_50_base, motif_1)
        motif_1_100_base = count_motif_in_sequence(seq_100_base, motif_1)

        motif_2 = "ATTAAA"
        motif_2_50_base = count_motif_in_sequence(seq_50_base, motif_2)
        motif_2_100_base = count_motif_in_sequence(seq_100_base, motif_2)

        # Add motif counts to their respective columns in the dataframe.
        polya_sites.loc[index,motif_1+"_50bp"] = motif_1_50_base
        polya_sites.loc[index,motif_2+"_50bp"] = motif_2_50_base

        polya_sites.loc[index,motif_1+"_100bp"] = motif_1_100_base
        polya_sites.loc[index,motif_2+"_100bp"] = motif_2_100_base
    
    # Change count columns from pandas default float to int
    for col in [motif_1+"_50bp", motif_1+"_100bp", motif_2+"_50bp", motif_2+"_100bp"]:
        polya_sites[col] = polya_sites[col].astype('int64')
    return polya_sites


args = get_args()

polya_sites = read_in_polya_site_bed_as_dataframe(args["input_file"])

genome_seq = read_in_genome_sequence(args["genome_file"])

polya_sites = find_poly_a_signals_for_each_site(polya_sites, genome_seq)

polya_sites.to_csv(args["output_file"], sep = '\t')


    

