from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import pandas as pd 
import sys
import re
import subprocess
import io
import os
from argparse import ArgumentParser


def get_args():
    parser = ArgumentParser()
    parser.add_argument("-r", "--reference", type =str, required = True, help = "The reference genome in .fasta format. Same genome used by ContextMap2")
    parser.add_argument("-i", "--input_file", type = str, required = True, help = "Tab separated polyA site .bed file output by ContextMap2. Or extended version of said file.")
    parser.add_argument("-o", "--output_file", type = str, required = True, help = "Output .tsv file to store results in.")
    parser.add_argument("-1", "--reads_1", type = str, required = True, help = "First paired-end reads file.")
    parser.add_argument("-2", "--reads_2", type = str, required = True, help = "Second paired-end reads file.")
    parser.add_argument("-t", "--threads", type = int, default = 1, help = "Number of threads to allocate to Jellyfish. Default: 1")
    parser.add_argument("-m", "--hash_size", type = int, default = 10, help = "Jellyfish hash size in millions of kmers. Larger hash size speeds up jellyfish but uses more memory. Default: 10")
    parser.add_argument("-l", "--transcript_length", type = int, default = 20, help = "Number of bases before a polyA site to extract when searching for polyA reads. Used by jellyfish. Default: 20")
    parser.add_argument("-a", "--polya_length", type = int, default = 6, help = "The number of As to add to the ends of transcripts when searching for poly A reads. Used by jellyfish. Default: 6")
    args = vars(parser.parse_args())
    return args


def construct_polyA_dataframe(reference, input_file, transcript_length, polya_length):
    ref_genome = SeqIO.read(reference, "fasta")

    col_names = ["Name","Chromosome", "Start", "End", "Score", "Strand"]
    polya_df = pd.read_csv(input_file, sep = '\t',  usecols = [0,1,2,3,4,5], index_col = 0, names = col_names, header = 0)

    polya_df.Start = polya_df.Start.astype('int64')
    polya_df.End = polya_df.End.astype('int64')
    polya_df.Score = polya_df.Score.astype('int64')

    for index, site in polya_df.iterrows():
        start_pos = polya_df.loc[index, "Start"]
        end_pos = polya_df.loc[index, "End"]
        pre_poly_a_seq = str(ref_genome.seq[start_pos-transcript_length:start_pos])
        if polya_df.loc[index, "Strand"] == '+':
            pre_poly_a_seq = str(ref_genome.seq[start_pos-transcript_length:start_pos])
        elif polya_df.loc[index, "Strand"] == '-':
            print("Antisense transcript")
            pre_poly_a_seq = ref_genome.seq[end_pos:end_pos+transcript_length]
            print("Sense: " + str(pre_poly_a_seq))
            pre_poly_a_seq = str(pre_poly_a_seq.reverse_complement())
            print("Antisense: "+ str(pre_poly_a_seq))

        a_tail = "A"*polya_length
        transcript_end = pre_poly_a_seq+a_tail

        polya_df.loc[index, "Sequence"] = transcript_end

    polya_df.loc[:,"Reads"] = 0

    return polya_df


def get_jellyfish_reads_index_filename(reads_1, reads_2):
  reads_1 = os.path.basename(reads_1)
  reads_2 = os.path.basename(reads_2)

  prefix = os.path.commonprefix([reads_1, reads_2])
  index_filename = prefix + "_polya_kmers.jf"
  directory = "Jellyfish/"

  # Create directory to store jellyfish index if it does not exist.
  if not os.path.isdir(directory):
      os.mkdir(directory)

  filepath = os.path.join(directory, index_filename)

  return filepath


# Creates jellyfish kmer index for a fasta or fastq file.
def create_jellyfish_index_for_reads(reads_1, reads_2, threads, jellyfish_index_filename, hash_size, kmer_length):
    # cat genome file
    cat_genome = subprocess.Popen(["cat", reads_1, reads_2], shell = False, stdin = subprocess.PIPE, stdout = subprocess.PIPE)

    # Pass genome file to jellyfish count via standard input to create kmer index.
    jellyfish_count = subprocess.Popen(["jellyfish", "count", "-t", str(threads), "-m", str(kmer_length), "-C", "-s", str(hash_size)+"M", "-o", jellyfish_index_filename, "/dev/fd/0"], shell = False, stdin = cat_genome.stdout)

    jellyfish_count.communicate()
    return 


def query_jellyfish_index(sequence_list, jellyfish_index, reference_sequence_list):
    # Create jellyfish query command which will query all sequences provided.
    query_cmd = ["jellyfish", "query", jellyfish_index] + sequence_list.values.tolist() # Convert from pd.Series to list for concatenation

    jellyfish_query = subprocess.Popen(query_cmd, shell=False, stdout = subprocess.PIPE)

    # Grab query results
    query_results = jellyfish_query.communicate()[0] # [0] takes only the STDOUT output

    # decode from binary
    query_results = query_results.decode().strip() # Removes last newline

    # Turn into list of strings, each element becomes "KMER COUNT"
    query_results = re.split("\n", query_results)

    # Revert reverse complementation of kmers done by Jellyfish
    # TODO: Change this to work with a separate list of accepted sequences instead of the one being operated on
    query_results = reverse_complement_kmers_modified_by_jellyfish(query_results, reference_sequence_list)

    # Turn results into dict, where key is the sequence and value is the count
    query_results_dict = {re.split(" ", kmer)[0]:int(re.split(" ", kmer)[1]) for kmer in query_results}


    return query_results_dict


# Checks list of kmers return by Jellyfish query
# Any kmers that don't match sequences extracted from the reference genome are
# reverse complemented. This reverts reverse complementation performed by the -C option in Jellyfish
def reverse_complement_kmers_modified_by_jellyfish(jellyfish_query_results, sequence_list):
    # Turn genome sequences into Seq objects for comparison
    sequence_list = [Seq(seq, IUPAC.unambiguous_dna) for seq in sequence_list]

    corrected_query_results = []

    for kmer in jellyfish_query_results:
        # Split kmer output into sequence and count
        kmer = re.split(" ", kmer) 
        sequence = kmer[0]
        count = kmer[1]

        # If kmer does not match any sequence extracted from reference, reverse complement
        if sequence not in sequence_list:
            sequence = str(Seq(sequence, IUPAC.unambiguous_dna).reverse_complement())
            if sequence not in sequence_list:
                print("Error, neither the kmer nor its reverse complement is found in the reference genome:")
                print(sequence)

        # Format sequence and count back into a single string and add to list.
        corrected_query_results.append(sequence + " " + count)

    return corrected_query_results


# Takes a dataframe containing splice junction sequences and their counts
# Takes a dictionary {k:v} where k=sequence and v=count
# Updates the counts for the splice junction sequences under column col
def update_counts(sequences_df, counts_dict):
    count_col = "Reads"

    for sequence in counts_dict:
        # Get indices of sequence in datafrane as boolean mask
        sequence_mask = sequences_df["Sequence"].isin([sequence]) 

        # Update count for sequence
        sequences_df.loc[sequence_mask, count_col] += counts_dict[sequence]
    return sequences_df


# Get counts for each splice junction's sequences
# Takes a dataframe containing sequences and the count of their occurrences 
def get_counts_for_polya_sites(polya_df, jellyfish_index):
    sequences = polya_df["Sequence"]

    # Query jellyfish index for counts
    seq_occurence_counts = query_jellyfish_index(sequences, jellyfish_index, sequences)

    # Use query results to update dataframe containing sequence counts
    polya_df = update_counts(polya_df, seq_occurence_counts)

    return polya_df


def main():
    args = get_args()

    # Extract splice junction sequences from reference genome using STAR output
    polya_df = construct_polyA_dataframe(args["reference"], args["input_file"], args["transcript_length"], args["polya_length"])

    jellyfish_kmer_length = args["transcript_length"] + args["polya_length"]

    # Create filename to store Jellyfish index under Jellyfish directory for first fastq.
    jellyfish_reads_index_filename = get_jellyfish_reads_index_filename(args["reads_1"], args["reads_2"])

    # Create jellyfish kmer index for first fastq.
    create_jellyfish_index_for_reads(args["reads_1"], args["reads_2"], args["threads"], jellyfish_reads_index_filename, args["hash_size"], jellyfish_kmer_length)

    polya_df = get_counts_for_polya_sites(polya_df, jellyfish_reads_index_filename)

    polya_df.to_csv(args["output_file"], sep = "\t")


if __name__ == '__main__':
    main()