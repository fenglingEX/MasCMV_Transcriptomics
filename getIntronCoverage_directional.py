from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import biomodule
import os
import io
import sys
import subprocess
import re
import pandas as pd
import argparse
from argparse import ArgumentParser


def check_is_even(value):
    int_value = int(value)
    if (int_value % 2) !=0:
         raise argparse.ArgumentTypeError("%s is an invalid even integer value" % value)
    else:
        return int_value 


def get_args():
    parser = ArgumentParser()
    parser.add_argument("-r", "--reference", type =str, required = True, help = "The reference genome mapped onto by STAR, in fasta format.")
    parser.add_argument("-i", "--input_file", type = str, required = True, help = "Tab separated splice junction file output by STAR.")
    parser.add_argument("-o", "--output_file", type = str, required = True, help = "Output .tsv file to store results in.")
    parser.add_argument("-s", "--transcript_sense_reads", type = str, required = True, help = "Paired reads file in fastq format for reads in the same orientation as the splice junctions.")
    parser.add_argument("-a", "--transcript_antisense_reads", type = str, required = True, help = "Paired reads file in fastq format for reads in the opposite orientation to the splice junctions.")
    parser.add_argument("-t", "--threads", type = int, default = 1, help = "Number of threads to allocate to Jellyfish. Default: 1")
    parser.add_argument("-m", "--hash_size", type = int, default = 10, help = "Jellyfish hash size in millions of kmers. Larger hash size speeds up jellyfish but uses more memory. Default: 10")
    parser.add_argument("-k", "--kmer_length", type = check_is_even, default = 20, help = "Length of splice sequences to count. Used by jellyfish to count kmers. Must be even. Default: 20")
    args = vars(parser.parse_args())
    return args


# Takes splice junction file ouput by STAR and uses it to find intron-related sequences in the reference genome.
# Stores sequences in a DataFrame; each row corresponds to a single splice junction detected by STAR
def construct_splice_site_dataframe(reference, input_file, kmer_length):

    # Create dictionary to store intron-related sequences as parallel lists.
    sequence_keys = ["Junction", "Intron5end", "Intron3end"]
    all_keys = ["Chromosome", "Intron First Base", "Intron Last Base", "Strand", "Intron Motif", "Maximum Overhang"] + sequence_keys
    value_lists = []
    for key in all_keys:
        value_lists.append([])
    sequences_dict = dict(zip(all_keys, value_lists))
    
    # Read STAR output file and extract each junction's related sequences.
    while True:
        line = input_file.readline().rstrip() # remove trailing newline
        if not line:
            break 
        fields = line.split("\t")

        # Store STAR output in dict
        sequences_dict["Chromosome"].append(fields[0])
        sequences_dict["Intron First Base"].append(fields[1])
        sequences_dict["Intron Last Base"].append(fields[2])

        # Convert STAR code for strand to strand symbol.
        strand_codes = {"0": "Undefined", "1": "+", "2": "-"}
        sequences_dict["Strand"].append(strand_codes[fields[3]])

        # Convert STAR code for intron motif to intron motif itself.
        intron_motif_codes = {"0": "non-canonical", "1": "GT/AG", "2": "CT/AC", "3": "GC/AG", "4":"CT/GC", "5":"AT/AC", "6": "GT/AT"}
        sequences_dict["Intron Motif"].append(intron_motif_codes[fields[4]])

        sequences_dict["Maximum Overhang"].append(fields[8])

        # Extract splice site sequences from reference genome and store them in dict
        # Get 0-based intron positions.
        intron_start = int(fields[1])-1
        intron_end = int(fields[2])
        nucleotides_to_extract = kmer_length
        junction_nucleotides_to_extract = int(kmer_length/2)

    
        # Extract splice junction sequence (junction_nucleotides_to_extract bases on each side)
        junction = str.upper(reference[intron_start - junction_nucleotides_to_extract:intron_start] + reference[intron_end:intron_end + junction_nucleotides_to_extract])
    
        # Extract 5' intron sequences (first nucleotides_to_extract bases of intron)
        intron5end = str.upper(reference[intron_start:intron_start + nucleotides_to_extract])

        # Extract 3' intron sequences (last nucleotides_to_extract bases of intron)
        intron3end = str.upper(reference[intron_end - nucleotides_to_extract:intron_end])

        # Add sequences to their respective lists.
        sequences_dict["Junction"].append(junction)
        sequences_dict["Intron5end"].append(intron5end)
        sequences_dict["Intron3end"].append(intron3end)

    # Create column headings for counts of each sequence.
    # Zip arranges the columns so that corresponding columns are right next to each other in the dataframe.
    expected_keys = zip([x+"_Expected" for x in sequence_keys], [x+"_Unexpected" for x in sequence_keys])
    read_keys = zip([x+"_sense_reads" for x in sequence_keys], [x+"_antisense_reads" for x in sequence_keys])
    count_keys = []
    for key_pair in expected_keys:
        count_keys = count_keys + list(key_pair)
    for key_pair in read_keys:
        count_keys = count_keys + list(key_pair)

    # Create column heading to store non-unique kmers 
    non_unique_key = ["Non-Unique Kmers"]

    # Use dictionary of lists to create dataframe. Count columns are NaN
    splice_sequences = pd.DataFrame(data = sequences_dict, columns = all_keys+count_keys+non_unique_key)

    # Replace count column NaNs with 0
    for key in count_keys:
        splice_sequences[key] = 0

    # Replace non-unique column NaNs with empty string
    splice_sequences[non_unique_key] = ""

    return splice_sequences


def get_jellyfish_index_filename(sequences_file, kmer_length):
    base_filename = os.path.basename(sequences_file)

    index_filename = re.split("\.", base_filename)[0] + "_" + str(kmer_length) + "mers.jf"

    directory = "Jellyfish/"

    # Create directory to store jellyfish index if it does not exist.
    if not os.path.isdir(directory):
        os.mkdir(directory)

    filepath = os.path.join(directory, index_filename)

    return filepath


# Creates jellyfish kmer index for a fasta or fastq file.
def create_jellyfish_index_for_file(sequences_file, threads, jellyfish_index_filename, hash_size, kmer_length):
    # cat genome file
    cat_genome = subprocess.Popen(["cat", sequences_file], shell = False, stdin = subprocess.PIPE, stdout = subprocess.PIPE)

    # Pass genome file to jellyfish count via standard input to create kmer index.
    jellyfish_count = subprocess.Popen(["jellyfish", "count", "-t", str(threads), "-m", str(kmer_length), "-s", str(hash_size)+"M", "-o", jellyfish_index_filename, "/dev/fd/0"], shell = False, stdin = cat_genome.stdout)

    jellyfish_count.communicate()
    return 


# Returns a list of kmers that appear more than once in the reference genome.
# Used as a blacklist when assessing splice junction kmer frequency.
def create_reference_genome_kmer_blacklist(jellyfish_reference_index_filename):
    # Retrieve list of all kmers found in reference genome and their count. One kmer per line in format "KMER COUNT"
    jellyfish_dump = subprocess.Popen(["jellyfish", "dump", "-c", jellyfish_reference_index_filename], shell = False, stdin = subprocess.PIPE, stdout = subprocess.PIPE)

    # Grab results and decode from binary
    reference_kmer_counts = jellyfish_dump.communicate()[0].decode().rstrip() # Removes last newline

    # Turn into list of strings, each element becomes "KMER COUNT"
    reference_kmer_counts = re.split("\n", reference_kmer_counts)

    non_unique_kmers = []
    # Add any kmers with count > 1 to list
    for row in reference_kmer_counts:
        kmer = re.split(" ", row) 
        sequence = kmer[0]
        count = int(kmer[1])
        if count > 1:
            non_unique_kmers.append(sequence)

    return non_unique_kmers


def check_for_blacklisted_kmers(sequences_df, kmer_blacklist):
    # If any of the splice-site sequences are in the blacklist, add them to the non-unique column
    for col in ["Junction", "Intron5end", "Intron3end"]:
        # Get boolean mask of all sequences in blacklist for this column
        isin_mask = sequences_df[col].isin(kmer_blacklist)
        # Add all non-unique sequences to the non-unique column of their respective row
        sequences_df.loc[isin_mask, "Non-Unique Kmers"] = sequences_df.loc[isin_mask, "Non-Unique Kmers"].str.cat(others = sequences_df[col][isin_mask], sep = ',')
        
    unique_kmers_mask = sequences_df.loc[:, "Non-Unique Kmers"] == ""
    sequences_df.loc[unique_kmers_mask, "Non-Unique Kmers"] = "0"
    return sequences_df


# Takes a list of sequences and the filepath to a jellyfish kmer index
# Calls jellyfish query to count occurences of each kmer.
# Returns count for each input sequence in a dictionary{k:v} where k=sequence and v=count
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
def update_counts(sequences_df, counts_dict, col, reads_file):
    count_col = col+"_Total"
    fq_col = col+"_"+reads_file

    for sequence in counts_dict:
        # Get indices of sequence in datafrane as boolean mask
        sequence_mask = sequences_df[col].isin([sequence]) 
        # Update count for sequence
        sequences_df.loc[sequence_mask, count_col] += counts_dict[sequence]
        sequences_df.loc[sequence_mask, fq_col] += counts_dict[sequence]
    return sequences_df


def update_expected_counts(sequences_df, counts_dict, col, reads_file):
    expected_count_col = col+"_Expected"
    fq_col = col+"_"+reads_file

    for sequence in counts_dict:
            # Get indices of sequence in datafrane as boolean mask
            sequence_mask = sequences_df[col].isin([sequence]) 
            # Update count for sequence
            sequences_df.loc[sequence_mask, expected_count_col] += counts_dict[sequence]
            sequences_df.loc[sequence_mask, fq_col] += counts_dict[sequence]
    return sequences_df


def update_unexpected_counts(sequences_df, counts_dict, col, reads_file):
    unexpected_count_col = col+"_Unexpected"
    fq_col = col+"_"+reads_file

    for sequence in counts_dict:
            # Get indices of sequence in datafrane as boolean mask
            sequence_mask = sequences_df[col].isin([sequence]) 
            # Update count for sequence
            sequences_df.loc[sequence_mask, unexpected_count_col] += counts_dict[sequence]
            sequences_df.loc[sequence_mask, fq_col] += counts_dict[sequence]
    return sequences_df


# Utility function to reverse complement pd.Series of DNA sequences
def reverse_complement_series_of_sequences(sequences_series):
    rc_sequences = [str(Seq(sequence, IUPAC.unambiguous_dna).reverse_complement()) for sequence in sequences_series]
    rc_sequences = pd.Series(rc_sequences)
    return rc_sequences


# Get counts for each sense (plus strand) splice junction's sequences
# Takes a dataframe containing sequences and the count of their occurrences 
def get_counts_for_sense_splice_sites(sequences_df, reads_jellyfish_index, reads_are_sense):
    # For each splice region (junction, 5' end of intron and 3' end of intron)
    # for col in sequences_df:
    sense_splice_site_mask = sequences_df.loc[:,"Strand"] == '+'
    sense_splice_site_df =sequences_df.loc[sense_splice_site_mask]
    for col in ["Junction", "Intron5end", "Intron3end"]:
        sequences = sense_splice_site_df[col]

        rc_sequences = reverse_complement_series_of_sequences(sequences)

        # Query jellyfish index for plus strand counts
        seq_occurence_counts = query_jellyfish_index(sequences, reads_jellyfish_index, sequences)

        # Query jellyfish index for minus strand counts
        rc_seq_occurence_counts = query_jellyfish_index(rc_sequences, reads_jellyfish_index, sequences)

        if reads_are_sense:
            sequences_df = update_expected_counts(sequences_df, seq_occurence_counts, col, "sense_reads")

            sequences_df = update_unexpected_counts(sequences_df, rc_seq_occurence_counts, col, "sense_reads")
        elif not reads_are_sense:
            sequence_df = update_expected_counts(sequences_df, rc_seq_occurence_counts, col, "antisense_reads")

            sequences_df = update_unexpected_counts(sequences_df, seq_occurence_counts, col, "antisense_reads")
    return sequences_df 


def get_counts_for_antisense_splice_sites(sequences_df, reads_jellyfish_index, reads_are_sense):
    # For each splice region (junction, 5' end of intron and 3' end of intron)
    # for col in sequences_df:
    antisense_splice_site_mask = sequences_df.loc[:,"Strand"] == '-'
    antisense_splice_site_df =sequences_df.loc[antisense_splice_site_mask]
    for col in ["Junction", "Intron5end", "Intron3end"]:
        sequences = antisense_splice_site_df[col]

        rc_sequences = reverse_complement_series_of_sequences(sequences)

        # Query jellyfish index for plus strand counts
        seq_occurence_counts = query_jellyfish_index(sequences, reads_jellyfish_index, sequences)

        # Query jellyfish index for minus strand counts
        rc_seq_occurence_counts = query_jellyfish_index(rc_sequences, reads_jellyfish_index, sequences)

        if reads_are_sense:
            sequences_df = update_expected_counts(sequences_df, rc_seq_occurence_counts, col, "sense_reads")

            sequences_df = update_unexpected_counts(sequences_df, seq_occurence_counts, col, "sense_reads")
        elif not reads_are_sense:
            sequences_df = update_expected_counts(sequences_df, seq_occurence_counts, col,  "antisense_reads")

            sequences_df = update_unexpected_counts(sequences_df, rc_seq_occurence_counts, col,  "antisense_reads")    
    return sequences_df        
   


def main():
    args = get_args()
    
    # Read in reference genome.
    reference = str(SeqIO.read(args["reference"], "fasta").seq)
    
    # Open STAR output file for reading splice junction info
    infile = open(args["input_file"])

    # Extract splice junction sequences from reference genome using STAR output
    splice_sequences = construct_splice_site_dataframe(reference, infile, args["kmer_length"])

    # Create filename to store Jellyfish index under Jellyfish directory for first fastq.
    jellyfish_transcript_sense_reads_index_filename = get_jellyfish_index_filename(args["transcript_sense_reads"], args["kmer_length"])
    # Create jellyfish kmer index for first fastq.
    create_jellyfish_index_for_file(args["transcript_sense_reads"], args["threads"], jellyfish_transcript_sense_reads_index_filename, args["hash_size"], args["kmer_length"])

    # Create filename to store Jellyfish index under Jellyfish directory for seconds fastq.
    jellyfish_transcript_antisense_reads_index_filename = get_jellyfish_index_filename(args["transcript_antisense_reads"], args["kmer_length"])
    # Create jellyfish kmer index for second fastq.
    create_jellyfish_index_for_file(args["transcript_antisense_reads"], args["threads"], jellyfish_transcript_antisense_reads_index_filename, args["hash_size"], args["kmer_length"])

    # Create filename to store Jellyfish index for the reference index, under Jellyfish directory
    jellyfish_reference_index_filename = get_jellyfish_index_filename(args["reference"], args["kmer_length"])
    # Create jellyfish kmer index for reference genome.
    create_jellyfish_index_for_file(args["reference"], args["threads"], jellyfish_reference_index_filename, args["hash_size"], args["kmer_length"])

    # Create kmer blacklist: list of all kmers that appear more than once in reference genome
    kmer_blacklist = create_reference_genome_kmer_blacklist(jellyfish_reference_index_filename)

    splice_sequences = check_for_blacklisted_kmers(splice_sequences, kmer_blacklist)

    for col in splice_sequences:
        isin_mask = splice_sequences[col].isin(kmer_blacklist)

    splice_sequences = get_counts_for_sense_splice_sites(splice_sequences, jellyfish_transcript_sense_reads_index_filename, True)
    splice_sequences = get_counts_for_sense_splice_sites(splice_sequences, jellyfish_transcript_antisense_reads_index_filename, False)

    splice_sequences = get_counts_for_antisense_splice_sites(splice_sequences, jellyfish_transcript_sense_reads_index_filename, True)
    splice_sequences = get_counts_for_antisense_splice_sites(splice_sequences, jellyfish_transcript_antisense_reads_index_filename, False)

    # # Query counts of reference genome splice junction sequences in Jellyfish index and update counts in dataframe.
    # for i, reads_index in enumerate([jellyfish_transcript_sense_reads_index_filename, jellyfish_transcript_antisense_reads_index_filename], 1):
    #   splice_sequences = get_counts_for_splice_sites(splice_sequences, reads_index, i) # Add extra argument which specifies which fastq column to add numbers for this index to. Still sum total as normal.

    # TODO: Compare splice sequences to blacklist.
    # TODO: Finalise output format
    
    splice_sequences.to_csv(args["output_file"], sep = "\t")


if __name__ == '__main__':
    main()