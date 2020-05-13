# -*- coding: utf-8 -*-

"""Console script for ROTLA."""

import click
import os

from ROTLA import ROTLA as _find_breakpoints
from compile_breakpoint_results import compile_breakpoints as _compile_breakpoints
from aligned_bases_from_psl import get_aligned_bases as _get_aligned_bases

@click.group()
def main(args=None):
    pass

@main.command()
@click.option('--length', type=int, help='Minimum required alignment length, default = 25',
              default=25)
@click.argument('read_1_fastq_file', type=str)
@click.argument('read_2_fastq_file', type=str)
@click.argument('reference_sequence', type=str)
@click.argument('output_prefix', type=str)
def find_breakpoints(read_1_fastq_file, read_2_fastq_file, reference_sequence,
                     output_prefix, length):
    '''
    Identify mitochondrial breakpoints.

    Given a set of paired-end FASTQ files and FASTA reference sequence,
    identify breakpoint coordinates and determine count of supporting
    reads.

    This command will produce the following output files, with each
    name below preceded by the provided output_prefix:

    .read_1.psl         Output of Read 1 FASTQ blat alignment in psl format
    .read_2.psl         Output of Read 2 FASTQ blat alignment in psl format
    .read_1.blat.out    Content written to STDOUT during Read 1 blat alignment
    .read_2.blat.out    Content written to STDOUT during Read 2 blat alignment
    .breakpoints.txt    Tab-delimited table of breakpoint start, end, counts
    '''
    args = { 'read_1_file_name':read_1_fastq_file,
             'read_2_file_name':read_2_fastq_file,
             'reference_sequence':reference_sequence,
             'output_prefix':output_prefix,
             'length':length }
    _find_breakpoints(**args)

@main.command()
@click.argument('list_file_name', type=str)
@click.argument('output_file_name', type=str)
def compile_breakpoint_results(list_file_name, output_file_name):

    '''
    Combine results from multiple samples.

    Given a list of breakpoint files created using find_breakpoints,
    create a composite table containing counts for all observed
    breakpoints in all files. The input list file must contain
    two tab-separated columns with no header line. Entries in column 1
    should identify the name of a breakpoint file and entries in
    column 2 should specify the corresponding name to be written to
    the header line in the output file.
    '''

    _compile_breakpoints(list_file_name, output_file_name)

@main.command()
@click.argument('input_file_prefix', type=str)
@click.argument('reference_sequence', type=str)
def get_aligned_bases(input_file_prefix, reference_sequence):

    '''
    Count bases aligned by find_breakpoints.

    Given a pair of PSL files produced using find_breakpoints and the FASTA
    reference sequenced used, this command will determine the total count
    of aligned bases and print this value to an output file named
    [input_prefix].aligned_bases.txt. To allow aligned base counts of many
    samples to be easily combined, this output file utlizes a two-column
    tab-delimited format where the first contains the input file prefix and
    the second contains the count itself.
    '''

    _get_aligned_bases(input_file_prefix, reference_sequence)

if __name__ == "__main__":
    main()
