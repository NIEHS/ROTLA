# ROTLA
ROTLA (Reader of the Lost Arcs) is a Python package that applies a split-read approach to detect deletions in mitochondrial genomes.

## Requirements
ROTLA was developed using Python 2.7.13. In addition to requirements specified in setup.py, ROTLA requires installation of the BLAT command-line alignment utility. BLAT binaries may be downloaded from the UCSC Genome Browser here:
* [UCSC Utilities Download Page](http://hgdownload.soe.ucsc.edu/downloads.html#source_downloads)

## Installation
The location of the BLAT executable must be specified prior to installation. To do this, manually edit the path in `paths.cfg` using a text editor.

After the path has been set, install ROTLA using:
```
python setup.py install
```
Users without administrative privileges may install ROTLA in their home directories by appending '--user' to the command above. The ROTLA package can then be executed as ${HOME}/.local/bin/ROTLA, or simply, ROTLA, provided the user adds the location ${HOME}/.local/bin to their PATH environment variable.

## Usage
Muliple functions are accessible using ROTLA's command line interface. General usage is as follows:
```
ROTLA COMMAND [OPTIONS] [ARGS]...
```
Available commands:
* [compile-breakpoint-results](#compile-breakpoint-results)
* [find-breakpoints](#find-breakpoints)
* [get-aligned-bases](#get-aligned-bases)

### compile-breakpoint-results
```
ROTLA compile-breakpoint-results [OPTIONS] LIST_FILE_NAME OUTPUT_FILE_NAME
```
Given a list of breakpoint files, create a composite table containing counts for all observed breakpoints in all files. The input list file must contain two tab-separated columns with no header line. Entries in column 1 should identify the name of a breakpoint file and entries in column 2 should specify the corresponding name to be written to the header line in the output file. See `example_list.txt` in the `docs` folder for an illustration of this format.

### find-breakpoints
```
ROTLA find-breakpoints [OPTIONS] READ_1_FASTQ_FILE READ_2_FASTQ_FILE REFERENCE_SEQUENCE OUTPUT_PREFIX
```
Given a set of paired-end FASTQ files and FASTA reference sequence, identify breakpoint coordinates and determine count of supporting reads.
This command will produce the following output files, with each name below preceded by the provided `OUTPUT_PREFIX`:

* `OUTPUT_PREFIX`.read_1.psl

Output of Read 1 FASTQ blat alignment in psl format

* `OUTPUT_PREFIX`.read_2.psl

Output of Read 2 FASTQ blat alignment in psl format

* `OUTPUT_PREFIX`.read_1.blat.out

Content written to STDOUT during Read 1 blat alignment

* `OUTPUT_PREFIX`.read_2.blat.out

Content written to STDOUT during Read 2 blat alignment

* `OUTPUT_PREFIX`.breakpoints.txt

Tab-delimited table of breakpoint start coordinates, end coordinates, and counts of supporting reads

#### Options
* `--length INTEGER`

Minimum required alignment length, default = 25

### get-aligned-bases
```
ROTLA get-aligned-bases [OPTIONS] INPUT_FILE_PREFIX REFERENCE_SEQUENCE
```
Given a pair of PSL files produced using find_breakpoints and the FASTA reference sequence, this command will determine the total count of aligned bases and print this value to an output file named `INPUT_PREFIX`.aligned_bases.txt. To allow aligned base counts of many samples to be easily combined, this output file utlizes a two-column tab-delimited format where the first contains the input file prefix and the second contains the count itself.

## Authors
ROTLA was conceptualized by Christopher Lavender and Scott Lujan. ROTLA was written by Christopher Lavender and Adam Burkholder.

## License
This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
