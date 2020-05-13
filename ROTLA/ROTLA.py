#!/usr/bin/env python

import sys
import os
import argparse
import copy

from subprocess import call
from collections import defaultdict
from __init__ import PATHS

class ROTLA(object):
    
    def __init__(self, **kwargs):
        
        # Set instance variables
        self.read_1_fn = kwargs['read_1_file_name']
        self.read_2_fn = kwargs['read_2_file_name']
        self.ref_fn = kwargs['reference_sequence']
        self.output_header = kwargs['output_prefix']
        self.blat_path = PATHS['blat']
        self.required_alignment_length = kwargs['length']
        
        # File checks
        for fn in [
            self.read_1_fn,
            self.read_2_fn,
            self.ref_fn,
            self.blat_path
        ]:
            os.path.exists(fn)
        
        self.alignment = dict()
        self.breakpoints = dict()
        self.break_count = defaultdict(int)
        
        self.execute()
        
    @staticmethod
    def readReference(fasta_file):
        reference = ""
        header_count = 0
        
        with open(fasta_file) as f:
            for line in f:
                if line[0] != ">":
                    reference += line.strip()
                else:
                    header_count += 1
                
                if header_count > 1:
                    raise StandardError('FASTA contains two headers.')
        
        return reference.upper()

    @staticmethod
    def makeFASTA(fastq_file, fasta_file):
        count = 0
        with open(fastq_file) as fastq, open(fasta_file, "w") as fasta:
            for line in fastq:
                if count % 4 == 0:
                    fasta.write("> " + line.strip() + "\n")
                if count % 4 == 1:
                    fasta.write(line)
                count += 1
    
    @staticmethod
    def makePaddedFASTA(fasta_file, output_file):
        reference_sequence = ROTLA.readReference(fasta_file)
        
        with open(output_file, "w") as OUTPUT:
            OUTPUT.write("> Padded reference\n")
            OUTPUT.write(reference_sequence + reference_sequence + "\n")

    def cleanFASTA(self):
        os.remove(self.output_header + ".padded_reference.fasta")
        os.remove(self.output_header + ".read_1.fasta")
        os.remove(self.output_header + ".read_2.fasta")
    
    def readAlignments(self, input_read_1_file, input_read_2_file):
        
        def readFile(input_file, read, add_if_present=False):
            with open(input_file) as f:
                
                # Go through BLAT header
                for i in range(5):
                    next(f)
                    
                # Iterate through lines
                for line in f:
    
                    [
                        matches,
                        misMatches,
                        repMatches,
                        nCount,
                        qNumInsert,
                        qBaseInsert,
                        tNumInsert,
                        tBaseInsert,
                        strand,
                        qName,
                        qSize,
                        qStart,
                        qEnd,
                        tName,
                        tSize,
                        tStart,
                        tEnd,
                        blockCount,
                        blockSizes,
                        qStarts,
                        tStarts,
                    ] = line.strip().split()
                    
                    qStarts = qStarts.split(",")[:-1]
                    tStarts = tStarts.split(",")[:-1]
                    blockSizes = blockSizes.split(",")[:-1]
                    
                    q_list = []
                    t_list = []
                    
                    if (len(blockSizes) > 1) or (add_if_present and qName in self.alignment):
                        
                        if qName not in self.alignment:
                            self.alignment[qName] = {'read_1': [], 'read_2': []}
                        
                        for q, t, size in zip(qStarts, tStarts, blockSizes):
                            if strand == "+":
                                q = [
                                    int(q) + 1,
                                    int(q) + int(size),
                                ]
                            if strand == "-":
                                q = [
                                    int(qSize) - int(q) - int(size) + 1,
                                    int(qSize) - int(q),
                                ]
                            
                            t = [
                                int(t) + 1,
                                int(t) + int(size),
                            ]
        
                            for i, value in enumerate(t):
                                if value > self.ref_seq_length:
                                    t[i] = value - self.ref_seq_length
                            
                            q_list.append(q)
                            t_list.append(t)
                        
                        mapped_regions = []
                        for q, t in zip(q_list, t_list):
                            mapped_regions.append({
                                'q': q,
                                't': t,
                            })
                        
                        read_dict = {
                            'mapped_regions': mapped_regions,
                            'strand': strand,
                            'matches': int(matches),
                            'tGap': int(tNumInsert),
                            'qGap': int(qNumInsert),
                        }
                        
                        if read_dict not in self.alignment[qName][read]:
                            self.alignment[qName][read].append({
                                'mapped_regions': mapped_regions,
                                'strand': strand,
                                'matches': int(matches),
                                'tGap': int(tNumInsert),
                                'qGap': int(qNumInsert),
                            })

        readFile(input_read_1_file, 'read_1')
        readFile(input_read_2_file, 'read_2')
        readFile(input_read_1_file, 'read_1', add_if_present=True)
        readFile(input_read_2_file, 'read_2', add_if_present=True)

    def findBreaks(self):
        for query in self.alignment:
            for read, alignments in self.alignment[query].items():

                for alignment in alignments:
                    
                    strand = alignment['strand']
                
                    if strand == '+':
                        sorted_alignments = sorted(alignment['mapped_regions'], key=lambda k: k['q'][0])
                    if strand == '-':
                        sorted_alignments = sorted(alignment['mapped_regions'], key=lambda k: k['q'][0], reverse = True)
                        
                    for i in range(len(sorted_alignments)-1):
                        
                        size_1 = abs(sorted_alignments[i]['q'][1] - sorted_alignments[i]['q'][0]) + 1
                        size_2 = abs(sorted_alignments[i+1]['q'][1] - sorted_alignments[i+1]['q'][0]) + 1
                        
                        if size_1 >= self.required_alignment_length and \
                                size_2 >= self.required_alignment_length:
                            
                            if query not in self.breakpoints:
                                self.breakpoints[query] = {'read_1': set(), 'read_2': set()}
                            
                            self.breakpoints[query][read].add((sorted_alignments[i]['t'][1], sorted_alignments[i+1]['t'][0]))
                            
    def compareBreaksAcrossReads(self, ref_seq):
        
        def leftAlign(breakpoints):
            output_breakpoints = {'read_1': set(), 'read_2': set()}
            padded_seq = ref_seq + ref_seq
            
            for read, breakpoint_list in breakpoints.items():

                for breakpoint in breakpoint_list:
                    breakpoint = list(breakpoint)
                    
                    if breakpoint[1] < breakpoint[0]:
                        breakpoint[1] += self.ref_seq_length
                    
                    repeat = True
                    while repeat:
                        repeat = False
                        
                        del_seq = padded_seq[breakpoint[0]:breakpoint[1]-1]
                        
                        index = 1
                        while index <= len(del_seq):
                            if del_seq[-index:] == padded_seq[breakpoint[0]-index:breakpoint[0]]:
                                breakpoint[0] -= index
                                breakpoint[1] -= index
                                repeat = True
                                break
                            index += 1
                    
                    if breakpoint[1] > self.ref_seq_length:
                        breakpoint[1] -= self.ref_seq_length

                    output_breakpoints[read].add(tuple(breakpoint))
            
            return output_breakpoints
        
        def findNestedDeletions(breakpoints):
            output_breakpoints = {'read_1': set(), 'read_2': set()}
            
            breaklist = []
            for read, breakpoint_list in breakpoints.items():
                for breakpoint in breakpoint_list:
                    breaklist.append([read, breakpoint[0], breakpoint[1]])

            repeat = True
            while repeat:
                repeat = False
                aligned_breaklist = copy.copy(breaklist)
                
                for i, break_1 in enumerate(breaklist):
                    for j, break_2 in enumerate(breaklist):
                        if break_1 != break_2 and break_1[1] >= break_2[1] and break_1[2] <= break_2[2] and (break_1[1] != break_2[1] or break_1[2] != break_2[2]):
                            aligned_breaklist[j] = [break_2[0], break_1[1], break_1[2]]
                            repeat = True
                
                breaklist = aligned_breaklist
            
            for breakpoint in breaklist:
                output_breakpoints[breakpoint[0]].add((breakpoint[1], breakpoint[2]))
            
            return output_breakpoints
        
        def removeConflictingDeletions(query, breakpoints):
            
            def overlapsOther(breakpoint, query, break_read):
                if break_read == 'read_1':
                    read = 'read_2'
                if break_read == 'read_2':
                    read = 'read_1'
                
                for alignment in self.alignment[query][read]:
                    strand = alignment['strand']
                    mapped_regions = alignment['mapped_regions']
                    
                    if strand == '+':
                        sorted_regions = sorted(mapped_regions, key=lambda k: k['q'][0])
                    if strand == '-':
                        sorted_regions = sorted(mapped_regions, key=lambda k: k['q'][0], reverse = True)
                    
                    t_range = [(sorted_regions[0]['t'][0], sorted_regions[-1]['t'][1])]

                    if t_range[0][0] > t_range[0][1]:
                        t_range = [(1, t_range[0][1]), (t_range[0][0], self.ref_seq_length)]
                    
                    for _range in t_range:
                        if breakpoint[0] >= _range[0] and breakpoint[0] <= _range[1] and \
                            breakpoint[1] >= _range[0] and breakpoint[1] <= _range[1]:
                                return True
                
                return False
            
            output_breakpoints = {'read_1': set(), 'read_2': set()}
            
            breaklist = []
            for read, breakpoint_list in breakpoints.items():
                for breakpoint in breakpoint_list:
                    breaklist.append([read, breakpoint[0], breakpoint[1]])
            cleaned_breaklist = copy.copy(breaklist)
            
            for break_1 in breaklist:
                if overlapsOther((break_1[1], break_1[2]), query, break_1[0]):
                    other_read = False
                    for break_2 in breaklist:
                        if break_1[0] != break_2[0] and break_1[1] == break_2[1] and break_1[2] == break_2[2]:
                            other_read = True
                    if not other_read:
                        cleaned_breaklist.remove(break_1)
            
            for breakpoint in cleaned_breaklist:
                output_breakpoints[breakpoint[0]].add((breakpoint[1], breakpoint[2]))
            
            return output_breakpoints
        
        for query in self.breakpoints:
            breakpoints = copy.copy(self.breakpoints[query])
            
            breakpoints = leftAlign(breakpoints)
            breakpoints = findNestedDeletions(breakpoints)
            breakpoints = removeConflictingDeletions(query, breakpoints)
            
            self.breakpoints[query] = breakpoints
    
    def compileBreaks(self):
        
        for query in self.breakpoints:
            break_set = set()
            
            for breakpoint_list in self.breakpoints[query].values():
                for breakpoint in breakpoint_list:
                    break_set.add(tuple(breakpoint))
        
            for breakpoint in break_set:
                self.break_count[breakpoint] += 1
    
    def compareAcrossAllBreaks(self, ref_seq):
        
        repeat = True
        while repeat:
            repeat = False
            
            for break_1 in self.break_count.keys():
                for break_2 in self.break_count.keys():
                    if break_1 != break_2 and break_1[0] < break_2[0] and break_1[1] > break_2[0] and break_1[1] < break_2[1]:
                        offset_1 = ref_seq[break_1[0]:break_2[0]]
                        offset_2 = ref_seq[break_1[1]-1:break_2[1]-1]
                        if offset_1 == offset_2:
                            self.break_count[break_1] += self.break_count[break_2]
                            self.break_count.pop(break_2, None)

                            repeat = True
    
    def printBreaks(self):
        
        def checkBreakPosition(position):
            if position == 0:
                return self.ref_seq_length
            elif position == self.ref_seq_length + 1:
                return 1
            return position
        
        break_list = []
        for breakpoint, count in self.break_count.items():
            if breakpoint[0]+1 != breakpoint[1]:
                start = checkBreakPosition(breakpoint[0]+1)
                end = checkBreakPosition(breakpoint[1]-1)
                break_list.append([start, end, count])
        
        with open(self.output_header + ".breakpoints.txt", "w") as OUTPUT:
            OUTPUT.write('Start\tEnd\tCount\n')
            for breakpoint in sorted(break_list, key=lambda k: (int(k[0]), int(k[1]), -int(k[2]))):
                OUTPUT.write(str(breakpoint[0]) + '\t' + str(breakpoint[1]) + "\t" + str(breakpoint[2]) + "\n")
    
    def execute(self):
        
        # Make FASTA files from DNA-seq
        self.makeFASTA(self.read_1_fn, self.output_header + ".read_1.fasta")
        self.makeFASTA(self.read_2_fn, self.output_header + ".read_2.fasta")
        
        # Make padded reference
        padded_fn = self.output_header + ".padded_reference.fasta"
        self.makePaddedFASTA(self.ref_fn, padded_fn)
        
        # Read in reference
        ref_seq = self.readReference(self.ref_fn).upper()
        self.ref_seq_length = len(ref_seq)
        
        # Perform alignments
        with open(self.output_header + ".read_1.blat.out", "w") as blat_out:
            call([self.blat_path, padded_fn, self.output_header + ".read_1.fasta", self.output_header + ".read_1.psl"], stdout=blat_out)
        with open(self.output_header + ".read_2.blat.out", "w") as blat_out:
            call([self.blat_path, padded_fn, self.output_header + ".read_2.fasta", self.output_header + ".read_2.psl"], stdout=blat_out)
        
        # Read alignments
        self.readAlignments(self.output_header + ".read_1.psl", self.output_header + ".read_2.psl")

        self.findBreaks()
        self.compareBreaksAcrossReads(ref_seq)
        self.compileBreaks()
        self.compareAcrossAllBreaks(ref_seq)
        self.printBreaks()
        self.cleanFASTA()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('read_1_file_name', type=str, help='Read 1 file')
    parser.add_argument('read_2_file_name', type=str, help='Read 2 file')
    parser.add_argument('reference_sequence', type=str, help='Reference sequence in FASTA format')
    parser.add_argument('output_prefix', type=str, help='Prefix for output file name')
    parser.add_argument('--length', type=int, help='Minimum required alignment length', default=25)
    args = parser.parse_args()

    ROTLA(**vars(args))
