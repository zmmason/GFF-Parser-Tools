# !/usr/bin/env python3
# Last Update: June 17, 2021
# Name: Zachary Mason (zmmason@ucsc.edu, zachmason97.zm@gmail.com)

# Copyright (C) 2021  Zachary M Mason

# This program is free software: you can redistribute it and/or modify it under the terms of
# the GNU General Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version. This program is distributed
# in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
# License for more details. You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>

import argparse
from argparse import RawTextHelpFormatter
import timeit
import sys
# python3 GFF_Parser.py -id Gene__ID[6718].txt -gff Mcap_reference_SAB.gff3 -fa Mcap_reference.fa


class CommandLine:
    """ Handle the command line, usage and help requests """

    def __init__(self, in_opts=None):
        """ CommandLine constructor: Implements a parser to interpret the command line argv string using argparse. """
        self.parser = argparse.ArgumentParser(
            description='GFF-Parser\n'
                        'This program is designed to take an input txt file containing significant\n '
                        'feature labels and search a reference GFF file for their corresponding parent sequence\n '
                        'ID, frame, strand, and region index. The target features are then extracted from the\n '
                        'reference FASTA file using the parent sequence titles and indexes associated with each\n '
                        'target feature.\n '
                        'NOTE: Only handles GFF entries with 1 transcript (t1). Additional transcripts will be ignored\n'
                        'unless the target feature is a gene. Start and stop codons will not always begin and end\n'
                        'a sequence. Sequences based on GFF locations given.\n',
            epilog='Copyright (C) 2021  Zachary M Mason\n'
                   'This program comes with ABSOLUTELY NO WARRANTY. This is free software, and you are welcome\n '
                   'to redistribute it under certain conditions. For more information on warranty and conditions,\n '
                   'see the programs associated licence, or go to <https://www.gnu.org/licenses/>\n',
            formatter_class=RawTextHelpFormatter,
            add_help=True,
            prefix_chars='-',
            usage='python3 %(prog)s  [-options]')
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 3.2')
        self.parser.add_argument('-id', action='store', help='input feature/GeneID text (.txt) file')
        self.parser.add_argument('-gff', action='store', help='input reference GFF (.gff/.gff3) file')
        self.parser.add_argument('-fa', action='store', help='input reference FASTA (.fa/.fas/.txt) file')
        self.parser.add_argument('-i', action='store', help='feature of interest',
                                 choices=['mrna', 'mRNA', 'exon', 'intron',
                                          'cds', 'CDS', 'trna', 'tRNA', 'transcript',
                                          'rrna', 'rRNA', 'match', 'gene'],
                                 default='gene')
        self.parser.add_argument('-out', action='store', help='output nucleotide FASTA file (.fa)', default='NT_FASTA.fa')

        if in_opts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(in_opts)


def file_check(feature_List, fnGFF, fnFA, target_feature, correctFile=True):
    """ Checking if appropriate files are given. """
    validFA = ['.fa', '.fas', '.fasta', '.txt']
    validGFF = ['.gff', '.gff3']
    validIDs = ['.txt']

    while correctFile:
        if not feature_List.lower().endswith(tuple(validIDs)):  # file compatibility checks
            print("Incompatible CDS/SNPs/gene file format...")
            correctFile = False
            break
        if not fnGFF.lower().endswith(tuple(validGFF)):
            print("Incompatible GFF file format...")
            correctFile = False
            break
        if not fnFA.lower().endswith(tuple(validFA)):
            print("Incompatible FA file format...")
            correctFile = False
        else:
            feature_list = []
            if target_feature != 'gene':  # prepping ID for specific feature search
                suffix = '.t1'
            else:
                suffix = ''
            for feature in open(feature_List, 'r'):
                if not feature.startswith('adi2'):
                    feature_list.append(feature.strip()+suffix)
                elif feature.startswith('adi2'):
                    feature_list.append(feature.strip())
                    
            if fnGFF.lower().endswith('3'):
                gff_type = 'gff_three'
            else:
                gff_type = 'gff_std'

            return feature_list, gff_type
    while not correctFile:
        print('FILE COMPATIBILITY CHECK: FAILED')
        quit()


class FastaReader:
    """ Define objects to read FastA files. """

    def __init__(self, fname=''):
        """contructor: saves attribute fname """
        self.fname = fname

    def do_open(self):
        """ Handle file opens, allowing STDIN."""
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)

    def read_fa(self):
        """ Read an entire FastA record and return the sequence header/sequence"""
        header = ''
        sequence = ''
        with self.do_open() as fileH:
            header = ''
            sequence = ''
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>'):
                line = fileH.readline()
            header = line[1:].rstrip()
            for line in fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()
        yield header, sequence


class Features:
    """
    Searches reference GFF file for each CDS's/intron's/gene's corresponding parent sequence ID and
    the region index of the significant feature.
    Input:  Annotated reference GFF file (that has its own reference FA file containing parent sequences) and the
            target CDS/intron/gene file.
    Output: list of feature specific data corresponding to the target feature/transcript/intron/gene. This includes the
            feature label, parent sequence number (to reference reference FA file), and feature region in parent sequence.
    COST:   designated in significant functions; M = length of GFF file, N = length of feature's/intron's/gene's list file.
    """
    def __init__(self, fnGFF, gff_type, feature_List, target_feature):
        """ Constructor: saves data from input files. """
        self.GFF = fnGFF
        self.gff_type = str(gff_type)
        self.feature_List = feature_List
        self.target_feature = target_feature.lower()

    def gff_parser(self):
        """
        Parse GFF file for feature specific parent sequence ID and index based on significant feature tags.
        Cost = O(MN+N)
        """
        gffContents = {}
        feature_options = ['mrna', 'exon', 'intron', 'cds', 'trna', 'transcript', 'mrna', 'match']
        frame_options = ['.', '-']
        match = {}
        open_GFF = open(self.GFF)

        for line in open_GFF:
            if not line.startswith('#'):
                transcript = line.strip().split('\t')
                feature = transcript[2].lower()
                if feature == self.target_feature:
                    feature_ID = transcript[-1]

                    if self.gff_type == 'gff_three':  # handling GF3
                        if self.target_feature == 'gene':
                            feature_ID = feature_ID.split(';')[1].split('=')[1]
                        elif self.target_feature in feature_options:
                            feature_ID = feature_ID.split('=')[1]

                    elif self.gff_type == 'gff_std':  # handling GFF
                        if self.target_feature in feature_options:
                            feature_ID.replace('"', '')
                            feature_ID = feature_ID.split(';')[0].split('transcript_id')[1]

                    if feature_ID in self.feature_List:
                        parent_ID = transcript[0]
                        match[feature_ID] = None
                        feature_start = int(transcript[3])
                        feature_stop = int(transcript[4])
                        strand = transcript[6]
                        if transcript[7] in frame_options:
                            frame = 0
                        else:
                            frame = int(transcript[7])
                        if parent_ID not in gffContents:
                            gffContents[parent_ID] = [[feature_ID, feature_start, feature_stop, strand, frame]]
                        elif parent_ID in gffContents:
                            gffContents[parent_ID].append([feature_ID, feature_start, feature_stop, strand, frame])

                    else:
                        continue
        open_GFF.close()
        if len(match) < len(self.feature_List):
            for feature in self.feature_List:
                if feature not in match:
                    print('FEATURE NOT FOUND: ' + feature)

        parent_count = len(gffContents)
        match_count = len(match)
        return gffContents, parent_count, match_count


def reverse_comp_strand(fwd_strand):
    """ Creates a reversed complimentary strand of a DNA sequence. """
    complimentStrand = fwd_strand.lower().replace('a', 'T').replace('t', 'A').replace('c', 'G').replace('g', 'C')
    return complimentStrand[::-1]  # Reverses and returns the new complimentary strand of DNA


class GeneSearch:
    """
    Parse feature data against the reference FA file containing parent sequences and
    generate FASTA feature specific file.
    Input:  Reference FASTA file and feature data (feature label, parent
            sequence number (to reference reference FA file), and feature index region in parent sequence).
    Output: Nucleotide FASTA formatted output containing target feature information and its corresponding sequence.
    COST:   designated in significant functions; M = length of FA file, N = length of feature data list file created
            in gff_parser.
    """
    def __init__(self, fnFA, gffContents, target_feature):
        """ Constructor: saves data from input files. """
        self.FA = fnFA
        self.target_feature = target_feature
        self.FA_title = fnFA.split('.')[0]
        self.gffContents = gffContents
        self.genFASTAList = []

    def gen_fasta(self, nt_name):
        """
        Parse the parent sequence reference FA file with the feature data (parent ID and feature indexes)
        obtained from the GFF file.
        Cost = O(MN)
        """
        myFastaReader = FastaReader(self.FA)
        newFAS = open(nt_name, 'w')
        for parent_ID, seq in myFastaReader.read_fa():
            if parent_ID in self.gffContents:
                feature_Points = self.gffContents.get(parent_ID)
                for feature in feature_Points:
                    feature_ID = feature[0]
                    feature_start = feature[1]
                    feature_stop = feature[2]
                    strand = feature[3]
                    frame = feature[4]
                    feature_title = '{}|{}|NT|{}|{}:{}'.format(feature_ID, parent_ID, self.target_feature, feature_start,
                                                               feature_stop)
                    if strand == '+':
                        newFAS.write(">" + feature_title + "\n" + seq[feature_start+frame-1:feature_stop+frame].upper() + '\n')
                    elif strand == '-':
                        newFAS.write(">" + feature_title + "\n" + reverse_comp_strand(seq[feature_start-frame-1:feature_stop-frame]) + '\n')

        newFAS.close()


def main(command=None):
    """
    Obtains the target feature/intron/gene list and the reference GFF/FA file from commandline and parses the
    information to create both a nucleotide-based FASTA formatted output file
    containing target features. The files are designed to be used in database searches such as BLAST or HMMER.
    Output: Nucleotide fasta containing all features within each significant feature <NT_FASTA.fa>
            *Nucleotide fasta seq titles: feature-title|parent-seq-name|NT|feature|start:stop,
    COST: O(MN+N) + O(MN)
    """

    if command is None:
        my_command = CommandLine()  # read options from the command line
    else:
        my_command = CommandLine(command)  # interpret the list passed from the caller of main

    start = timeit.default_timer()

    print('CHECKING FILE COMPATIBILITY...')
    feature_list, gff_type = file_check(my_command.args.id, my_command.args.gff, my_command.args.fa, my_command.args.i)  # file check
    print('---FILE COMPATIBILITY CHECK: PASSED---')

    print('SEARCHING {} FEATURES FOR CORRESPONDING SEQUENCE LOCATION...' .format(len(feature_list)))
    my_feature = Features(my_command.args.gff, gff_type, feature_list, my_command.args.i)
    gffContents, parent_count, match_count = my_feature.gff_parser()

    if len(feature_list) - match_count >= 0:
        unmatched = len(feature_list) - match_count
    else:
        unmatched = 0

    print("---NUMBER OF {}s MATCHED: {}---\n---NUMBER OF PARENT SEQUENCES: {}---\n---NUMBER OF UNMATCHED FEATURES {}---"
          .format(my_command.args.i.upper(), match_count, parent_count, unmatched))

    if match_count == 0:
        stop = timeit.default_timer()
        print('NO DATA...')
        print('####################################################\nRun Time: {}'
              '\n####################################################'.format(stop - start))
        quit()

    print('GENERATING NUCLEOTIDE FASTA FILE...')
    myGene = GeneSearch(my_command.args.fa, gffContents, my_command.args.i)
    myGene.gen_fasta(my_command.args.out)
    print('---NT FASTA COMPLETE: "' + my_command.args.out + '"---')

    stop = timeit.default_timer()
    print('####################################################\nRun Time: {}'
          '\n####################################################' .format(stop - start))


if __name__ == "__main__":
    main()
