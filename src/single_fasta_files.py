#!/usr/bin/python
import sys
from Bio import SeqIO
import getopt
import os
from BCBio import GFF
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import random
import pandas as pd
import comparativesrna as srna

help = '''
    {script_name} -c com_port [-o output_file] [--loglevel level]

    Reads the temperature data from a radio.  The temperature data is output in csv form.

    examples:
        Read table from radio attached to com4 and write the table to the file
        output.csv.

            {script_name} -c com4 -o output.csv

        Read table from radio attached to com3 and write the table to stdout. 
        You can use IO redirection to send the contents where every you want.

            # just print to the terminal 
            {script_name} -c com3

            # redirect to another file
            {script_name} -c com3 > somefile.csv

            # filter out temperatures that are -100
            {script_name} -c com3 | grep -v '^-100' 


    -c com_port
    --com_port comport
        Name of the COM port attached to the radio

    -o output_file
    --output output_file
        If specified write the table data to the given file.  If not specified
        the data will be written to stdout.

    --loglevel critical | error | warning | info | debug | notset
        Control the verbosity of the script by setting the log level.  Critical
        is the least verbose and notset is the most verbose.

        The default loglevel is {default_loglevel}.

        These values correspond directly to the python logging module levels. 
        (i.e. https://docs.python.org/3/howto/logging.html#logging-levels)


    -h 
    --help 
        print this message

'''

def usage():
    print help

def rungetopts():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:f:qh", ["infile", "folder", "quiet", "help"])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    file = ""
    for o, a in opts:
            if o in ("-h", "--help"):
                usage()
                sys.exit()
            elif o in ("-i", "--infile"):
                file = a
            elif o in ("-f", "--folder"):
                folder = a
            else:
                assert False, "unhandled option"
    if file == "":
        print "-f <file> missing. For more help use -h"
        sys.exit(2)
    return(file, folder)

# def single_fasta(fastaFile, folder):
#     for seq in fastaFile:
#         id = seq.id
#         outname = id.split("[")
#         outname = outname[0]
#         my_seq = seq.seq
#         outFile = open("/Users/thomasnicholson/phd/RNASeq/srna_seqs/version_1/%s/%s.fna" % (folder, outname), "w")
#         outFile.write(">%s\n%s\n" % (id, my_seq))



def main():

    filename, folder = rungetopts()
    print filename, folder
    print "Reading files"
    try:
        fastaFile = list(SeqIO.parse("/Users/thomasnicholson/phd/RNASeq/srna_seqs/version_1/%s" % filename, "fasta"))
    except IOError:
        print "/Users/thomasnicholson/phd/RNASeq/srna_seqs/version_1/%s not found" % filename
        sys.exit(2)

    srna.single_fasta(fastaFile, folder)


if __name__ == "__main__":
    main()
