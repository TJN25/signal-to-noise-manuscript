#!/usr/bin/python

'''
file paths are hard coded
'''


import sys
from Bio import SeqIO
import getopt
import os
from BCBio import GFF
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import random
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
        opts, args = getopt.getopt(sys.argv[1:], "a:sqh", ["accession", "shuffle", "quiet", "help"])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    accession = ""
    shuffled = False
    for o, a in opts:
            if o in ("-h", "--help"):
                usage()
                sys.exit()
            elif o in ("-a", "--accession"):
                accession = a
            elif o in ("-s", "--shuffle"):
                shuffled = True
            else:
                assert False, "unhandled option"
    if accession == "":
        print "-a <accession> missing. For more help use -h"
        sys.exit(2)
    return(accession, shuffled)



def main():

    accession, shuffled = rungetopts()
    print "Reading files"
    try:
        inFile = open("/Users/thomasnicholson/phd/RNASeq/new_calls/%s_new_calls.txt" % accession, 'r')

        fileLength = file_len("/Users/thomasnicholson/phd/RNASeq/new_calls/%s_new_calls.txt" % accession)
    except IOError:
        print "/Users/thomasnicholson/phd/RNASeq/new_calls/%s_new_calls.txt not found" % accession
        sys.exit(2)
    try:
        fastaFile = list(SeqIO.parse("/Users/thomasnicholson/phd/RNASeq/sequences/%s.fna" % accession, "fasta"))
    except IOError:
        print "/Users/thomasnicholson/phd/RNASeq/sequences/%s.fna not found" % accession
        sys.exit(2)

    print "Combining contigs"
    my_seq = srna.concatenateSequence(fastaFile)



    print "Getting intergenic sequence"
    random_seq = srna.intergenicSequence(accession, my_seq, shuffled)


    print "Getting intergenic positions"
    positions = srna.intergenicPositions(accession)

    print "Selecting random sRNAs"
    srna.selectRandomLocation(inFile, positions,fileLength, random_seq, accession)





if __name__ == "__main__":
    main()




