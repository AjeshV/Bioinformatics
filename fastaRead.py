#!/usr/local/bin/python3
#
#  fastaRead.py fastafile
#               Read a FastA file representing genome sequence data,
#                   read 1 sequence header and its data at a time
#                   do <something> with the sequence
#                   output <something>
#               This version computes the length of the sequence and outputs
#               the sequence header line with the length info appended.
# Purpose:
# provides summary, statistics information from reading a fastafile
# contacts processSequence after finding header to get information for particular sequence
# uses map (key, value pairs) data structure to get maximum number and percentage of non-gcat chars

import os    # module that handles OS features
import sys   # sys module needed to access command line information
import statistics
import operator
from collections import Counter

#----------- global variables -----------------------
usageMsg = '''
Usage: fastaRead fastafile
          Extract each sequence from a fastafile into a header string and
          a sequence string.
              <do something to the sequence>
          This version computes its length and prints the seqId and length
          to standard output.
'''
#-------------------------- usage() ---------------------------------------
#----- usage test
#      If we have no arguments or the first argument is -h
#      print some "usage" information and quit.
#
# Note that the sys.argv array includes the "programName" in positions 0
#      and command line arguments begin at position 1.

def usage():
    if len( sys.argv ) < 2 or sys.argv[ 1 ] == "-h":  # command
        print( usageMsg )
        exit( 0 )
#---------------------------------------------------------------------
#----------------------- processSequence -----------------------------
# Given a header and a sequence as a single character string with no line feeds,
#   do whatever processing is desired for that sequence.
# The starting code just computes the sequence size and prints the seq id and
#   size to standard output.
#


def processSequence( header, sequence ):
    # ------------------------------------------
    #  Supplement or replace the lines below with the 
    #  sequence processing and output you want to do.
    #
    #basesCount = len( sequence )
    sequences.append(sequence)
    basesCount.append(len(sequences[i]))
    x = Counter(sequences[i])
    #all = Counter(''.join(sequences))
    gc1 = x['G'] + x['C']
    gc2 = x['A'] + x['T']
    tmp =  (len(sequences[i]) - (gc1+gc2))
    gcContent = gc1 / (gc1+gc2)
    nonGcat.append((basesCount[i] - (gc1+gc2)))
    nonGcat2 = nonGcat[i]/basesCount[i]

    # extract the sequence id from the header
    header = header + " " # add a space at end; guarantees that index finds one
    spacePos = str.find( header, ' ' ) # find first space in header
    seqId = header[ 1 : spacePos ]     # [ start : end ] is python "slice" operator
                                       #  end is 1st position NOT in the slice
                                       # slice works for both strings and arrays.
    #
    # print automatically adds newline appropriate for current OS
    print( 'seqId=' + seqId  + ' length=' + str( basesCount[i] ))
    print( 'summary of characters in sequence = '+ str(x))
    print( 'gc-Content=' + str(gcContent))
    print( 'non-Gcat=' + str(nonGcat2))
    gcContentAll.append(gcContent)
    nonGcatAll.append(nonGcat2)
    if tmp >= 1:
        newdict[tmp] = str(basesCount[i])

#-----------------------------------------------------------------------
#--------------------------- main --------------------------------------
usage()

seqFile = sys.argv[ 1 ]

#----- open file
try:
    inFile = open( seqFile, 'r' )
except ( OSError, IOError ) as e:
    print( 'Unable to open ', seqFile )

#----- read file and process sequences
# first line better be a sequence header
header = inFile.readline()

if not header.startswith( '>' ):
    print( "*** ERROR: Expecting sequence header, got:\n{0}".format( header ),
           end="", file=sys.stderr )
    print( "****** is this a fasta file??? ", file=sys.stderr )
    sys.exit()

#----- read and process the sequences until done
count = 0
basesCount = []
sequences = []
gcContentAll = []
nonGcat = []
nonGcatAll = []
seqwnonGcat = []
newdict = {}
global i
global j
global max_value
global max_key
i = 0
#j = 0

while header != '':
    header = header.rstrip( os.linesep )   # delete line separator for any OS
    #-------------------------------------------------
    # Non-regex code for extracting the sequence Id from the header

    #
    seq = ''
    line = inFile.readline()  # returns empty string on eof

    while  line != '' and not line.startswith( '>' ):
        line = str.rstrip( line ) # delete trailing white space including lf
        seq += line               # append line sequence data to seq
        line = inFile.readline()  # read next line
    count += 1
    processSequence( header, seq )
    i += 1

    #------------------------------------------
    header = line    # last line read is either next header or null
print( 'Statistics:' )
print( '1. Number of Sequences = ' + str(count))
avgSequenceLength = statistics.mean(basesCount)
print( '2. Average Sequence Length = ' +  str(avgSequenceLength))
stddevSequenceLength = statistics.pstdev(basesCount)
print( '3. Standard Deviation = ' +  str(stddevSequenceLength))
gcContentAllTotal = sum(map(float,gcContentAll))
print( '4. GC-content of all sequences in file = ' + str(gcContentAllTotal))
avgSequencegcContent =  statistics.mean(gcContentAll)
print( '5. Average sequence GC-content = ' + str(avgSequencegcContent))
stddevgcContent = statistics.pstdev(gcContentAll)
print( '6. Standard Deviation of GC-content = ' +  str(stddevgcContent))
#nonGcatAll2 = sum(map(float,nonGcatAll))
nonGcatAll2 = sum(nonGcat)
print( '7. Total number of non-Gcat Characters = ' + str(nonGcatAll2))
seqwnonGcat2 =  int((len(newdict)/count)*100)
print( '8. Percentage of sequences with 1 or more non-Gcat Characters = ' + (str(seqwnonGcat2)))
max_key = max(newdict)
for key, value in newdict.items():
    if key == max_key:
        max_value = int(value)
maxPercentage = int((max_key/max_value)*100)
print( '9. Max number of non-Gcat character in a sequence = ' + str(max_key))
print( '10. Max percentage of non-Gcat character in a sequence = ' + str(maxPercentage))
print( "-------------------------------------------------------------------" )
#-------------------- end main -------------------------------------
