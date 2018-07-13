#!/usr/local/bin/python3
#
#   dna2aa.py
#
#   The dictionary object in Python looks a lot like an array except instead 
#   of accessing a value in the dictionary by its integer index in the 
#   dictionary, values can be stored and accessed using an (almost) arbitrary 
#   index, such as a string or array or almost any "object".  

#   The implementation uses hash tables.  This demo uses strings as the access 
#   "key".  The "values" in the dictionary can also be nearly anything.

import os    # module that handles OS features
import sys   # sys module needed to access command line information

#--------- Start of a Dictionary for codon to amino acid conversion
#----------- global variables -----------------------
usageMsg = '''

Usage: fastaRead fastafile
          Extract each sequence from a fastafile into a header string and
          a sequence string.
              <do something to the sequence>
          This version computes its length and prints the seqId and length
          to standard output.

'''

def usage():
    if len( sys.argv ) < 2 or sys.argv[ 1 ] == "-h":  # command
        print( usageMsg )
        exit( 0 )

# subroutine that translates DNA sequence to amino acid
def dna2aa( seqs ):
    x = 0
    while x < len(seqs):
        seq1 = seqs[ x ]
        seqLen = len( seq1 )
        # assigning value if multiple of 3 and if not accordingly
        # to capture separate nucleotides at the end.
        if (seqLen % 3 == 0) :
            nCodons = int( seqLen / 3 ) # get codon count
        else :
            nCodons = int(seqLen / 3) + 1# get codon count

        aaSeq = ''
        for codonIndex in range( nCodons ) :
            dnaIndex = codonIndex * 3    # convert codon index to a dna index
            codon = seq1[ dnaIndex : dnaIndex + 3 ]  # get the codon from the dna sequence
            if codon in dna_aa:
                aa = dna_aa[ codon ]
            else :
                # using '-' at the end to replace 1 or 2 separate nucleotides
                if (1 <= len(codon) < 3) :
                    aa = '-'
                # using '?' if mapping is not found for the codon
                else :
                    aa = '?'
            aaSeq += aa
        # writing the derived amino acids without header since not a fasta file now
        print( "\n" + aaSeq + "\n")
        x+=1
    print( "-------------- Done -----------------" )

# Genetic Code dictionary
dna_aa = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
    }

print( "------- DNA to AA ---------" )
#***------ Starts: process file
usage()

seqFile = sys.argv[ 1 ]

#----- open file
try:
    inFile = open( seqFile, 'r' )
except ( OSError, IOError ) as e:
    print( 'Unable to open ', seqFile )

#----- read file and process sequences
# first line better be a sequence header
global header
header = inFile.readline()

if not header.startswith( '>' ):
    print( "*** ERROR: Expecting sequence header, got:\n{0}".format( header ),
           end="", file=sys.stderr )
    print( "****** is this a fasta file??? ", file=sys.stderr )
    sys.exit()

#----- read and process the sequences until done

seqs = []
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
    seqs.append(seq)
    #------------------------------------------
    header = line    # last line read is either next header or null
dna2aa(seqs)
