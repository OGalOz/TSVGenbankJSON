#!python3

import os
import sys
import logging
import json
import itertools as itrt
from Bio import SeqIO

'''
Supported parsing formats: (https://biopython.org/wiki/SeqIO#file-formats)
    abi, abi-trim, ace, cif-atom, cif-seqres, clustal,
    embl, fasta, fasta-2line, fastq-sanger or fastq, 
    fastq-solexa, fastq-illumina, gck, genbank or gb,
    ig, imgt, nexus, pdb-seqres, pdb-atom, phd, phylip,
    pir, seqxml, sff, sff-trim, snapgene, stockholm, 
    swiss, tab, qual, uniprot-xml, xdna 

SeqRecord:
    SeqIO always returns SeqRecord objects. SeqRecord objects have the following 
    methods: 'annotations', 'dbxrefs', 'description', 'features', 'format', 
            'id', 'letter_annotations', 'name', 'seq'
    if not more methods.
    Attributes:
        seq: A Seq Object
        id: (str)
        [name]: (str)
        [description]: (str)
        [dbxref]: list<str>
        [features]: list<SeqFeature>
        annotations: dict 
            Not clearly defined
        letter_annotations: dict Per letter/symbol annotation
                annotation_type -> values (?)
                    Each value holds Python sequences (lists, strings, or 
                    tuples) whose length matches that of the sequence.
                    A typical example is a list of integers representing
                    sequence quality scores
SeqFeature:
    Attributes:
        location - FeatureLocation
            (location of the feature on the sequence)
        type: <str> type of the feature
        location_operator: (str) complex
        strand: int [1,0,-1]
            1 -> plus strand
            -1 -> minus strand
            0 -> stranded but unknown
        id: <str>
        ref: <str?>
        ref_db: A differnt database for the reference accession number
        qualifiers: dict
            qualifier_name (str) -> qualifier value (?)
    

FeatureLocation:
     start: + int
     end: + int
     [strand]: int 1 (+), -1 (-), 0 (?)   
     [ref]: <str>
     [ref_db]: <str> (?)

     creating a FeatureLocation returns a loc which is an
        iterator over all positions 

'''


def parse_tsv(inp_tsv_fp, headers=True):
    """
    Returns:
            "total_line_number": int, Number of lines in file
            "num_columns": int, Total Number of columns
            "header_d": header_d, dict: column -> index
            "matrix": matrix the file in list format (list<list>)
    """

    tsv_fh = open(inp_tsv_fp, "r")
    matrix = []
    line_num = 1

    if headers:
        headers = tsv_fh.readline().rstrip().split('\t')
        num_col = len(headers)
        header_d = {headers[i]:i for i in range(num_col)}
    else:
        first_line = tsv_fh.readline().split('\t')
        num_col = len(first_line)
        header_d = None
        first_line[-1] = first_line[-1].rstrip()
        matrix.append(first_line)
        line_num = 2

    c_line = tsv_fh.readline()
    while c_line != "":
        c_list = c_line.split('\t')
        c_list[-1] = c_list[-1].rstrip()
        num_val = len(c_list)
        if  num_val != num_col:
            raise Exception(f"Incorrect number of values on line {line_num}.\n" + \
                    f"Expected {num_col} but got {num_val}. File: {inp_tsv_fp}")
        matrix.append(c_list)
        line_num += 1
        c_line = tsv_fh.readline()

    tsv_fh.close()



    return {
            "total_line_number": line_num - 1,
            "num_columns": num_col,
            "header_d": header_d,
            "matrix": matrix
    }


def parseFASTA(fasta_fp):
    """
    Args:
        fasta_fp: filepath to FASTA file
    Returns:
        id2seq: (dict)
            Goes from sequence ID/ name (str) -> sequence (str)
    """
    id2seq = {}
    seq_generator = SeqIO.parse(fasta_fp, "fasta")
    for sequence in seq_generator:
        id2seq[sequence.id] = str(sequence.seq)

    return id2seq



def parseGenBankExample(gbk_fp):
    """
    We want to take the genbank file and move it to a format we can easily work with.
    Genbank JSON (?)

    Multiple record names?
    Multiple feature lists?
    What is included in each record, what in each feature?
    """
    gbk_record_iterator = SeqIO.parse(gbk_fp, "genbank")

    copy_1, copy_2 = itrt.tee(gbk_record_iterator)
  
    first_SeqRecord = next(copy_1)
    print(dir(first_SeqRecord))
    ftrs = first_SeqRecord.features
    print(ftrs)
    print("first type:")
    ftr1 = ftrs[0]
    print('qualifiers')
    print(ftr1.qualifiers)
    

    '''
    for i in range(2):
        crnt_rec = next(gbk_record_iterator)
        print(crnt_rec.features)
    '''
    '''
    i = 0
    for x in copy_1:
        print(i)
        i += 1
    '''
    return gbk_record_iterator





def test(crnt_fp):
    my_d = parseFASTA(crnt_fp)
    print(my_d.keys())
    print(my_d[list(my_d.keys())[0]])
    return None


def main():
    
    args = sys.argv
    if args[-1] not in ["1"]:
        print("Incorrect args. Use the following:\n")
        help_str = "python3 FileName.py genes.gff genome.fna 1"
        print(help_str)
        sys.exit(0)
    else:
        if args[-1] == "1":
            gff_fp = args[1]
            fasta_fp = args[2]

    return None

if __name__ == "__main__":
    main()




