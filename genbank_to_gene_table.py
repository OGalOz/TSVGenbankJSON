#python3
import logging
import os
import sys
import subprocess
from bioparsing import parseFASTA
#from BCBio import GFF
from Bio import SeqIO



'''
The gene table we'll build will look like:

            locusId (str):sysName (str):type (int):scaffoldId (str):begin (int):end (int):
                strand (str +/-):name (str):desc (str):GC (float [0,1]):nTA (int)
'''


def convert_genbank_to_gene_table(genbank_filepath, gt_filepath, gffPrlScriptPath):
    """
    All inputs are str paths
    """



    genbank_filepath = os.path.abspath(genbank_filepath)
    gt_filepath = os.path.abspath(gt_filepath)

    in_handle = open(genbank_filepath, "r")
    out_handle = open(gt_filepath +".gff", "w")
    
    GFF.write(SeqIO.parse(in_handle, "genbank"), out_handle)
    
    in_handle.close()
    out_handle.close()

    gffToGeneTable_command = "perl " + gffPrlScriptPath + " < " + gt_filepath + ".gff > " + gt_filepath
    subprocess.call(gffToGeneTable_command, shell=True)
    logging.info("Wrote Gene Table file to " + gt_filepath)
    return 0




def genbank_and_genome_fna_to_gene_table(gbk_fp, gnm_fp, op_fp):
    """
    Args:
        gbk_fp: (str) Path to genbank file.
        gnm_fp: (str) Path to genome fna file
        op_fp: (str) Path to write genes table to

    We use GenBank Records and Features

    """

    
    id2seq = parseFASTA(gnm_fp)

    out_FH = open(op_fp, 'w')
    # This is the output file start line:
    out_FH.write("locusId\tsysName\ttype\t" + \
                    "scaffoldId\tbegin\tend\tstrand\t" + \
                    "name\tdesc\tGC\tnTA\n")
    
    gb_record_generator = SeqIO.parse(open(gbk_fp, "r"), "genbank")

    # TYPE - Note that we don't like misc_feature or gene
    # May need to skip anything besides values under 10
    types_dict = {"CDS" : 1, "rRNA" : 2, "tRNA" : 5, 
                   "RNA" : 6, "transcript" : 6,
                   "pseudogene": 7, "misc_feature": 20, 
                   "gene": 21}

    for gb_record in gb_record_generator:

        locus_tag = gb_record.name
        #Genome sequence:
        g_seq = str(gb_record.seq)

        scaffold = findSeqInId2Seq(g_seq, id2seq)

        g_len = len(g_seq)

        #Genome features (list of features):
        g_features = gb_record.features
        #DEBUG
        #print(g_features[0])
        g_feat_len = len(g_features)

        """
        scaffoldId_exists= False
        if "scaffold_name" in config_dict:
            scaffold_id = config_dict["scaffold_name"]
            scaffoldId_exists = True
        """

        try:
            for i in range(g_feat_len):


                current_feat = g_features[i]

                if current_feat.type == "source":
                    print("Feature is of type 'source':")
                    print(current_feat)
                    continue

                
                #scaffold already set above
                begin = "null"
                end = "null"
                strand = "null"
                desc = "null"
                typ = "null"
                #locus_tag = "null"
                sysName = "null"
                name = "null"
                GC = "null"
                nTA = "null"


                """
                # Scaffold Id
                if scaffoldId_exists:
                    if scaffold_id in g_features[i].qualifiers:
                        scaffold = g_features[i].qualifiers[scaffold_id]
                    else:
                        logging.debug("Could not find scaffold id "
                                "{} in qualifiers:".format(scaffold_id))
                        logging.debug(g_features[i].qualifiers)
                        scaffold = "1"
                else:
                    scaffold = "1"

                """

                # Begin
                begin = str(current_feat.location.start + 1)
                # End
                end = str(current_feat.location.end + 1)

                # Strand
                if current_feat.strand == 1:
                    strand = "+"
                elif current_feat.strand == -1:
                    strand = "-"
                else:
                    logging.critical("Could not recognize strand type.")
                    raise Exception("Parsing strand failed.")

                # Desc (Description)
                if "product" in current_feat.qualifiers.keys():
                    desc = str(current_feat.qualifiers['product'][0])
                else:
                    desc = "Unknown function" 
                    logging.critical("Could not find description in current_feat: ")
                    logging.critical(current_feat)
                    continue

                typ_str = current_feat.type.strip()
                if typ_str in types_dict:
                    typ = str(types_dict[typ_str])
                else:
                    logging.info("Could not recognize type from feature: " \
                            + typ_str)
                    typ = "0"
                if typ == "1":
                    out_FH.write("\t".join([locus_tag, sysName, typ, scaffold,
                            begin, end, strand, name, desc, GC, nTA]) + "\n")
                else:
                    logging.info(f"Did not write annotation for gene with type {typ}")
        except:
            logging.critical("Could not parse all features in genbank file.")
            out_FH.close()
            raise Exception("Parsing genbank file into gene table failed")
    
    out_FH.close()

    return None




def findSeqInId2Seq(seq_str, id2seq):

    seq_id = None

    for s in id2seq.keys():
        if seq_str == id2seq[s]:
            seq_id = s
            break

    if seq_id == None:
        raise Exception("Couldn't find Sequence ID for sequence:\n" + \
                        f" '{seq_str}'\n\nin id2seq")

    return seq_id


"""
    

The genes table must include the fields
   scaffoldId, begin, end, strand, desc, type
        'type' with type=1 for protein-coding genes,
    other fields:
        sysName, name, GC, nTA
"""
def OLD_convert_genbank_to_gene_table(genbank_filepath, output_filepath, gff_fasta_gbk=True,
                                      config_dict={}):

    """
    Args:
        genbank_filepath: (str) Path to genbank file.
        output_filepath: (str) Path to write genome table to
        gff_fasta_gbk: bool True/False
        config_dict:
            [gff_fasta_gbk]: True/False
            [scaffold_name]: (str) The scaffold id in the genbank file
            [keep_types]: (list) If you only want specific types
    """
    
    out_FH = open(output_filepath, 'w')

    # This is the output file start line:
    output_header = "locusId\tsysName\ttype\t" + \
                    "scaffoldId\tbegin\tend\tstrand\t" + \
                    "name\tdesc\tGC\tnTA\n"

    out_FH.write(output_header)

    # We use BioPython SeqIO to parse genbank file:
    # https://biopython.org/DIST/docs/api/Bio.SeqIO-module.html
    gb_record_generator = SeqIO.parse(open(genbank_filepath, "r"), "genbank")

    gb_record = next(gb_record_generator)


    record_name = gb_record.name
    scaffold = record_name 
    #Genome sequence:
    g_seq = gb_record.seq
    g_len = len(g_seq)

    #Genome features (list of features):
    g_features = gb_record.features
    #DEBUG
    print(g_features[0])
    g_feat_len = len(g_features)

    """
    scaffoldId_exists= False
    if "scaffold_name" in config_dict:
        scaffold_id = config_dict["scaffold_name"]
        scaffoldId_exists = True
    """

    try:
        for i in range(g_feat_len):
            
            #scaffold already set above
            begin = "null"
            end = "null"
            strand = "null"
            desc = "null"
            typ = "null"
            locus_tag = "null"
            sysName = "null"
            name = "null"
            GC = "null"
            nTA = "null"

            current_feat = g_features[i]

            print(current_feat)
            print(current_feat.qualifiers)
            
            """
            # Scaffold Id
            if scaffoldId_exists:
                if scaffold_id in g_features[i].qualifiers:
                    scaffold = g_features[i].qualifiers[scaffold_id]
                else:
                    logging.debug("Could not find scaffold id "
                            "{} in qualifiers:".format(scaffold_id))
                    logging.debug(g_features[i].qualifiers)
                    scaffold = "1"
            else:
                scaffold = "1"

            """

            # Begin
            begin = str(current_feat.location.start)
            # End
            end = str(current_feat.location.end)

            # Strand
            if current_feat.strand == 1:
                strand = "+"
            elif current_feat.strand == -1:
                strand = "-"
            else:
                logging.critical("Could not recognize strand type.")
                raise Exception("Parsing strand failed.")

            # Desc (Description)
            if "product" in current_feat.qualifiers.keys():
                desc = str(current_feat.qualifiers['product'][0])
            else:
                if gff_fasta_gbk and 'locus_tag' in current_feat.qualifiers:
                    desc = " ".join(current_feat.qualifiers['locus_tag'][
                        0].split(' ')[1:]).split('=')[-1]
                else:
                    desc = "Unknown function" 
                    logging.critical("Could not find description in current_feat")


            # Locus ID:
            if "locus_tag" in current_feat.qualifiers.keys():
                if gff_fasta_gbk:
                    locus_tag = str(current_feat.qualifiers['locus_tag'][0].split(' ')[0])
                else:
                    locus_tag = str(current_feat.qualifiers['locus_tag'][0])
            else:
                locus_tag = "Unknown_Locus_tag."
                logging.critical("Could not find locus tag in current_feat")

            # TYPE - Note that we don't like misc_feature or gene
            # May need to skip anything besides values under 10
            types_dict = {"CDS" : 1, "rRNA" : 2, "tRNA" : 5, 
                            "RNA" : 6, "transcript" : 6,
                            "pseudogene": 7, "misc_feature": 20, 
                            "gene": 21}
            typ_str = current_feat.type.strip()
            if typ_str in types_dict:
                typ = str(types_dict[typ_str])
            else:
                logging.info("Could not recognize type from feature: " \
                        + typ_str)
                typ = "0"
            if typ == "1":
                out_FH.write("\t".join([locus_tag, sysName, typ, scaffold,
                        begin, end, strand, name, desc, GC, nTA]) + "\n")
    except:
        logging.critical("Could not parse all features in genbank file.")
        raise Exception("Parsing genbank file into gene table failed")

    out_FH.close()

    logging.info("Wrote Gene Table to " + output_filepath)

    '''
    below is mainly broken
    output_file_string = open(output_filepath,'r').read()
    #We remove duplicate gene lines and remove the last new line symbol
    output_file_string = unduplicate_gene_table(output_file_string)
    '''

    '''
    if "keep_types" in config_dict:
        types_to_keep = config_dict["keep_types"]
        output_file_string = keep_types_gene_table(output_file_string, 
                                                    types_to_keep)
    '''
    '''
    with open(output_filepath, "w") as f:
        f.write(output_file_string)
    '''
    return output_filepath



"""
This function removes duplicates from the gene table
Inputs:
    gene_table_string: (str) A string of the entire gene table file
Outputs:
    gene_table_string: (str) A string of the entire gene table file
Process: 
    Compares the location of features and if they are the same removes
        one of the two.
"""
def unduplicate_gene_table(gene_table_string):


    #first we split the gene_table into lines:
    split_list = gene_table_string.split("\n")
    header_line = split_list[0] + "\n" 
    gt_lines = split_list[1:]
    logging.info("Total number of lines besides headers: " + str(len(gt_lines)))

    #Then for each line we check if it's a duplicate or not.
    #We add the indices of duplicate lines and then remove the lines in reverse order.
    # 'loc' means begin to end in sequence
    splitLine = gt_lines[0].split("\t")
    print(splitLine)
    existing_loc = splitLine[1:3]; existing_typ = splitLine[2]
    previous_index = 0
    # We create a set, duplicate_line_indices
    duplicate_line_indices = set()
    print(len(gt_lines))
    for i in range(1, len(gt_lines) - 1):
        splitLine = gt_lines[i].split("\t")
        print(i)
        current_loc = splitLine[1:3]; crnt_typ = splitLine[-1]
        if (current_loc[0] == existing_loc[0]) or \
                (current_loc[1] == existing_loc[1]):
            if crnt_typ == '1':
                if existing_typ == '1':
                    logging.warning("Two overlapping location Protein "
                            "Features: loc: {}{},{}{} types: {},{}".format(
                                existing_loc[0], existing_loc[1],
                                current_loc[0], current_loc[1],
                                existing_typ, crnt_typ))
                duplicate_line_indices.add(previous_index)
                previous_index = i
            else:
                if existing_typ == '1':
                    duplicate_line_indices.add(i)
                else:
                    duplicate_line_indices.add(previous_index)

        else:
            existing_loc = current_loc
            previous_index = i
    logging.info("Duplicate Lines: " + str(len(duplicate_line_indices)))

    # Sorting indices so they ascend 
    duplicate_line_indices = sorted(list(duplicate_line_indices))
    #removing the indeces in reverse order:
    duplicate_line_indices.reverse()
    for i in range(len(duplicate_line_indices)):
        del gt_lines[duplicate_line_indices[i]]

    logging.info("New total line number (after duplicate line removal): " \
            + str(len(gt_lines)))



    #Converting list into string again:
    gene_table_string = header_line + "\n".join( gt_lines)


    return gene_table_string


"""
Inputs:
    gene_table_string: (str) The gene table string
    types_to_keep: list<str> Each string in list is a type we want.
Outputs:
    gene_table_string: (str) The gene table string.
"""
def keep_types_gene_table(gene_table_string, types_to_keep):

    split_list = gene_table_string.split("\n")
    header_line = split_list[0] + "\n" 
    gt_lines = split_list[1:]

    non_good_type_indices = []

    #For each line, we check if its type is 1. If not, we remove it later.
    for i in range(len(gt_lines)):
        current_type = gt_lines[i].split("\t")[-1]
        if current_type not in types_to_keep:
            non_good_type_indices.append(i)

    #removing the indeces in reverse order:
    non_good_type_indices.reverse()
    for i in range(len(non_good_type_indices)):
        del gt_lines[non_good_type_indices[i]]


    logging.info("New total line number (after type 1): " + str(len(gt_lines)))

    #Converting list into string again:
    gene_table_string = header_line + "\n".join(gt_lines)

    return gene_table_string



def test(args):
    logging.basicConfig(level=logging.DEBUG)
    gb_fp = args[1]
    op_fp = args[2]
    #prlScript = args[3]

    #config_dict = {"keep_types": ["1","5","6"]}
    OLD_convert_genbank_to_gene_table(gb_fp, 
                                    op_fp)

def main():
    """
    args should be genbank_to_gene_table.py genbank, output, gffPrlScript
    """
    args = sys.argv
    if args[-1] != "1":
        print("python3 genbank_to_gene_table.py genbank_fp genome_fna output")
        sys.exit(0)
    else:
        args = sys.argv
        gbk_fp = args[1]
        gnm_fp = args[2]
        op_fp = args[3]
        genbank_and_genome_fna_to_gene_table(gbk_fp, gnm_fp, op_fp)
        print(f"Wrote gene table to {op_fp}")


if __name__ == "__main__":
    main()
