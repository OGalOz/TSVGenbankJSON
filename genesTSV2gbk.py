#!python3

import os
import sys
import logging
import json
import itertools as itrt
import datetime
from Bio import SeqFeature
from Bio.GenBank import Record
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from bioparsing import parseFASTA, parse_tsv
from print_genbank_features import *
from Bio import SeqIO


"""
Info about BioPython objects:

    SeqFeature.FeatureLocation
        start
        end
        strand
        ref
        ref_db

    SeqFeature.SeqFeature
    BioPython SeqFeature:
        location - the location of the feature on the sequence (FeatureLocation)
        type - the specified type of the feature (ie. CDS, exon, repeat…)
        location_operator - a string specifying how this SeqFeature may be related to others. For example, in the example GenBank feature shown below, the location_operator would be “join”. This is a proxy for feature.location.operator and only applies to compound locations.
        strand - A value specifying on which strand (of a DNA sequence, for instance) the feature deals with. 1 indicates the plus strand, -1 indicates the minus strand, 0 indicates stranded but unknown (? in GFF3), while the default of None indicates that strand doesn’t apply (dot in GFF3, e.g. features on proteins). Note this is a shortcut for accessing the strand property of the feature’s location.
        id - A string identifier for the feature.
        ref - A reference to another sequence. This could be an accession number for some different sequence. Note this is a shortcut for the reference property of the feature’s location.
        ref_db - A different database for the reference accession number. Note this is a shortcut for the reference property of the location
        qualifiers - A dictionary of qualifiers on the feature. These are analogous to the qualifiers from a GenBank feature table. The keys of the dictionary are qualifier names, the values are the qualifier values. As of Biopython 1.69 this is an ordered dictionary.



    Record.Record
    BioPython Record:
    https://biopython.org/docs/1.75/api/Bio.GenBank.Record.html
    Genbank Record Attributes:
        locus: (str)
        accession: list<>
        base_counts: (str)
        comment: (str)
        contig: (str)
        data_file_division: (str)
        nid: (str) Nucleotide identifier number
        db_source: (str) Information about the database record came from
        gi: (str) The NCBI gi identifier for the record
        keywords: list<str> A list of keywords related to the record
        source: (str) The source of material where the sequence came from
        organism: (str) The genus and species of the organism (i.e. Homo sapiens)
        taxonomy: list<str> A listing of the taxonomic classification of the organism,
                        starting general and getting more specific

"""

def ConvertTSVToGBK(tsv_fp, fasta_fp, op_gbk_fp, type2ix=None,
                    accession_list=[], locus_base="", json_fp=None):

    TSV_Info =  parse_tsv(tsv_fp, headers=True)
    print(TSV_Info)
    print(TSV_Info["header_d"])
    check_headers(TSV_Info["header_d"])
    id2seq = parseFASTA(fasta_fp)
    id2seqfeat = {x:[] for x in list(id2seq)}
    print(id2seqfeat)

    """
    full_seq = ""
    for k in list(id2seq):
        full_seq += id2seq[k]
    """
    

    for i in range(len(TSV_Info["matrix"])):
        #for i in range(10):
        row = TSV_Info["matrix"][i]
        seq_feat, scaffoldId = CreateSeqFeatureFromTSVRow(id2seqfeat, 
                                            row, 
                                            TSV_Info["header_d"], 
                                            id2seq, 
                                            str(i))
        if scaffoldId in id2seqfeat:
            id2seqfeat[scaffoldId].append(seq_feat)
        else:
            raise Exception(f"scaffoldId {scaffoldId} not in id2seqfeat")

    # We add the required 'source' feature to all these feature lists
    ids_l = list(id2seqfeat)
    for seq_id in ids_l:
        source_feature = create_source_feature(seq_id, id2seq[seq_id])
        id2seqfeat[seq_id].insert(0, source_feature)

    f_str_list = []
    dt = datetime.datetime.today()
    date_str = f"{dt.day}-{dt.month}-{dt.year}" 

    records = []
    features_file_str = ""
    for i in range(len(ids_l)):
        x = ids_l[i]
        gbk_record = createGenBankRecord(locus_base + "_" + str(i+1), 
                                        id2seqfeat[x],
                                        id2seq[x],
                                        x,
                                        date_str,
                                        [x] 
                                        )
        # func from print_genbank_features
        #features_file_str += get_features_str(gbk_record)
        records.append(gbk_record)
        f_str_list.append(get_gbk_str(gbk_record))

    with open(op_gbk_fp, "w") as g:
        g.write("\n".join(f_str_list))

    if json_fp is not None:
        gbk_json = []
        for rc in records:
            gbk_json.append(createGenBankRecordDict(rc))
        with open(json_fp, 'w') as g:
            g.write(json.dumps(gbk_json, indent=2))


    # my_fh2 = open("test_gbk_feats.gbk", "w")
    # my_fh2.write(features_file_str)
    # my_fh2.close()
    
    return id2seqfeat




def createGenBankRecord(locus_name, 
                        features_list, 
                        seq_str,
                        nid,
                        date_str,
                        accession_list,
                        contig="",
                        taxonomy = [],
                        molecule_type="DNA",
                        residue_type="DNA",
                        data_file_division="BCT"):
    """
    Args:
        features_list: list<seqfeat>
                        seqfeat: SeqFeature.SeqFeature from BioPython (desc top)

    data_file_division: "BCT" - bacteria, "VRL": viral sequences, "PHG": bacteriophage
                        "SYN" - synthetic sequences, "UNA": unannotated sequences


    GenBank record description at top of file.
    """
    
    # initializing genbank record
    gbk_record = Record.Record()

    # Adding attributes
    gbk_record.accession = accession_list 
    gbk_record.base_counts = ""
    gbk_record.comment = ""
    gbk_record.contig = contig
    gbk_record.data_file_division = data_file_division
    gbk_record.date = date_str
    gbk_record.db_source = ""
    gbk_record.dblinks = []
    gbk_record.definition = ""
    gbk_record.features = features_list 
    gbk_record.gi = ""
    gbk_record.keywords = []
    gbk_record.locus = locus_name
    gbk_record.molecule_type = molecule_type 
    # is NID deprecated?
    gbk_record.nid = nid
    gbk_record.organism = ""
    gbk_record.origin = ""
    gbk_record.pid = ""
    gbk_record.primary = []
    gbk_record.projects = []
    gbk_record.references = []
    gbk_record.residue_type = residue_type
    gbk_record.segment = ""
    gbk_record.sequence = seq_str
    gbk_record.size = str(len(seq_str))
    gbk_record.source = ""
    gbk_record.taxonomy = []
    gbk_record.topology = ""
    gbk_record.version = ""
    gbk_record.wgs = ""
    gbk_record.wgs_scafld = []

    return gbk_record

def createGenBankRecordDict(gbk_record):
    """
    Part of the attempt to createa genbank JSON
        object
    """
    record_d = {
        "accession" : gbk_record.accession,
        "base_counts" :  gbk_record.base_counts,
        "comment" :  gbk_record.comment,
        "contig" :  gbk_record.contig,
        "data_file_division" : gbk_record.data_file_division,
        "date" :  gbk_record.date,
        "db_source" : gbk_record.db_source,
        "dblinks" :  gbk_record.dblinks,
        "definition" :  gbk_record.definition,
        "features" :  features_to_json_list(gbk_record.features),
        "gi" :  gbk_record.gi,
        "keywords" :  gbk_record.keywords,
        "locus" :  gbk_record.locus,
        "molecule_type" :  gbk_record.molecule_type,
        "nid" :  gbk_record.nid,
        "organism" :  gbk_record.organism,
        "origin" :  gbk_record.origin,
        "pid" :  gbk_record.pid,
        "primary" :  gbk_record.primary,
        "projects" :  gbk_record.projects,
        "references" :  gbk_record.references,
        "residue_type" :  gbk_record.residue_type,
        "segment" :  gbk_record.segment,
        "sequence" :  gbk_record.sequence,
        "size" :  gbk_record.size,
        "source" :  gbk_record.source,
        "taxonomy" :  gbk_record.taxonomy,
        "topology" :  gbk_record.topology,
        "version" :  gbk_record.version,
        "wgs" :  gbk_record.wgs,
        "wgs_scafld" :  gbk_record.wgs_scafld
    }
    
    return record_d

def features_to_json_list(features_list):
    # Args: features_list: list<feature>
    # where feature is a BioPython SeqFeature Object

    json_features_list = []

    for feat in features_list:
        json_features_list.append({
            "location": location_to_json(feat.location),
            "type": feat.type,
            "location_operator": feat.location_operator,
            "strand": feat.strand,
            "id": feat.id,
            "ref": feat.ref,
            "ref_db": feat.ref_db,
            "qualifiers": feat.qualifiers
    })

    return json_features_list

        

def location_to_json(BioLoc):
    """
    FeatureLocation:
        start
        end
        strand
        ref
        ref_db
    """
    return {
        "start": BioLoc.start,
        "end": BioLoc.end,
        "strand": BioLoc.strand,
        "ref": BioLoc.ref,
        "ref_db": BioLoc.ref_db
    }


def check_headers(headers_d):
    """
    headers_d:
        type -> loc within row list
    """
    headers = list(headers_d)

    recognized_list = [
            "locusId",
            "sysName",
            "type",
            "scaffoldId",
            "begin",
            "end",
            "strand", 
            "name",
            "desc",
            "GC",
            "nTA"
            ]

    # check if headers are recognizable
    for h in headers:
        if h not in recognized_list:
            raise Exception(f"header from TSV file {h} not in recognized_list")

def create_source_feature(seq_id, seq_str):

    feature_loc = SeqFeature.FeatureLocation(1,
                                            len(seq_str),
                                            strand=1,
                                            ref=seq_id
                                            ) 

    new_seq_feat = SeqFeature.SeqFeature(
                                       location=feature_loc, 
                                       type="source",
                                       strand=1,
                                       id=seq_id,
                                       qualifiers={},
                                       ref="")
    return new_seq_feat


def CreateSeqFeatureFromTSVRow(id2seqfeat, TSV_row, type2ix, id2seq, feat_id, special_type=None):
    """


    We create a BioPython SeqFeature object from a TSV row
    Args:
        id2seqfeat: dict mapping ids -> list<SeqFeature>
        TSV_row: list<str>
        type2ix: (dict) type -> index within TSV row e.g.:
                "locusId" : 0
                "sysName" : 1
                "type" : 2
                "scaffoldId",
                "begin",
                "end",
                "strand", 
                "name",
                "desc",
                "GC",
                "nTA"
        special_type: int or None
            If there are weird cases, refer to special_type int



    """

    for v in ["scaffoldId", "begin", "end", "strand", "type", "locusId", "name"]:
        if v not in type2ix:
            raise Exception(f"type value {v} not in headers of TSV input")


    current_seq = id2seq[TSV_row[type2ix['scaffoldId']]]
    qualifier_d = {}
    if special_type == 1:
        qualifier_d = parseSpecialString1(TSV_row[type2ix['qualifiers']])
    else:
        if "desc" in type2ix:
            qualifier_d["product"] = TSV_row[type2ix["desc"]]
        if "name" in type2ix and TSV_row[type2ix["name"]] != "":
            qualifier_d["name"] = TSV_row[type2ix["name"]]
        if "locusId" in type2ix and TSV_row[type2ix["locusId"]] != "":
            qualifier_d["locus_tag"] = TSV_row[type2ix["locusId"]]
   
    strand = 1 if TSV_row[type2ix['strand']] == "+" else -1
    type_num_str = TSV_row[type2ix['type']]
    
    if type_num_str in prog_cfg["type_num_to_str"]:
        feat_type =  prog_cfg["type_num_to_str"][type_num_str]
    else:
        feat_type = "gene"
        

    feature_loc = SeqFeature.FeatureLocation(int(TSV_row[type2ix['begin']]),
                                            int(TSV_row[type2ix['end']]),
                                            strand=strand,
                                            ref=TSV_row[type2ix['scaffoldId']]
                                            ) 

    new_seq_feat = SeqFeature.SeqFeature(
                                       location=feature_loc, 
                                       type=feat_type,
                                       strand=strand,
                                       id=TSV_row[type2ix["locusId"]] + feat_id,
                                       qualifiers=qualifier_d,
                                       ref=TSV_row[type2ix['scaffoldId']])


    #print(new_seq_feat)
    

    return [new_seq_feat, TSV_row[type2ix["scaffoldId"]]]


def parseSpecialString1(s):
    # input s is special string
    # with format x=y z='l s d' f='k' etc..

    arg_d = {}
    
    # First we connect all strings within quotations with |&| 
    ap_l = findOccurrences(s, "'")
    if len(ap_l) % 2 != 0:
        raise Exception("number of ' not divisible by 2 \n locs:" + ",".format(ap_l))

    ap_l.reverse()

    x = createSpecialReplace(s)

    for i in range(int(len(ap_l)/2)):
        end_loc = ap_l[2*i]
        start_loc = ap_l[2*i + 1]
        s = s[:start_loc] + s[start_loc:end_loc].replace(" ", x) + s[end_loc:]

    split_args = s.split(' ')
    for a in split_args:
        k, v = a.split('=')
        if "'" in v:
            v = v[1:-1]
        arg_d[k] = v.replace(x, ' ')

    return arg_d

def createSpecialReplace(s):
    replace = "|&*"
    while replace in s:
        replace += "_<$"

    return replace



def findOccurrences(st, ch):
    return [i for i, letter in enumerate(st) if letter == ch]




def main():
    # This is where the program starts
    args = sys.argv
    if args[-1] not in ["1"]:
        print("Incorrect args. Use the following:\n")
        help_str = "python3 genesTSV2gbk.py genes.tsv genome.fna op_gbk_fp accession_name_1[,accession_name_2] " + \
                    "locus_name [gbk_json_fp] 1"
        print(help_str)
        sys.exit(0)
    else:
        if args[-1] == "1":
            tsv_fp = args[1]
            fasta_fp = args[2]
            op_gbk_fp = args[3]
            accession_list = args[4].split(',')
            locus_name = args[5]
            global prog_cfg 
            prog_cfg = json.loads(open("config.json").read())
            if len(args) == 7:
                ConvertTSVToGBK(tsv_fp, fasta_fp, op_gbk_fp, "", 
                                accession_list=accession_list,
                                locus_base=locus_name
                               )
                print(f"Wrote GBK to {op_gbk_fp}")
            elif len(args) == 8:
                # WRITING JSON FILE AS WELL
                json_fp = args[6]
                ConvertTSVToGBK(tsv_fp, fasta_fp, op_gbk_fp, "", 
                                accession_list=accession_list,
                                locus_base=locus_name,
                                json_fp=json_fp
                               )
                print(f"Wrote GBK to {op_gbk_fp} and JSON file to {json_fp}")
            sys.exit(0)



    return None

if __name__ == "__main__":
    main()


