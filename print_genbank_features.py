
import textwrap

"""
General info


        3.4.12.3 Feature Qualifiers
          Qualifiers provide additional information about features. They take
          the form of a slash (/) followed by a qualifier name and, if
          applicable, an equal sign (=) and a qualifier value. Feature
          qualifiers begin at column 22.
"""

def create_first_line(feat_type, feature_loc):
    """
        3.4.12.1 Feature Key Names
          The first column of the feature descriptor line contains the feature
          key. It starts at column 6 and can continue to column 20

    Args:
        feat_type: (str)
        feature_loc: BioPython's FeatureLocation Object
    """
    if len(feat_type) > 14:
        feat_type = feat_type[:14]
    feat_loc_str = get_genbank_loc_string(feature_loc)
    # Adding first line to feature string 
    first_line = " "*5 + feat_type + (16 - len(feat_type))*" " + feat_loc_str + \
                    (58 - len(feat_loc_str))*" " + "\n"
    return first_line


def get_genbank_loc_string(feature_loc):
    """
        3.4.12.2 Feature Location
          The second column of the feature descriptor line designates the
          location of the feature in the sequence. The location descriptor
          begins at position 22. Several conventions are used to indicate
          sequence location.
    Args: 
        feature_loc (Feature in BioPython's FeatureLocation object)
        start (int) (zero based, should add one)?
        end (int) (zero based, should add one)?
        strand: int (1, -1, 0)
    """
    
    if feature_loc.strand not in [1, -1]:
        raise Exception(f"Does not recognize feature strand {feature_loc.strand}")

    loc_str = ""
    if feature_loc.start != feature_loc.end:
        loc_str = f"{feature_loc.start}..{feature_loc.end}"
    else:
        loc_str = str(feature_loc.start)
    if feature_loc.strand == -1:
        loc_str = "complement(" + loc_str + ")"

    if len(loc_str) > 58:
        raise Exception(f"location string length is too long (line above 80 chars). Length:{len(loc_str)}")

    return loc_str


def create_qualifier_str(qualifier_name, qualifier_value):
    """
        3.4.12.3 Feature Qualifiers
          Qualifiers provide additional information about features. They take
          the form of a slash (/) followed by a qualifier name and, if
          applicable, an equal sign (=) and a qualifier value. Feature
          qualifiers begin at column 22.

        We use the library textwrap
    """

    wrap_text  =  f'/{qualifier_name}="{qualifier_value}"'
    text_list = textwrap.wrap(wrap_text, width=58)
    qual_str = ""
    for l in text_list:
        qual_str += " "*21 + l + "\n"

    return qual_str
    

def get_gbk_str(gbk_rec):
        """Provide a GenBank formatted output option for a Record.

        The objective of this is to provide an easy way to read in a GenBank
        record, modify it somehow, and then output it in 'GenBank format.'
        We are striving to make this work so that a parsed Record that is
        output using this function will look exactly like the original
        record.

        Much of the output is based on format description info at:

        ftp://ncbi.nlm.nih.gov/genbank/gbrel.txt
        """
        output = gbk_rec._locus_line()
        output += gbk_rec._definition_line()
        output += gbk_rec._accession_line()
        output += gbk_rec._version_line()
        output += gbk_rec._project_line()
        output += gbk_rec._dblink_line()
        output += gbk_rec._nid_line()
        output += gbk_rec._pid_line()
        output += gbk_rec._keywords_line()
        output += gbk_rec._db_source_line()
        output += gbk_rec._segment_line()
        output += gbk_rec._source_line()
        output += gbk_rec._organism_line()
        for reference in gbk_rec.references:
            output += str(reference)
        output += gbk_rec._comment_line()
        output += gbk_rec._features_line()
        for feature in gbk_rec.features:
            output += get_feature_string(feature)
        output += gbk_rec._base_count_line()
        output += gbk_rec._origin_line()
        output += gbk_rec._sequence_line()
        output += gbk_rec._wgs_line()
        output += gbk_rec._wgs_scafld_line()
        output += gbk_rec._contig_line()
        output += "//"
        return output




def get_feature_string(feat):
    """
    feat is a BioPython SeqFeature object
    """

    crnt_ft_str = ""

    # Adding first line to feature string (func from print_genbank_features) 
    crnt_ft_str +=  create_first_line(feat.type, feat.location)
    feat_qualifiers_d = feat.qualifiers
    for qlf in feat_qualifiers_d.keys():
        # (func from print_genbank_features)
        crnt_ft_str += create_qualifier_str(qlf, feat_qualifiers_d[qlf])

    return crnt_ft_str



def get_all_features_str(gbk_record):
    """
    We print out features in the style of a GenBank file.
    https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html#FeaturesSourceB

    Below, columns are numbered 1-based
    Info from file gbrel.txt (3.4.12)
        3.4.12.1 Feature Key Names

          The first column of the feature descriptor line contains the feature
          key. It starts at column 6 and can continue to column 20

        3.4.12.2 Feature Location

          The second column of the feature descriptor line designates the
          location of the feature in the sequence. The location descriptor
          begins at position 22. Several conventions are used to indicate
          sequence location.

        3.4.12.3 Feature Qualifiers

          Qualifiers provide additional information about features. They take
          the form of a slash (/) followed by a qualifier name and, if
          applicable, an equal sign (=) and a qualifier value. Feature
          qualifiers begin at column 22.



    features_list: list<seqfeat>
                        seqfeat: SeqFeature.SeqFeature from BioPython
                        SeqFeature.SeqFeature at top of file
    
    """

    genbank_features_str = ""
    # getting the features list of this record
    features_list = gbk_record.features
    for feat in features_list:
        genbank_features_str += get_feature_string(feat)

    return genbank_features_str

