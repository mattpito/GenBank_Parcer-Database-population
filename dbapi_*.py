#!/usr/bin/python3
"""
	---------Created by Manthos Pitoulias---------------
	This API is performing 2 major tasks. 1) Conneting to db and fetching all the data
                                              2) Converting the data to callable pandas fo the middle layer

    Create locus_entries list for the middle layer to use in dataframe format
    The format as of now 20/03/2021 is for the middle layer to call 8 pieces of 	
    information:
    0:Genebank_accession,
    1:gene_identifier,
    2:chromosomal_location,
    3:protein_product_name,
    4:coding_sequence_boundaries,
    5:complete_DNA_sequence,
    6:aminoacid_sequence,
    7:is_reverse_copliment
    This is the database API - it needs to access the MySQL database for real values.
"""

# Add the directory above to the module path to import the config file
#import sys
#sys.path.insert(0, "../")
import sys
sys.path.insert(0, "../db/")
sys.path.insert(0, "../")

#import dbapi   # Import the database api
#import config  # Import configuration information (if needed)

#import config  # Import configuration information (e.g. database connection)

                       #####################################   FILTERING FUNCTIONS   ##############################################

def dframe_abnormalCDS_clean(df,column_name):
    """Accepts a dataframe as argument.
    Some entries (very few) contain abnormal cds boundaries
    such as (1-7|9) (one nucleotide exon?).
    These entries are cleaned. However if the user wants to keep them, simply don't run the function"""
    #indexes_cds_boundaries = []
    column_name 
    for i in df.index:
        x = df.loc[i,column_name]
        if x.count('-') != x.count('|') + 1:
            df = df[df.index != i]
        else:
            df = df
            #indexes_cds_boundaries.append(i)
    #print(indexes_cds_boundaries)
    return(df)

                      ###################################################################################

def dframe_non3exons_clean(df,column_name):
    """Accepts a dataframe.
    CDS must be in number-number or number-number|number-number..... format
    The combines length of exons should be a multiple of 3 to allow aminoacid conversion. Any other
    value is perhaps the result of alternative splicing or some other biological modification which
    we will not account for. The user can simply not-apply this function if he wishes to or change the CDS
    format"""
    import re

    #indexes_triplet_exons = []
    for i in df.index:
        x = df.loc[i,column_name]
        m1 = re.sub(r'(\d+)-(\d+)',lambda m: str(int(m.group(2)) - int(m.group(1)) + 1),x)
        m2 = (re.findall(r'\d+',m1))
        m2 = sum([int(ele) for ele in m2])
        if m2 % 3 != 0:
            df = df[df.index != i]
        else:
            df = df
            #indexes_triplet_exons.append(i)
    #print(indexes_triplet_exons))
    return(df)

                       #####################################   CONNECT TO DB   ##############################################
import pymysql as mdb
def connect2db():
    db_conn = mdb.connect(host ='****', user = '****', password = '****', db = '****')
    cursor = db_conn.cursor()
    return(cursor)

                        ######################################  CALLABLE FUNCTIONS FOR BL   #############################################
    
def get_all():
    """
    This function will create a list of locus_entries and generate a pandas dataframe   
    out of it named df.
    You can assign it to a vector with the form: x = get_all()
    x will become the df dataframe
    """
    import re
    x = connect2db()
    query = 'SELECT gbk_accession, gene_name, chrom_loc, protein_product, CDS, DNA_seq , translation, complement_status FROM genebank_entries;'
    x.execute(query)
    row = x.fetchall()
    standard_bases = []
    for ele in row:
        atcg_only = re.search(r'[^ATCG]',ele[5])
        if not atcg_only:
            standard_bases.append(ele)
    locus_entries = standard_bases
    #locus_entries.append(['AC001564','ann3','southwest','anna dehydroxylase','20-52|80-202','CGCCCACTCTCCCATTTGTCCGCGCAAGCGGATGCGATGCGATTGCCCGCTAAGATATTCTTACGCGTAACGCAGCTGAGTATTCTACAGAGCTGGCGTACGCGTTGAACACTTCACGGATGATAGGGATTCGGGTAAAGAGCGCGTCATTGGGGGCTTATACAGGCGTAGACTACAACGGGCCCGGCTCAATCACAGCTCGAGCGCCTTGAATAACATACTCATCTCTATACATTCTCGACAGTCTATCGAGCGAGTCGATTACCAACGGGTGCGTTGCAGTTTTAATCTCTCGCCAGCACTGTAATAGCCACCAAGAGACTGATGATAGTCATGGGTGCTGAGCTGAG','QLSTHDYHQSLGGYYSAGERLKLQRTRW-STRSIDCRECIEMSMLFKALEL-LSRARCSLRLYKPPMTRSLPESLSSVKCSTRTPAL-NTQLRYA-EYLSGQSHRIRLRGQMGEWA',1])
    #locus_entries.append(['AB000321','har2','center','harry kinase','19-88','ATTGTTATCTAGGTTGAGATAGTCATGTGCTCACGGAATTTATTGTATGAGTAGTGATTTGAAAGAATTGTTAGTTTGCTGGTTCAAGTAAAAATTCTTT','LLSRLR-SCAHGIYCMSSDLKELLVCWFK-KFF',1])
    #locus_entries.append(['AB000123','man1','northeast','manthase','14-30|70-150','AATTCTTACATTTAAGACCCTAATATCACATCATTAGTGATTAATTGCCACTGCCAAAATTCTGTCCAGAAGCGTTTTAGTTCGCTCCACTAAAGTTGTTTAAAACGACTACCAAATCCGCATGTTAGGGGATTTCTTATTAATTCTTTTATCGTGAGGAACAGCGGATCTTAATGGATGGCCGCAGGTGGTATGGAAGCTAATAGCGCGGGTGAGAGGGTAATCAGCCGTGTTCACCTACACAACGCTAACGGGCGATTCTATAAGATTCCGCATTGCGTCTACTTATAAGATGTCTCAACGGTATCCGCAACTTGTGAAGTGCCTACTATCCTTAAACGCATATCTCGCCCAGTAGCTTCCCAATATGTGAGCATCA','NSYI-DPNITSLVINCHCQNSVQKRFSSLH-SCLKRLPNPHVRGFLINSFIVRNSGS-WMAAGGMEANSAGERVISRVHLHNANGRFYKIPHCVYL-DVSTVSATCEVPTILKRISRPVASQYVSI',0])

    import pandas as pd
    df = pd.DataFrame(locus_entries)
    df.columns = ['Genbank_accession','gene_identifier','chromosomal_location','protein_product_name','coding_sequence_boundaries','complete_DNA_sequence','aminoacid_sequence','is_reverse_complement']
    df = dframe_abnormalCDS_clean(df,'coding_sequence_boundaries')
    df = dframe_non3exons_clean(df,'coding_sequence_boundaries')
    df = df.reset_index(drop=True)
    return(df)
            
                        ###################################################################################
def get_summary(n,par = ''): 
    """
Returns a data frame for the rows in the database that satisfy a certain condition. The columns have the following names and data types:

    Genbank_accession: string
    gene_identifier: string
    chromosomal_location: string
    protein_product_name: string

If it is called with paramter 0 it returns the above information for all the rows in the database. If it is called with parameter 1, 2, 3, and 4, it returns respectively the data that have a specific Genbank accession, gene identifier, chromosomal location, and protein product name.
    """
    if n == 0:
        df_slice = get_all().iloc[:,0:4]
        return(df_slice)
    if n ==1:
        df_slice = get_all().iloc[:,0:4]
        df_slice_ga = df_slice[df_slice['Genbank_accession'] == par]
        return(df_slice_ga)
    if n==2:
        df_slice = get_all().iloc[:,0:4]
        df_slice_gi = df_slice[df_slice['gene_identifier'] == par]
        return(df_slice_gi)
    if n==3:
        df_slice = get_all().iloc[:,0:4]
        df_slice_cl = df_slice[df_slice['chromosomal_location'] == par]
        return(df_slice_cl)
    if n==4:
        df_slice = get_all().iloc[:,0:4]
        df_slice_ppn = df_slice[df_slice['protein_product_name'] == par]
        return(df_slice_ppn)
    
                        ###################################################################################
def get_details(user_input):
    """ returns a dictionary for a single row in the database. The row has the Genbank accession given as input to the function. The keys are the following and the values have the following data types:

complete_DNA_sequence: string
coding_sequence_boundaries: string (for example: "1-92|168-320")
aminoacid_sequence: string
is_reverse_complement: integer (value 0 or 1)
GENEBANK ACCESSION NUMBERS MUST BE UNIQUE FOR THE DICTIONARY TO WORK, otherwise turn to dataframe
    """
    df_slice = get_all().iloc[:,[0,4,5,6,7]]
    df_slice_get_det = df_slice[df_slice['Genbank_accession'] == user_input].iloc[:,[1,2,3,4]]
    d = df_slice_get_det.to_dict('list')
    d = {i:str(j[0]) for i,j in d.items() }
    d = {k:int(v) if v.isdigit() else v for k,v in d.items()}
    return(d)
                        ###################################################################################


