#!/usr/bin/env python

##WRITTEN BY MATTHAIOS PITOULIAS

##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################

##############################################    SOME FUNCTIONS    ##################################################################


def genebank_to_list(filepath,replaceFROM,replaceTO,spliter):
    """Input .gbk filepath, this script creates a single string that is splitted
    into a list of strings based on the splliter argument.
    replaceFROM and replaceTO are the arguments of the replace(a,b) function that help to adjust the string to the
    desired output.
    !!! IN THE COURSEWORK WE MAINLY USE replace(" ","") TO ELIMINATE SPACES!!! and then split on '//LOCUS'
    This format allows to iterate through
    each entry and use regex to find the first or all matches per a single accession (spliter='//LOCUS')
    While // is at the end of each entry, it also exists in text inside some entries. Therefore //LOCUS
    is used. The first entry will start with LOCUS and all subsequent ones will start from the LOCUS ACCESSION
    (e.g., AB12345). By eliminating the  aberrant empty spaces we also increase the speed manytimes fold"""
    gbk = ""
    with open(f'{filepath}') as f:
        gb = f.read().splitlines()
        for line in gb:
            gbk += line
            
    gbk_nospace = gbk.replace(replaceFROM,replaceTO)
    gbk_aslist = list(gbk_nospace.split(spliter))
    return(gbk_aslist)

##########################################################################################################

def clean_unmatched(gbk_aslist,x):
    """Removing all unmatched entries, e.g., entries without CDS. gbk_aslist = name of the list variable.\
    x = regex exrpession in rhe form r'' """
    import re
    gbk_aslist_toremove = []
    for entry in gbk_aslist:
        first_match = re.search(x,entry)  #remove entries without CDS
        if not first_match:
            #print(LOCUS.group())
            gbk_aslist_toremove.append(entry)
    for x in gbk_aslist_toremove:
        if x in gbk_aslist:
            gbk_aslist.remove(x)
            
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################

##############################################    DATA EXTRACTION    ##################################################################


#Working example: From the downloaded X chromosome genebankfile. The file was downloaded and unzipped via the cmd with gunzip filename.gz at the appropriate dictionary.

#access home directory as if in bash
from os.path import expanduser
home = expanduser('~')
#call function genebank_to_list assigned to z
z = genebank_to_list(filepath = home+'/biocomp2/x',replaceFROM=' ',replaceTO='',spliter='//LOCUS') #// is also found in-entry and messes
#the split. Therefore I split at //LOCUS which is only found in entry junctions (assuming no whitespace string).
#call function to clean z from unmatched CDS entries
clean_unmatched(z,x=r'CDS\w+[\(\w+]*[\(*\d+\.+\d+\,]+\)|CDS\d+\.+\d+')
#check number of entries to currently work with
#len(z) # 733 Entries. All downstream call will have to be 733.

                         #########GET CDS############

#First filtering will be based on CDS. Only matches that meet the CDS criteria will be displayed
#Gene identifiers that are not present in an entry are allowed
#however CDS is mandatory and the first filtering will be based on entries that match the CDS criteria
#Find first CDS match on condition (not allowing weird matches <> or fusions) and insert in list
import re
CDS_list = []
for entry in z:
    CDS = re.search(r'CDS\w+[\(\w+]*[\(*\d+\.+\d+\,]+\)|CDS\d+\.+\d+',entry)  
    if CDS:
        CDS_list.append(CDS.group())

                         #########GET ACCESSION############

#Find ACCESSION numbers and insert them in a list
ACCESSION_list = []
for entry in z:
    ACCESSION = re.search(r'ACCESSION([A-Z]+\d+)',entry)
    if ACCESSION:
        ACCESSION_list.append(ACCESSION.group(1))

                         #########GET Chromosomal location############
#Find chromosomal MAP locations and insert them in a list
MAP_list = []
for entry in z:
    MAP = re.search(r'/map="(.+?)"',entry)
    if MAP:
        MAP_list.append(MAP.group(1)) 
        #print(MAP.group())
    if not MAP: #Not every gene is properly mapped.
        MAP_list.append('N/A')    #if no map information provided, create a NOT APPLICABLE entry.

                        ############# GET gene name######################

#Find gene identifier (GI) entries and insert them in a list
GI_list = []
for entry in z:
    GI = re.search(r'/gene="(.+?)"',entry)
    if GI:
        GI_list.append(GI.group(1)) 
        #print(GI.group())
    if not GI: #Not every gene is properly identified.
        GI_list.append('N/A')    #if no GI provided, create a NOT APPLICABLE entry.


                       ############# GET DNA SEQ############################
#Find DNA sequences for entries and insert them in a list
DNAseq_list = []
for entry in z:
    DNAseq = re.search(r'ORIGIN(.+)',entry)
    if DNAseq:
        DNAseq_bases = re.sub(r"[0-9]", "", DNAseq.group(1))
        DNAseq_list.append(DNAseq_bases)
DNAseq_list = [k.upper() for k in DNAseq_list]

                        ############GET AMINO ACID ######################
#Find translation (amino acid sequence) of entries and insert them in a list
AAseq_list = []
for entry in z:
    AAseq = re.search(r'/translation="(.+?)"',entry)
    if AAseq:
        AAseq_list.append(AAseq.group(1)) 
        #print(aa_seq.group())
    if not AAseq: #This can happen for example in the case of pseudogenes
        AAseq_list.append('N/A')

                     ################# GET PROTEIN PRODUCTS ###############
# Slightly altered. This time we want to have whitespaces to get back meaningful text
#Extract act protein products for the filtered entries based on the 733 accession numbers
c = genebank_to_list(filepath = home+'/biocomp2/x',replaceFROM='',replaceTO='',spliter='//LOCUS')
gg = [' '.join(x.split()) for x in c]
protprod_pool = []
for ele in ACCESSION_list:
    for i,v in enumerate(gg):
        m1=re.search('ACCESSION +' + ele,gg[i])
        if m1:
            protprod_pool.append(v)
PROTprod_list = []
for x in protprod_pool:
    PROTprod = re.search(r'/product="(.+?)"',x)
    if PROTprod:
        PROTprod_list.append(PROTprod.group(1)) 
        #print(aa_seq.group())
    if not PROTprod: #This can happen for example in the case of pseudogenes
        PROTprod_list.append('N/A')

                     ################# ENCODE normal or complement CDS to 0,1#######################
is_reverse_complement = []
for element in CDS_list:
    match = re.search(r'complement',element)
    if match:
        is_reverse_complement.append(1)
    else:
        is_reverse_complement.append(0)

	###################### CLEAN CDS to match requested format for the BLAYER##################
CDS_clean1= []
for x in CDS_list:
    clean1 = re.sub(r'[a-zA-Z]','',x)
    CDS_clean1.append(clean1)
CDS_clean2 = [w.replace('(', '') for w in CDS_clean1]
CDS_clean3 = [w.replace(')', '') for w in CDS_clean2]
CDS_clean4 = [w.replace('..','-') for w in CDS_clean3]
CDS_formatted = [w.replace(',','|') for w in CDS_clean4]

                       ############# Conjoint all to a list of lists #########################
entry_locus = [list(x) for x in zip(ACCESSION_list,GI_list,MAP_list,PROTprod_list,CDS_formatted,DNAseq_list,AAseq_list,is_reverse_complement)]
#print(entry_locus[0:3]) test

##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################

##############################################    ACCESS DB AND POPULATE    ##################################################################
#ACCESS
import pymysql as mdb

db_conn = mdb.connect('localhost', user = '***', password = '****', db = '****')
cursor = db_conn.cursor()

#POPULATE. TO POPULATE UNHASH AND RUN. DB IS  ALREADY POPULATED
#cursor.execute('SHOW TABLES')
#cursor.executemany("INSERT INTO genebank_entries(gbk_accession, gene_name, chrom_loc, protein_product, CDS, DNA_seq , translation, complement_status) VALUES (%s, %s, %s, %s, %s, %s, %s, %s)", entry_locus)
#db_conn.commit()










