import json
from Bio import Entrez
from urllib.request import urlopen
import csv

def get_ncbi_ids(email, database, search_term, maximum_returned_items):
    Entrez.email = email 
    search_results_xml = Entrez.esearch(db=database, term=search_term, retmax=maximum_returned_items)
    ncbi_dict = Entrez.read(search_results_xml)
    count = ncbi_dict["Count"]
    id_list = ncbi_dict['IdList']

    return (count, id_list) 

def get_search_cache_keys(email, database, search_term, maximum_returned_items):
    count, id_list = get_ncbi_ids(email, database, search_term, maximum_returned_items)

    search_post_request = Entrez.read(Entrez.epost(database,id=",".join(id_list)))
    webenv = search_post_request["WebEnv"]
    query_key = search_post_request["QueryKey"]

    return (count, id_list, webenv, query_key)

def get_csv(email, database, search_term, maximum_returned_items=100000, download=False):
    count, id_list, webenv, query_key = get_search_cache_keys(email, database, search_term, maximum_returned_items)

    url = 'https://www.ncbi.nlm.nih.gov/portal/utils/file_backend.cgi?Db=sra&HistoryId=' + webenv + '&QueryKey=' + query_key + '&Sort=&Filter=all&CompleteResultCount=' + count + '&Mode=file&View=docsumcsv&p$l=Email&portalSnapshot=%2Fprojects%2FSequences%2FSeqDbRelease%401.90&BaseUrl=&PortName=live&FileName=&ContentType=csv'
    response = urlopen(url)

    csv_string = response.read().decode()
    csv_data = csv.reader(csv_string.splitlines(), delimiter=",", quotechar='"')    # this is a generator, so you need to do the next line
    csv_list = [row for row in csv_data] 

    return csv_list
