from Bio import Entrez
from urllib.request import urlopen
from urllib.error import HTTPError
import csv

def get_ncbi_ids(email, database, search_term, maximum_returned_items):
    Entrez.email = email 
    search_results_xml = Entrez.esearch(db=database, term=search_term, retmax=maximum_returned_items)
    ncbi_dict = Entrez.read(search_results_xml)
    count = ncbi_dict["Count"]
    id_list = ncbi_dict['IdList']

    return (count, id_list) 

def get_search_cache_keys(database, id_list):
    search_post_request = Entrez.read(Entrez.epost(database,id=",".join(id_list)))
    webenv = search_post_request["WebEnv"]
    query_key = search_post_request["QueryKey"]

    return (webenv, query_key)

def get_query_keys(email, database, search_term, maximum_returned_items):
    count, id_list = get_ncbi_ids(email, database, search_term, maximum_returned_items)
    webenv, query_key = get_search_cache_keys(database, id_list)

    return (count, id_list, webenv, query_key)

def get_csv_stream(email, database, search_term, maximum_returned_items):
    count, id_list, webenv, query_key = get_query_keys(email, database, search_term, maximum_returned_items)

    attempts = 1
    while attempts <= 5:
        try:
            url = 'https://www.ncbi.nlm.nih.gov/portal/utils/file_backend.cgi?Db=' + database + '&HistoryId=' + webenv + '&QueryKey=' + query_key + '&Sort=&Filter=all&CompleteResultCount=' + count + '&Mode=file&View=docsumcsv&p$l=Email&portalSnapshot=%2Fprojects%2FSequences%2FSeqDbRelease%401.90&BaseUrl=&PortName=live&FileName=&ContentType=csv'
            print("Attempt %i of 5" % attempts)
            break
        except HTTPError as err: 
            # error handling in case NCBI server is down.
            if 500 <= err.code <= 599:
                print("Received error from server %s" % err)
                print("Attempt %i of 5" % attempts)
                attempts += 1
                time.sleep(15)
            else:
                raise

    response = urlopen(url)
    return response

def get_csv(email, database, search_term, maximum_returned_items=100000, download=False):
    # FIXME: add additional funciton to write csv to file, include logic here based on "download"
    response = get_csv_stream(email, database, search_term, maximum_returned_items)

    csv_string = response.read().decode()
    csv_data = csv.reader(csv_string.splitlines(), delimiter=",", quotechar='"')    # this is a generator
    csv_list = [row for row in csv_data] 

    return csv_list

def get_dict_w_xml(email, database, search_term, maximum_returned_items=100000):
    _, _, webenv, query_key = get_query_keys(email, database, search_term, maximum_returned_items)

    summary_results_xml = Entrez.esummary(db=database,webenv=webenv,query_key=query_key)
    dict_w_xml = Entrez.read(summary_results_xml)

    return dict_w_xml
