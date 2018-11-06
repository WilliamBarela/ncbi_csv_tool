from Bio import Entrez
from urllib.request import urlopen
from urllib.error import HTTPError
import json
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

def get_csv_stream(email, database, search_term, maximum_returned_items, report=""):
    count, id_list, webenv, query_key = get_query_keys(email, database, search_term, maximum_returned_items)

    url = ""
    if report == "summary":
        url = 'https://www.ncbi.nlm.nih.gov/portal/utils/file_backend.cgi?Db=' + database + '&HistoryId=' + webenv + '&QueryKey=' + query_key + '&Sort=&Filter=all&CompleteResultCount=' + count + '&Mode=file&View=docsumcsv&p$l=Email&portalSnapshot=%2Fprojects%2FSequences%2FSeqDbRelease%401.90&BaseUrl=&PortName=live&FileName=&ContentType=csv'
    elif report == "runinfo":
        url = 'https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&rettype=runinfo&db=' + database + '&WebEnv=' + webenv + '&query_key=' + query_key
    else:
        print("URL cannot be built because allowable report type was not specified!")

    attempts = 1
    while attempts <= 5:
        try:
            # print("Attempt %i of 5" % attempts)
            response = urlopen(url)
            return response
        except HTTPError as err: 
            # error handling in case NCBI server is down.
            if 500 <= err.code <= 599:
                print("Received error from server %s" % err)
                print("Attempt %i of 5" % attempts)
                attempts += 1
                time.sleep(15)
            else:
                raise


def get_csv(email, database, search_term, maximum_returned_items=100000, download=False, report="summary"):
    # FIXME: add additional funciton to write csv to file, include logic here based on "download"
    response = get_csv_stream(email, database, search_term, maximum_returned_items, report)

    csv_string = response.read().decode()
    csv_data = csv.reader(csv_string.splitlines(), delimiter=",", quotechar='"')    # this is a generator
    csv_list = [row for row in csv_data] 

    return csv_list

def get_dict_w_xml(email, database, search_term, maximum_returned_items=100000):
    _, _, webenv, query_key = get_query_keys(email, database, search_term, maximum_returned_items)

    summary_results_xml = Entrez.esummary(db=database,webenv=webenv,query_key=query_key)
    dict_w_xml = Entrez.read(summary_results_xml)

    return dict_w_xml

def print_pretty_json(dict_w_json):
    print(json.dumps(dict_w_json, sort_keys=True, indent=4, separators=(',',':')))


def get_dict_w_json(email, database, search_term, maximum_returned_items=100000, pretty_print=False):
    _, _, webenv, query_key = get_query_keys(email, database, search_term, maximum_returned_items)

    summary_results_json = Entrez.esummary(db=database,webenv=webenv,query_key=query_key,retmode="json")
    raw_json = summary_results_json.read()
    dict_w_json = json.loads(raw_json)
   

    if pretty_print == True:
        print_pretty_json(dict_w_json)

    return dict_w_json
