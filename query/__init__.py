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

def open_url(url):
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

def get_report_urls(email, database, search_term, maximum_returned_items):
    count, id_list, webenv, query_key = get_query_keys(email, database, search_term, maximum_returned_items)

    summary_url = 'https://www.ncbi.nlm.nih.gov/portal/utils/file_backend.cgi?Db=' + database + '&HistoryId=' + webenv + '&QueryKey=' + query_key + '&Sort=&Filter=all&CompleteResultCount=' + count + '&Mode=file&View=docsumcsv&p$l=Email&portalSnapshot=%2Fprojects%2FSequences%2FSeqDbRelease%401.90&BaseUrl=&PortName=live&FileName=&ContentType=csv'
    run_url = 'https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&rettype=runinfo&db=' + database + '&WebEnv=' + webenv + '&query_key=' + query_key

    return (summary_url,run_url)



def get_csv_list(response):
    csv_string = response.read().decode()
    csv_data = csv.reader(csv_string.splitlines(), delimiter=",", quotechar='"')    # this is a generator
    csv_list = [row for row in csv_data] 

    return csv_list

def get_csv_lists(email, database, search_term, maximum_returned_items=100000):
    """
        This function has been added to ensure that only one post request is made to the NCBI API.
        This is not only to ensure efficiency of API calls, but also to ensure consistency of the reports.
        If this is not done this way, there will be two webenv calls, one for each report.
        If either contains a more recent dataset, this could break other code.
        As a result, it is best to use one webenv for both the summary and the run reports.
    """
    summary_url, run_url = get_report_urls(email, database, search_term, maximum_returned_items)

    summary_csv_list = get_csv_list(open_url(summary_url))
    run_csv_list = get_csv_list(open_url(run_url))

    return (summary_csv_list, run_csv_list)

def save_to_csv(csv_list=None, download_directory=None):
    # FIXME: add additional funciton to write csv to file, include logic here based on "download"
    print("Use os from the standard library to set the download directory; use csv.write() to write the csvs to the specified location.")

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
