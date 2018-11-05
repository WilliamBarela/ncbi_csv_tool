import json
from Bio import Entrez
from urllib.request import urlopen
import csv

Entrez.email = "william.barela@ttu.edu"
handle = Entrez.esearch(db="sra", term="heterodera[Organism]",retmax=1000)
ncbi_dict = Entrez.read(handle)
count = ncbi_dict["Count"]
id_list = ncbi_dict['IdList']

search_post_request = Entrez.read(Entrez.epost("sra",id=",".join(id_list)))
webenv = search_post_request["WebEnv"]
query_key = search_post_request["QueryKey"]

url = 'https://www.ncbi.nlm.nih.gov/portal/utils/file_backend.cgi?Db=sra&HistoryId=' + webenv + '&QueryKey=' + query_key + '&Sort=&Filter=all&CompleteResultCount=' + count + '&Mode=file&View=docsumcsv&p$l=Email&portalSnapshot=%2Fprojects%2FSequences%2FSeqDbRelease%401.90&BaseUrl=&PortName=live&FileName=&ContentType=csv'
response = urlopen(url)

csv_string = response.read().decode()
csv_data = csv.reader(csv_string.splitlines(), delimiter=",", quotechar='"')    # this is a generator, so you need to do the next line
csv_list = list(csv)

