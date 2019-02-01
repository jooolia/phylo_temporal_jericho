## script to downloap cyanophage genomes from NCBI
## then want to compare the gp20 regions of them. 

from Bio import Entrez
from datetime import date
Entrez.email = "jgustavsen@eos.ubc.ca"

today = date.today()
print today

handle = Entrez.einfo(db="pubmed")
result = Entrez.read(handle)
handle.close()

print result["DbInfo"].keys()

for field in result["DbInfo"]["FieldList"]:
    print("%(Name)s, %(FullName)s, %(Description)s" % field)

search_handle = Entrez.esearch(db="pubmed", term="Picornavirales[ORGN] AND marine", usehistory="y")
result = Entrez.read(search_handle)
search_handle.close()
print result["Count"]
print result["IdList"]
print result.keys()
print result


list_of_ids = result["IdList"]
print list_of_ids


gi_list = result["IdList"]
count = int(result["Count"])
print count, len(gi_list)
webenv = result["WebEnv"]
query_key = result["QueryKey"] 


batch_size = 3
out_handle = open("../data/MPL_sequence_references/picornavirales_environmental_downloaded_"+ str(today) + ".fasta", "w")
print(out_handle)
for start in range(0, count, batch_size):
    end = min(count, start + batch_size)
    print("Fetching records %i thru %i" % (start + 1, end))
    fetch_handle = Entrez.efetch(db="pubmed", rettype="gb", retmode="text",
                                 retstart=start, retmax=batch_size,
                                 webenv=webenv, query_key=query_key)
    data = fetch_handle.read()
    fetch_handle.close()
    out_handle.write(data)
out_handle.close()