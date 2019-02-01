###

## get isolates from Entrez

from Bio import Entrez
Entrez.email = 'j.gustavsen@gmail.com'

handle = Entrez.einfo()
result = handle.read()

print(result)

handle = Entrez.einfo()
record = Entrez.read(handle)

record.keys()

handle = Entrez.einfo(db="genome")
record = Entrez.read(handle)
print(record)
print(record["DbInfo"]["Description"])
#'PubMed bibliographic record'
print(record["DbInfo"]["Count"])
#'17989604'
print(record["DbInfo"]["LastUpdate"])
#'2008/05/24 06:45'

handle = Entrez.esearch(db="nucleotide", term="Phycodnavirus AND 1000:[Sequence Length]")
record = Entrez.read(handle)
print(record["IdList"])
print(record)

#handle = Entrez.esearch(db="nucleotide", term="Cypripedioideae[Orgn] AND matK[Gene]")
#record = Entrez.read(handle)
#record["Count"]

id_list = record["IdList"]

## get the records, so use the record ids to get the records and then pull out the relevant sections. 
handle = Entrez.efetch("nucleotide", id=",".join(id_list),rettype="gb", retmode="text")
print(handle.read())
