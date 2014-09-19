import sys
import json
from datetime import *

def encode_ott_json(uid,name,rank):
    x = {"@context": {"name": "http://schema.org/name","prov": "http://www.w3.org/ns/prov#","oa": "http://www.w3.org/ns/oa#"},"@type": "oa:Annotation","oa:annotatedBy": {"@type": "prov:Entity","name": "blackrim"}}
    x["oa:annotatedAt"] = str(datetime.now())
    x["oa:hasTarget"] = { "@type":"node" ,"included_ids":[uid], "error_checks":[], "warning_checks":[]}
    x["oa:hasBody"] = {"@type" : "taxonomy label", "@id" : "IRI","name":name,"rank":"","source":"ott","unique id":uid}
    return x

"""
this presumes that there will be the file
taxonomy.tsv
within the ott_dir
"""
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print "python "+sys.argv[0]+" tax_file"
        sys.exit(0)
    tax_file = open(sys.argv[1],"r")
    tax_file.readline()
    of = open("ott_taxonomy_annotations.json","w")
    for i in tax_file:
        spls = i.strip().split("\t|")
        uid = spls[0].strip()
        name = spls[2].strip()
        rank = spls[3].strip()
        x = encode_ott_json(uid,name,rank)
        json.dump(x,of)
        of.write("\n")
    of.close()
    tax_file.close()
