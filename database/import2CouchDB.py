# reads json files into couchdb
# takes a directory and / or a single file
# does not error check for existance of duplicate documents
import sys
import os
import couchdb
import json
import glob
import argparse


parser = argparse.ArgumentParser(description="Import files into couchdb")
parser.add_argument('couchdb_url',help="location of the couch database")
parser.add_argument('database_name',help="name of the database")
parser.add_argument('-d','--source_dir',help="directory containing json files to import")
parser.add_argument('-f','--input_file',help="a single file to import")
args = parser.parse_args()

print "couchDB at",args.couchdb_url
database_url=args.couchdb_url+'/'+args.database_name
print "database at",database_url

# open connection to server and define database
couch = couchdb.Server()
db = couch[args.database_name]

# get the list of files to import
jsons=[]
if (args.source_dir):
	os.chdir(args.source_dir)
	for file in glob.glob("*.json"):
		print file
		jsons.append(file)

if (args.input_file):
	print args.input_file
	jsons.append(file)

nfiles = len(jsons)
print "putting",nfiles,"documents into couchDB"

# get n UUIDs, where n=number of json documents
# put them in a list
uuids = couch.uuids(nfiles)

# add {_id:uuid} to contents of each file
# and load into couchdb
# todo: should check to see if file exists / is modified
for doc in jsons:
    print doc
    data = json.load(open(doc))
    uuid = uuids.pop()
    data['_id'] = uuid
    #docid,docrev=db.save(data)
    #print docid,docrev





    
