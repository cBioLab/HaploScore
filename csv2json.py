import csv
import json

haplolist = []
seqdic = {}

with open('output.csv', 'r') as f:
  reader = csv.reader(f)

  for row in reader:
    if row == ['0'] :
      break
    haplolist.append(row)

  for row in reader:
    if row == ['0'] :
      break
    seqdic[row[0]]=row[1];

f = open('output.json', 'w')

for haplo in haplolist:
  seqstr = ""
  for hid in haplo:
    seqstr += seqdic[hid]
  result = {}
  result['sequence'] = seqstr
  result['name'] = 'hoge'
  result['path'] = {}
  result['path']['name'] = 'fuga'
  result['path']['mapping'] = []
  rank = 1
  for hid in haplo:
    mappingdict = {}
    mappingdict['position'] = {}
    mappingdict['position']['node_id'] = hid
    mappingdict['edit'] = []
    lendict = {}
    lendict['from_length'] = len(seqdic[hid])
    lendict['to_length'] = len(seqdic[hid])
    mappingdict['edit'].append(lendict)
    mappingdict['rank'] = str(rank)
    result['path']['mapping'].append(mappingdict)
    rank += 1
  json.dump(result, f)
  f.write("\n")
