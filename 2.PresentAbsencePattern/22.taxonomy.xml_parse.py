# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 11:37:55 2024

@author: xuyu
"""
from lxml import etree

from bs4 import BeautifulSoup
import os

file = "test.xml"
def parse_file(file):
    # Load and parse the XML file
    with open(file, "r") as fh:
        xml_content = fh.read()
    
    # Parse the XML content
    soup = BeautifulSoup(xml_content, 'xml')
    
    # Extract taxonomy details
    taxonomy_list = []
    for taxon in soup.find_all('Taxon'):
        name = taxon.find('ScientificName').text if taxon.find('ScientificName') else None
        #lineage = taxon.find('Lineage').text if taxon.find('Lineage') else None
        rank = taxon.find('Rank').text if taxon.find('Rank') else "Unknown"
        
        taxonomy_list.append({
            "name": name,
            #"lineage": lineage,
            "rank": rank
        })
    return taxonomy_list
    
    # Print the taxonomy details
    for item in taxonomy_list:
        print(f"name: {item['name']}, Rank: {item['rank']}")

    
def parse_list(in_list, data_dir, out_tsv):
    ifh = open(in_list, 'r')
    ranks = set()
    result_dict = dict()
    result_str = dict()
    name_list = list()
    for name in ifh:
        name = name.strip()
        file = os.path.join(data_dir, name + ".taxonomy_info.xml")
        result = parse_file(file)
        tax_str=";".join(((x['rank'] + ":" + x['name'] for x in result)))
        result_dict[name] = dict((x['rank'], x['name']) for x in result)
        result_str[name] = tax_str
        name_list.append(name)
        ranks.update(x['rank'] for x in result)
    ifh.close()
    
    ranks.remove("no rank")
    ranks.remove("clade")
    rank_list = list(ranks)
    ofh = open(out_tsv, 'w')
    print("name", *rank_list, "tax_str", file=ofh, sep="\t")
    for name in name_list:
        print(name, 
              *[result_dict[name].get(rank, "NA") for rank in rank_list], 
              result_str[name],
              file = ofh, sep="\t")
    ofh.close()
        

if __name__ == '__main__':
    import sys
    if len(sys.argv) - 1 != 3:
        print("Usage: python <script> <genes.list> <data.dir> out.tsv")
        sys.exit(1)
    in_list, data_dir, out_tsv = sys.argv[1:]
    parse_list(in_list, data_dir, out_tsv)
        
