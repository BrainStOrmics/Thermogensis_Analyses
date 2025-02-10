# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 17:55:07 2025

@author: xuyu
"""
import pandas as pd
import json
import sys
import os

#js_file = "test.pair_align.json"

def js2df(js_file, df_file):
    with open(js_file) as f:
        js_data = json.load(f)
        
    # with open('test.fmt.json', 'w') as file:
    #     json.dump(js_data, file, indent=4)  # Adjust the indent level as needed
    
    #data = js_data
    data = js_data['data'][0]['homologies']
    #type(data)
    
    #data.keys()
    #len(data)
    
    df = pd.json_normalize(data)
    df.to_csv(df_file, index=False)
    #df
    #df[1]
    
def walk(in_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    for filename in os.listdir(in_dir):
        if filename.endswith('.json'):
            print(filename)
            in_file = os.path.join(in_dir, filename)
            basename = os.path.basename(filename)
            out_file = os.path.join(out_dir, basename + ".csv")
            js2df(in_file, out_file)
            
if __name__ == '__main__':
    if len(sys.argv) - 1 !=2 :
        print("convert homology rest API json to df, from folder to folder")
        print("Usage: python <script> <in.json.dir> <out.csv.dir>")
        sys.exit(1)
    in_dir, out_dir = sys.argv[1:]
    walk(in_dir, out_dir)
