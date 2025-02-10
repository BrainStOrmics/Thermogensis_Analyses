import requests
import xml.etree.ElementTree as ET
import pandas as pd
import time
import os

# 读取基因名称
# with open('gene_names.txt', 'r') as f:
#     gene_names = [line.strip() for line in f.readlines()]
gene_df = pd.read_csv("homo.TMEM41B.txt", delimiter="\t")

server = "https://rest.ensembl.org"

# 遍历每个基因名称，获取同源基因
for i, rec in enumerate(gene_df.itertuples()):
    species = rec.target_species
    gene_id = rec.target_gene_stable_id
    prefix = f"{i}.{species}.{gene_id}"
    # ext = f"/homology/symbol/human/{gene_name}?format=condensed;type=orthologues"
    # ext = f"/homology/symbol/{species}/{gene_id}?content-type=application/json"
    ext = f"/homology/id/{species}/{gene_id}?type=orthologues"
    print(prefix)
    
    # 重试机制
    for attempt in range(5):  # 尝试最多5次
        try:
            #r = requests.get(server + ext, headers={"Content-Type": "text/xml"})
            r = requests.get(server + ext, headers={"Content-Type": "application/json"})
            r.raise_for_status()  # 检查请求是否成功
            with open(os.path.join("pair_align.json", prefix+".json"), 'w') as ofh:
                ofh.write(r.text)

            
            
#            root = ET.fromstring(r.text)
#            data = []
#            
#            for homology in root.findall('.//homologies'):
#                info = {
#                    'Homology ID': homology.attrib['id'],
#                    'Method Link Type': homology.attrib['method_link_type'],
#                    'Protein ID': homology.attrib['protein_id'],
#                    'Species': homology.attrib['species'],
#                    'Taxonomy Level': homology.attrib['taxonomy_level'],
#                    'Type': homology.attrib['type']
#                }
#                data.append(info)
#
#            # 如果找到同源基因，保存到CSV文件
#            if data:
#                df = pd.DataFrame(data)
#                output_file = f"{gene_name}_homologies.csv"
#                df.to_csv(output_file, index=False)
#                print(f"同源基因数据已保存到 {output_file}")
#            break  # 如果成功则退出重试循环
        
        except requests.exceptions.RequestException as e:
            print(f"请求失败: {e}")
            time.sleep(2 ** attempt)  # 指数退避

    time.sleep(1)  # 每次请求间隔1秒

print("所有基因的同源基因数据已处理完毕。")
