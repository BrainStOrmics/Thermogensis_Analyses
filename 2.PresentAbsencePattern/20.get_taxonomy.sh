#/bin/bash
export PATH=~/software/ncbi-edirect:$PATH
for i in $(cat homo.species.txt); do 
    echo $i;
    esearch -db taxonomy -query "$i"  | efetch -format xml > tax.d/$i.taxonomy_info.xml
done
