import json
import requests
import pandas as pd
import sys

species = sys.argv[1]

post_experiment = {
    "query": {
        "bool": {
            "must": [
                {"match": {"standardMet" : "FAANG"}},
                {"exists": {"field" : "ChIP-seq DNA-binding"}}
            ]
        }
    }
}

post_input = {
    "query": {
        "bool": {
            "must": [
                {"match": {"standardMet" : "FAANG"}},
                {"exists": {"field" : "ChIP-seq input DNA"}}
            ]
        }
    }
}

post_files = {
    "query": {
        "bool": {
            "must": [
                {"match": {"experiment.standardMet" : "FAANG"}},
                {"match": {"experiment.assayType" : "ChIP-seq"}},
                {"match": {"species.text" : species }}
            ]
        }
    }
}

def get_from_faang_api(url, post_data, verbose = False):
    response = requests.post(url, json=post_data)
    result = response.json()
    if verbose: print(f"Number of hits: {result['hits']['total']}")
    return result['hits']['hits']

url = 'http://data.faang.org/api/experiment/_search/?size=20000'
experiments = get_from_faang_api(url, post_data=post_experiment, verbose=True)
input_dna = get_from_faang_api(url, post_data=post_input, verbose=True)

files = get_from_faang_api('http://data.faang.org/api/file/_search/?size=20000', post_data=post_files, verbose=True)

inputs = {d["_id"]: "" for d in input_dna} 
exps = {d["_id"]: d["_source"]["ChIP-seq DNA-binding"]["controlExperiment"] for d in experiments}
all_chip = inputs
all_chip.update(exps)  

paired = []
single = []

for f in files:
    fname = f["_id"]
    print(fname)
    furl = "ftp://" + f["_source"]["url"]
    fexp = f["_source"]["experiment"]["accession"]
    if fname[-2:] != "_2":
        design = {"group" : fexp, "replicate": 1}
        design["fastq_1"] = furl
        design["fastq_2"] = ""
        design["antibody"] = f["_source"]["experiment"]["target"]
        design["control"] = all_chip[fexp]
        if fname[-2] == "_":
            design["fastq_2"] = furl[:-10] + "2" + furl[-9:]
            paired.append(design)
        else: 
            single.append(design)

# PAIRED

if paired:
    df = pd.DataFrame(paired)
    df.replace("input DNA", "", inplace=True)
    
    print(df)
    
    with open(f"paired.csv", 'w') as f:
        f.write(df.to_csv(index=False))
    
    nExpsPerControl = df['control'].value_counts()
    controls = nExpsPerControl.index
    nGroups = df['group'].value_counts()
    allGroups = nGroups.index
    
    for inputDNA in controls:
        if inputDNA in allGroups:
            csv = df[df.control == inputDNA]
            csv = csv.append(df[df.group == inputDNA])
            with open(f"paired_{inputDNA}.csv", 'w') as f:
                f.write(csv.to_csv(index=False))

# SINGLE

if single:
    df = pd.DataFrame(single)
    df.replace("input DNA", "", inplace=True)
    
    print(df)
    
    with open(f"single.csv", 'w') as f:
        f.write(df.to_csv(index=False))
    
    nExpsPerControl = df['control'].value_counts()
    controls = nExpsPerControl.index
    nGroups = df['group'].value_counts()
    allGroups = nGroups.index
    
    for inputDNA in controls:
        if inputDNA in allGroups:
            csv = df[df.control == inputDNA]
            csv = csv.append(df[df.group == inputDNA])
            with open(f"single_{inputDNA}.csv", 'w') as f:
                f.write(csv.to_csv(index=False))

