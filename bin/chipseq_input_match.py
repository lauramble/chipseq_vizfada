#!/usr/bin/env python

import sys
import os
from ftplib import FTP
import time

"""
import subprocess

installed = set(sys.modules.keys())
required = {"flatten-json", "requests"}
missing = required - installed

if missing:
    python = sys.executable
    for pkg in missing:
        subprocess.check_call(
            [python, '-m', 'pip', 'install', pkg], stdout=subprocess.DEVNULL)
"""

import requests
import flatten_json as fj
import pandas as pd

ONTOLOGY_URL = "http://purl.obolibrary.org/obo/"
KEEP_COLUMNS = ["derivedFrom","experimentalProtocol","extractionProtocol",
                "files","libraryPreparationDate","libraryPreparationLocation",
                "material","organism_spec","organization","paperPublished",
                "project","publishedArticles","releaseDate","RNA-seq",
                "ChIP-seq DNA-binding","runs","sampleStorage","sampleStorageProcessing",
                "samplingToPreparationInterval","secondaryProject_exp",
                "sequencingDate","sequencingLocation","species","study",
                "submission","updateDate","accession","assayType","biosampleId",
                "cellSpecimen","cellType","experiment","poolOfSpecimens"]


def check_paired_reads(s):
    if len(s.split("_"))>1:
        return True
    else:
        return False


def reformat(arg, level=0):
    #print("{} - {}".format(level, arg))

    if type(arg)==list:
        l = []
        for i in arg:
            l.append(reformat(i, level=level+1))
        return l

    if type(arg)==dict:
        if "ontologyTerms" in arg.keys():
            #print("{} {}".format(((type(arg["ontologyTerms"])==str) and (ONTOLOGY_URL not in arg["ontologyTerms"])), arg))
            if (type(arg["ontologyTerms"])==str):
                if (ONTOLOGY_URL not in arg["ontologyTerms"]):
                    arg["ontologyTerms"] = ONTOLOGY_URL + arg["ontologyTerms"]
                return arg
            elif (type(arg["ontologyTerms"])==list):
                l = []
                for term in arg["ontologyTerms"]:
                    if (ONTOLOGY_URL not in term):
                        l.append(ONTOLOGY_URL + term)
                return l
            else:
                return reformat(arg, level=level+1)
        else:
            d = {}
            for k,v in arg.items():
                #print("{}: {}".format(k, v))
                d.setdefault(k, reformat(v, level=level+1))
            return d

    return arg


def json_to_df(json, flatten=True, root_keys_to_ignore=set()):
    ignore = root_keys_to_ignore
    if (flatten):
        for i in range(len(json)):
            ignore = ignore.union({key for key in json[i]["_source"].keys(
            ) if isinstance(json[i]["_source"][key], list)})
        print(ignore)
        temp = [dict({"id": json[i]["_id"]}, **fj.flatten_json(json[i]
                     ["_source"], root_keys_to_ignore=ignore)) for i in range(len(json))]
    else:
        temp = [dict({"id": json[i]["_id"]}, **json[i]["_source"])
                for i in range(len(json))]
    temp = pd.DataFrame(temp)
    temp.index = temp.id
    return temp


def get_from_faang_api(url, post_data, size=20000, verbose=False):
    newURL = url + '/?size=' + str(size)
    response = requests.post(newURL, json=post_data)
    result = response.json()
    if verbose:
        print(f"Number of hits: {result['hits']['total']}")
    if (result['hits']['total'] > size):
        print(
            f"WARNING: more hits than expected ({size}).\nA new request will be made to include every hit.")
        return(get_from_faang_api(url, post_data, size=result['hits']['total'], verbose=verbose))
    else:
        return result['hits']['hits']


def match_exps_to_files(experiments, files, verbose=False):
    
    exps = {d["_id"]: {} for d in experiments}
    noExp = []

    ftp = FTP("ftp.sra.ebi.ac.uk")
    ftp.login()
    
    t0 = time.time()
    
    for i in range(len(files)):
        f = files[i]
        x = f["_source"]["experiment"]["accession"]
        r = f["_source"]["run"]["accession"]
        if x in exps.keys():
            exps[x].setdefault(r, [])
            furl = f["_source"]["url"]
            try:
                furl = 'ftp://' + furl[furl.index('ftp.sra.ebi'):]
            except ValueError:
                if verbose: print(f"Unusual file url: {furl}")
            
            dirs = furl.split("/")
            dir = "/".join(dirs[3:-1])
            fileName = dirs[-1]
            try:
                ftp.cwd(dir)
            except Exception as e:
                print(f"FTP directory not found: {dir}\n{e}")
                exps[x].pop(r, None)
            
            
            if fileName not in ftp.nlst():
                print(f"File not found: {fileName}")
                exps[x].pop(r, None)
            else:
                exps[x][r].append(furl)
            
            ftp.cwd("/")

            """
            cmd = "wget -q --spider '{}'".format(furl)
            dlExitStatus=os.WEXITSTATUS(os.system(cmd))
            if (dlExitStatus!=0 and dlExitStatus!=8):
                print(cmd)
                print(dlExitStatus)
                #print("Invalid link: {}. Run {} has been removed.".format(furl, r))
                name = f["_id"].split("_")
                exps[x].pop(r, None)
                noExp.append(x)
            else:
                exps[x][r].append(furl)
            """
        else:
            if verbose: print(f"WARNING: {x} not found in RNA-Seq experiments !")
            noExp.append(x)
        
        if verbose and (i % 10 == 0 and i != 0):
            dt = time.time() - t0
            est = dt/i*len(files)
            print("{}/{} - Estimated time remaining: {}".format(i, len(files), time.strftime("%H:%M:%S", time.gmtime(est))), end='\r')
            
        #if verbose: print(f"{i}/{len(files)}")
    noExp = list(set(noExp))
    if verbose and noExp:
        print("WARNING: The following experiments were not included in the list of experiments from FAANG :")
        print("\n\t".join(noExp))
    
    exps = {k:
            {
                r: exps[k][r] for r in exps[k].keys()
                if (len(exps[k][r])==2 and
                    list(map(check_paired_reads, exps[k][r]))==[True, True])
                or (len(exps[k][r])==1 and
                    list(map(check_paired_reads, exps[k][r]))==[False])
            }
            for k in exps.keys()}
    exps = {k: exps[k] for k in exps.keys() if len(exps[k]) != 0 and k not in noExp}
    return(exps)


# MAIN
if __name__ == "__main__":
    species = sys.argv[1]

    if (len(sys.argv) == 3):
        ids = sys.argv[2]
        
        try:
            with open(ids, "r") as f:
                accList = [l.strip() for l in f.readlines()]
        except FileNotFoundError:
            print("File not found ! ({})".format(ids) )
            sys.exit(1)

        shouldList = [{"match": {"accession": acc}} for acc in accList]

        post_experiment = {
            "query": {
                "bool": {
                    "should": shouldList,
                    "minimum_should_match":1,
                }
            }
        }
        
        shouldList = [{"match": {"experiment.accession": acc}} for acc in accList]
        
        post_files = {
            "query": {
                "bool": {
                    "should": shouldList,
                    "minimum_should_match":1,
                }
            }
        }
        
    else:
        post_experiment = {
            "query": {
                "bool": {
                    "must": [
                        {"match": {"standardMet": "FAANG"}},
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

        post_specimens = {
        "query": {
            "bool": {
                "must": [
                    {"match": {"organism.organism.text": species}},
                    {"match": {"standardMet": "FAANG"}}
                ]
            }
        }
    }

    # GET DATA
    url = 'https://data.faang.org/api/experiment/_search'
    experiments = get_from_faang_api(url, post_data=post_experiment, verbose=True)
    input_dna = get_from_faang_api(url, post_data=post_input, verbose=True)

    files = get_from_faang_api('https://data.faang.org/api/file/_search', post_data=post_files, verbose=True)
    specimens = get_from_faang_api('https://data.faang.org/api/specimen/_search', post_data=post_specimens, verbose=True)

    # MATCH INPUT AND DATA

    inputs = {d["_id"]: "" for d in input_dna}
    exps = {d["_id"]: d["_source"]["ChIP-seq DNA-binding"]["controlExperiment"] for d in experiments}
    all_chip = inputs
    all_chip.update(exps)

    x = experiments.copy()
    x.extend(input_dna)
    exps = match_exps_to_files(x, files, verbose=True)

    paired = []
    single = []
    #done = []

    for f in files:
        fname = f["_id"]
        #print(fname)
        furl = "ftp://" + f["_source"]["url"]
        fexp = f["_source"]["experiment"]["accession"]
        if fname[-2:] != "_2" and fexp in exps.keys():
            #done.append(fexp)
            design = {"group" : fexp, "replicate": "1"}
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
        df=df.sort_values(["control", "antibody",  "group", "replicate"])

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
        df=df.sort_values(["control", "antibody",  "group", "replicate"])

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

    ############# GET METADATA

    filesDf = json_to_df(files, flatten=False)
    filesDf["experiment_accession"]=[exp["accession"] for exp in filesDf["experiment"]]
    experimentsDf = json_to_df(experiments, flatten=False)
    #inputDf = json_to_df(input_dna, flatten=False)
    specDf = json_to_df(specimens, flatten=False)

    #allChipDf = inputDf.append(experimentsDf)
    meta = pd.merge(filesDf, experimentsDf, how='inner', left_on="experiment_accession", right_on=experimentsDf.index, suffixes=["_file", "_exp"])
    meta = pd.merge(meta, specDf, how="left", left_on="specimen", right_on=specDf.index, suffixes=["", "_spec"])
    meta.dropna(how='all', axis=1, inplace=True)

    meta.index = meta["experiment_accession"]
    meta["run_accession"] = [run["accession"] for run in meta["run"]]
    runs = meta["run_accession"].groupby(by="experiment_accession").agg(["unique"]).applymap(list)["unique"]
    files_name = meta["id_file"].groupby(by="experiment_accession").agg(["unique"]).applymap(list).applymap(sorted)["unique"]
    files_url = meta["url"].groupby(by="experiment_accession").agg(["unique"]).applymap(list).applymap(sorted)["unique"]
    files = [[{"text": text, "url":url} for text, url in zip(texts, urls)] for texts, urls in zip(files_name, files_url)]
    files = pd.Series(files, index=runs.index, name="files")
    runs = pd.Series(runs, index=files.index, name="runs")
    
    metadata=meta.copy()
    metadata.index=metadata["experiment_accession"]
    metadata=metadata.applymap(str)
    excludedCols=[]


    for c in metadata.columns:
        #breakpoint()
        counts=metadata[c].groupby(by=metadata.index).nunique()
        if max(counts)>1:
            excludedCols.append(c)
    
    metadata=meta.drop(excludedCols, axis=1)
    metadata["runs"]=runs
    metadata["files"]=files
    
    a = pd.DataFrame([metadata.loc[index,].iloc[0,] for index in metadata.index.drop_duplicates() if type(metadata.loc[index]) == type(pd.DataFrame())])
    b = pd.DataFrame([metadata.loc[index] for index in metadata.index.drop_duplicates() if type(metadata.loc[index]) == type(pd.Series(dtype='object'))])
    m = pd.concat([a, b])
    m.dropna(how='all', axis=1, inplace=True)
    m = m.astype({"paperPublished": "bool"})

    mTrimmed=m[[c for c in KEEP_COLUMNS if c in m.columns]]
    mTrim = mTrimmed.applymap(reformat)
    
    with open(f"{species}_chipseq.tsv", "w") as f:
        f.write(mTrim.to_csv(index=False, sep="\t"))
        
    with open(f"{species}_chipseq.json", "w") as f:
        f.write(mTrim.to_json(orient="index"))
