import numpy as np
import pandas as pd
import os
import fnmatch

def tuple2float(x):
    return x[1] if x[1] != '.' or x[1] != '' or x[1] else np.nan

def check_required_files(dir,analysis):
    fname_benign = [filename for filename in os.listdir(dir) if fnmatch.fnmatch(filename,"*benign*gz")]
    fname_deleterious = [filename for filename in os.listdir(dir) if fnmatch.fnmatch(filename,"*deleterious*gz") or fnmatch.fnmatch(filename,"*pathogenic*gz")]
    if fname_benign and fname_deleterious:
        benign = os.path.join(dir, fname_benign[0])
        deleterious = os.path.join(dir,fname_deleterious[0])
        return benign,deleterious

    print("The required files for {} analysis are missing.".format(analysis))
    exit(1)


def get_min_if_multiple(v):
    if type(v) == str:
        if v == 'None' or not v or all(x == '.' for x in v.split("&")):
            return np.nan
        try:
            return min([float(x) for x in v.split("|") if x.strip() != 'None'])
        except ValueError:
            if len(v.split("&")) > 1:
                return min([float(x) for x in v.split("&") if x.strip() != '.' ]) #if from VEP
            else:
                return min([float(x) for x in v.split(',') if x.strip() != '.'])
    else:
        return v

def get_max_if_multiple(v):
    if type(v) == str:
        if v == 'None' or not v or all(x == '.' for x in v.split("&")):
            return np.nan
        try:
            return max([float(x) for x in v.split("|") if x.strip() != 'None'])
        except ValueError:
            if len(v.split("&")) > 1:
                return max([float(x) for x in v.split("&") if x.strip() != '.']) #if from VEP
            else:
                return max([float(x) for x in v.split(',') if x.strip() != '.']) #if from vcfanno
    else:
        return v

def custom_range(series, steps):
    mi = series.min()
    ma = series.max()
    step = (ma-mi) / float(steps)
    return np.arange(mi, ma, step)
    
def simple_merge(scores, preds):
    final = []
    for pred, score in zip(preds, scores):
        if pred == np.nan or not pred or score == np.nan or not score:
            v = np.nan
  #      elif pred == 'N':
  #          v = -score
        else:
            v = float(score)
        final.append(v)
    return pd.to_numeric(final)

def polyphen_merge(scores):
    final = []
    for score in scores:
        if score == np.nan or not score or score == ".":
            v = np.nan
        elif type(score) == float:
            v = score
        else:
            highest_score = -100
            sep="|" if "|" in score else "&"
            for s in score.split(sep):
                if s == 'None' or s == ".":
                    continue
                s = float(s)
                if s > highest_score:
                    highest_score = s
            v=highest_score
        final.append(v)
    return pd.to_numeric(final)


def fathmm_merge(scores,preds):
    final = []
    for pred, score in zip(preds, scores):
        if pred == np.nan or not pred or score == np.nan or not score or pred == "." or score == ".":
            v = np.nan
        elif type(score) == float:
            v = score
        else:
            lowest_score = 100
            sep = "|" if "|" in score else "&"
            for p, s in zip(pred.split(sep), score.split(sep)):
                try:
                    s=float(s)
                except ValueError:
                    continue
                if s < lowest_score:
                    lowest_score = s
            v = lowest_score
        final.append(v)

    return pd.to_numeric(final)

def condel_carol(scores):
    final=[]
    for score in scores:
        if not score:
            s=np.nan
        elif score.lower().startswith('d') or score.lower().startswith('n'):
            s=float(score.split("(")[1].split(")")[0])
        final.append(s)
    return pd.to_numeric(final)

def dbscSNV_merge(ada_score,rf_score):
    pairs=zip(ada_score,rf_score)
    final=[]
    for (a,r) in pairs:
        try:
            if float(a) >= 0.6 and float(r) >= 0.6:
                final.append(max(float(a),float(r)))
            elif any(float(i) >= 0.6 for i in (a,r)):
                final.append(min(float(a),float(r)))
            elif isinstance(a, float) and isinstance(r,float): #if nan
                final.append(a)
            else:
                final.append(max(float(a),float(r))) #if both smaller than 0.6
        except ValueError:
            final.append(np.nan)
    return pd.to_numeric(final)

def maxentscan_merge(alt,ref,diff):
    final=[]
    for a, b, c in zip(alt, ref, diff):
        if isinstance(c,float):
            final.append(c)
        elif c:
            final.append(abs(float(c)))
        else:
            final.append(np.nan)
    return pd.to_numeric(final)

def spliceRegion_merge(spliceregion):
    final=[]
    [final.append(1) if isinstance(v,str) and v.startswith("splice") else final.append(0) for v in spliceregion]
    return pd.to_numeric(final)


def genesplicer_test(genesplicer): #Emailed author
    final=[]
    for v in genesplicer:
        if isinstance(v,float):
            final.append(v)
        elif v:
            print(v)
    #print(len(final))

def format_spliceai_fields(record):

    gname = record.INFO.get("gene_name")

    fields = [record.ALT[0], gname, record.INFO.get("DS_AG"), record.INFO.get("DS_AL"),
     record.INFO.get("DS_DG"), record.INFO.get("DS_DL"),
     record.INFO.get("DP_AG"), record.INFO.get("DP_AL"),
     record.INFO.get("DP_DG"), record.INFO.get("DP_DL")]

    if gname and ',' in gname: #if multiple predictions on different genes
        gene_from_id = record.ID.split(":")[0]
        try:
            idx = gname.split(",").index(gene_from_id)
            fields = [record.ALT[0], gname.split(",")[idx], record.INFO.get("DS_AG")[idx], record.INFO.get("DS_AL")[idx],
                      record.INFO.get("DS_DG")[idx], record.INFO.get("DS_DL")[idx],
                      record.INFO.get("DP_AG")[idx], record.INFO.get("DP_AL")[idx],
                      record.INFO.get("DP_DG")[idx], record.INFO.get("DP_DL")[idx]]

        except ValueError:
            idx = 0
            fields = [record.ALT[0], gname.split(",")[idx], record.INFO.get("DS_AG")[idx], record.INFO.get("DS_AL")[idx],
                      record.INFO.get("DS_DG")[idx], record.INFO.get("DS_DL")[idx],
                      record.INFO.get("DP_AG")[idx], record.INFO.get("DP_AL")[idx],
                      record.INFO.get("DP_DG")[idx], record.INFO.get("DP_DL")[idx]]

    return '|'.join(map(str, [0 if x is None else x for x in fields]))


def get_spliceai_loc(x):
    if x.SpliceAI.split("|")[1] == "0" and x.location == "unknown" or x.location == "mithocondria":
        x.SpliceAI=None

def process_spliceai(v):
    if v is None:
        return np.nan
    else:
        return max(map(float, v.split("|")[2:6]))

def process_spidex_remm_dann(x):

    if isinstance(x,float):
        return abs(x)
    elif x == "." or not x:
        return np.nan
    else:
        return abs(float(x))