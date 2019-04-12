import os.path
from .utils import *
from .vep import *
from .vep_cleaning import vep_cleaning
DATASET_FOLDER = "../datasets/"

def preprocess(loc,analysis):
    dfs=[]
    benign, deleterious = check_required_files(os.path.join(DATASET_FOLDER, analysis), analysis)
    for f in [benign,deleterious]:
        isbenign= False if "pathogenic" in f or "deleterious" in f else True
        dfs.append(get_df_ready(f, False, isbenign, loc))

    df = pd.concat(dfs)
    df = vep_cleaning(df)
    return df
