from .utils import *

def vep_cleaning(df):
    #pd.options.display.float_format = '{:,.6f}'.format
    df['GERP'] = df['GERP'].apply(get_max_if_multiple)
    df['phyloP'] = df['phyloP'].apply(get_max_if_multiple)
    df['phastCons'] = df['phastCons'].apply(get_max_if_multiple)
    df['SiPhy'] = pd.to_numeric(df['SiPhy'])
    df['fitCons'] = pd.to_numeric(df['fitcons'])
    df['LINSIGHT'] = pd.to_numeric(df['LINSIGHT'])
    df['MutationAssessor'] = df['MutationAssessor_score'].apply(get_max_if_multiple)
    df['LRT'] = simple_merge(df['LRT_score'],df['LRT_pred'])
    df['Sift'] = df['SIFT_score'].apply(get_min_if_multiple)
    df['Polyphen2HVAR'] = df['Polyphen2_HVAR_score'].apply(get_max_if_multiple)
    df['Polyphen2HDIV'] = df['Polyphen2_HDIV_score'].apply(get_max_if_multiple)
    df['MutationTaster'] = df['MutationTaster_score'].apply(get_max_if_multiple)
    df['FATHMM-MKL'] = pd.to_numeric(df['FATHMM-MKL'])
    df['FATHMM'] = df['FATHMM_score'].apply(get_min_if_multiple)
    df['Provean'] = df['PROVEAN_score'].apply(get_min_if_multiple)
    df['Mutpred'] = pd.to_numeric(df['MutPred_score'])
    df['GWAVA'] = df['GWAVA'].apply(get_max_if_multiple)
    df['VEST4'] = df['VEST4_score'].apply(get_max_if_multiple)
    df['ExACpLI'] = pd.to_numeric(df['ExACpLI'])
    df['CADD'] = pd.to_numeric(df['CADD_PHRED'])
    df['DANN'] = df['DANN'].apply(process_spidex_remm_dann)
    df['MetaLR'] = simple_merge(df['MetaLR_score'],df['MetaLR_pred'])
    df['MetaSVM'] = simple_merge(df['MetaSVM_score'],df['MetaSVM_pred'])
    df['Revel'] = pd.to_numeric(df['REVEL_score'])
    df['Eigen'] = pd.to_numeric(df['Eigen'])
    # df['Eigen-PC'] = pd.to_numeric(df['Eigen-PC'])
    df['ReMM'] = df['ReMM'].apply(process_spidex_remm_dann)
    df['CAROL'] = condel_carol(df['CAROL'])
    df['Condel'] = condel_carol(df['Condel'])
    df['M-CAP'] = pd.to_numeric(df['M-CAP_score'])
    df['traP'] = df['traP'].apply(get_max_if_multiple)
    df['SPIDEX'] = df['dpsi_zscore'].apply(process_spidex_remm_dann)
    df['dbscSNV'] = dbscSNV_merge(df['ada_score'], df['rf_score'])
    df['MaxEntScan'] = maxentscan_merge(df['MaxEntScan_alt'], df['MaxEntScan_ref'], df['MaxEntScan_diff'])
    df.apply(get_spliceai_loc, axis=1)
    df['SpliceAI'] = df['SpliceAI'].apply(process_spliceai)
    return df