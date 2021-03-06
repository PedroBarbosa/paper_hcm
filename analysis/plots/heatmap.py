import os.path

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
plt.switch_backend('agg')
sns.set(style="white")

GT_LABEL = '(*) Ground Truth'

def plot_heatmap_unlabelled(df):
    f, ax = plt.subplots(figsize=(0.5 * len(df.columns), 6))
    cmap = sns.diverging_palette(10, 220, as_cmap=True)
    sns.heatmap(df, cmap=cmap,  linewidths=0, cbar=False)
    plt.show()
def prepare_dataset_for_heatmap(df, thresholds):
    df['blank'] = pd.Series(np.nan, index = np.arange(df.shape[0]))
    df = df[['class'] + [ tool + "_prediction" for tool, *args in thresholds]].copy()
    df.columns = [GT_LABEL] + [ tool for tool, *args in thresholds]
    f = lambda x: pd.isna(x) and np.nan or (1-int(x)) 
    for col in df.columns:
        df[col] = df[col].apply(f)
        if df[col].isnull().all():
            del df[col]
    df = df.reset_index(drop=True)
    df = df.sort_values(by=GT_LABEL)
    return df

def hide_ylabel():
    fr = plt.gca()
    fr.axes.yaxis.set_ticklabels([])

def plot_heatmap(df, filtername, thresholds, fname):
    df = prepare_dataset_for_heatmap(df, thresholds)
    plt.clf()
    f, ax = plt.subplots(1,2,figsize=(0.5 * len(df.columns), 6))
    cmap="PiYG"
    cmap = sns.diverging_palette(10, 220, as_cmap=True)
    sns.heatmap(df[[GT_LABEL]], ax=ax[0], cmap=cmap, linewidths=0, cbar=False)
    hide_ylabel()
    sns.heatmap(df.drop(GT_LABEL,axis=1), ax=ax[1], cmap=cmap,  linewidths=0, cbar=False, mask=df.drop(GT_LABEL,axis=1).isnull())
    hide_ylabel()
    plt.savefig(fname, bbox_inches='tight')

def generate_heatmap(df_original, filters, thresholds, name, folder):
    for filtername, filterfunction in filters:
        df = filterfunction(df_original).copy()
        if df.shape[0] > 10:
            plot_heatmap(df,filtername, thresholds, os.path.join(folder, "figures", "heatmap_{}_{}.pdf".format(name, filtername)))