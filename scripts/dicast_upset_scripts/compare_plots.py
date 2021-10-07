import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import pandas as pd
import numpy as np
import sys
sys.path.insert(0,'../')
import UpSetPlot.upsetplot as up
import pyranges as pr
from collections import defaultdict
import seaborn as sns
import itertools as it
import warnings
from matplotlib.axes._axes import _log as matplotlib_axes_logger
matplotlib_axes_logger.setLevel('ERROR')
import os 
import argparse
import matplotlib

def arg_parser():
    parser = argparse.ArgumentParser('DICAST plots')
    parser.add_argument('--dir', dest='dir', help='Pass directory name contains all the unified.out')
    parser.add_argument('--outputdir', dest='outputdir', help='Directory path for outputs')

    return parser

def find_all_unified(dir):
    '''
    return a dict of paths for upset plot
    e.g. 
     {'S1': defaultdict(list,
                         {'50M': ['/nfs/scratch/AS_benchmarking/S1/AS_results/50M/unified/unified.out/SplAdder_filtered.unified.out',
                           '/nfs/scratch/AS_benchmarking/S1/AS_results/50M/unified/unified.out/MAJIQ_filtered.unified.out',
                           '/nfs/scratch/AS_benchmarking/S1/AS_results/50M/unified/unified.out/EventPointer_filtered.unified.out',
                           '/nfs/scratch/AS_benchmarking/S1/AS_results/50M/unified/unified.out/IRFinder_filtered.unified.out',
                           '/nfs/scratch/AS_benchmarking/S1/AS_results/50M/unified/unified.out/SGSeq_Anno_filtered.unified.out',
                           '/nfs/scratch/AS_benchmarking/S1/AS_results/50M/unified/unified.out/Whippet_filtered.unified.out',
                           '/nfs/scratch/AS_benchmarking/S1/AS_results/50M/unified/unified.out/SGSeq_denovo_filtered.unified.out'],
                          '200M': ['/nfs/scratch/AS_benchmarking/S1/AS_results/200M/unified/unified.out/SplAdder_filtered.unified.out',
                           '/nfs/scratch/AS_benchmarking/S1/AS_results/200M/unified/unified.out/MAJIQ_filtered.unified.out',
                           '/nfs/scratch/AS_benchmarking/S1/AS_results/200M/unified/unified.out/EventPointer_filtered.unified.out',
                           '/nfs/scratch/AS_benchmarking/S1/AS_results/200M/unified/unified.out/Whippet.unified.out',
                           ......
    '''

    ####find compare files
    compare_paths = []
    for dirname,_,z in os.walk(dir):
        if any("comparison.txt" in filename for filename in z):
            for filename1 in z:
                if "comparison.txt" in filename1:
                    compare_paths.append(os.path.join(dirname,filename1)) 

    compare_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for path in compare_paths:
        #if "M" in path.split('/')[6]:

        compare_dict[path.split("/")[-1].split("_")[0]][path.split("/")[-1].split('_')[1]][path.split("/")[-1].split('_')[2]].append(path)

    return compare_dict

def compare_each_plots(compare_paths, save_path, leg=True, fs=20):
    tooltype = {"ASGAL":"RNA-seq-driven", 
    "ASPLI":"Annotation-driven",
    "EVENTPOINTER":"RNA-seq-driven",
    "IRFINDER":"Annotation-driven",
    "MAJIQ":"RNA-seq-driven",
    "SGSEQ_ANNO":"Annotation-driven",
    "SGSEQ_DENOVO":"RNA-seq-driven",
    "SGSEQ":"not specified",
    "SPLADDER":"RNA-seq-driven",
    "WHIPPET":"Annotation-driven"}

    tools = []
    comparedf = []
    for file in compare_paths:
        if "comparison" in file:
            tmp = pd.read_csv(file, sep="\t")
            
            tmp['tool'] = list(tooltype.keys())[np.where([True if file.split("/")[-1].upper().count(too)==1 else False for too in tooltype.keys()])[0][0]]           
            tmp['tooltype'] = tooltype[tmp['tool'][0].upper()]
            
            comparedf.append(tmp)

    comparedf = pd.concat(comparedf)

    sns.set_style("whitegrid")
    sns.relplot(data=comparedf, x='recall', y='precision', hue="tool", style="tooltype", s=200)
    plt.xlabel("Recall",fontsize=20)
    plt.ylabel("Precision",fontsize=20)
    plt.xlim(right=1.05)
    plt.ylim(top=1.05)
    plt.savefig(os.path.join(save_path,"overall_compare.png"), dpi=500, bbox_inches="tight")
    plt.close()

    matplotlib.rcParams["legend.markerscale"] = 2.5
    markers = ['o','*','X','p','D','P',">"]
    plt.style.use('seaborn-colorblind')
    color=sns.color_palette("Set2", len(comparedf.tool.unique()))
    color=np.array(color)

    for e,typ in enumerate(comparedf.type.unique()):
        plt.figure(figsize=(8,7))

        tmp = comparedf[comparedf['type']==typ].reset_index(drop=True)
        tmp['tool']=pd.Categorical(tmp['tool'])
        g = sns.scatterplot(data=tmp, x='recall', y='precision', hue="tool", style="tooltype", s=200, legend=leg)
        if leg:
            plt.legend(bbox_to_anchor=(1.05, 1),
                borderaxespad=0, markerscale=2.5, ncol=7, title=None)
        plt.title(typ, y=1, fontsize=fs)
        plt.xlim(right=1.05)
        plt.ylim(top=1.05)
        plt.xlabel("Recall",fontsize=fs)
        plt.ylabel("Precision",fontsize=fs)
        plt.savefig(os.path.join(save_path,f"{typ}_compare.png"), bbox_inches="tight")
        plt.show()
        plt.close()

if __name__=="__main__":
    args=arg_parser().parse_args()

    if not os.path.exists(args.outputdir):
        os.mkdir(args.outputdir)

    compare_dict = find_all_unified(args.dir)
    
    for sample, v in compare_dict.items():
        for depth, vv in v.items():
            for mapper, paths in vv.items():
            
                print(f"Preparing compare plot for {sample}, {depth}, {mapper}")

                if not os.path.exists(os.path.join(args.outputdir,sample)):
                    os.mkdir(os.path.join(args.outputdir,sample))
                if not os.path.exists(os.path.join(args.outputdir,sample,depth)):
                    os.mkdir(os.path.join(args.outputdir,sample,depth))
                if not os.path.exists(os.path.join(args.outputdir,sample,depth, mapper)):
                    os.mkdir(os.path.join(args.outputdir,sample,depth, mapper))

                save_path = os.path.join(args.outputdir, sample, depth, mapper)
                compare_each_plots(paths, save_path)

                print("done")