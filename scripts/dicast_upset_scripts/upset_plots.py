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
matplotlib.use('Agg')
import re
pd.options.mode.chained_assignment = None
from pandas.errors import EmptyDataError
plt.ioff()

def arg_parser():
    parser = argparse.ArgumentParser('DICAST plots')
    parser.add_argument('--dir', dest='dir', help='Pass directory name contains all the unified.out')
    parser.add_argument('--outputdir', dest='outputdir', help='Directory path for outputs (exists or not also works)')

    return parser

def find_all_unified(dir):
    '''
    return a dict of paths for upset plot
    e.g. 
     {'S1': defaultdict(list,
                         {'50M': ['/nfs.../.../.../../S1/AS_results/50M/unified/unified.out/SplAdder_filtered.unified.out',
                           '/nfs.../.../.../../S1/AS_results/50M/unified/unified.out/MAJIQ_filtered.unified.out',
                           '/nfs.../.../.../../S1/AS_results/50M/unified/unified.out/EventPointer_filtered.unified.out',
                           '/nfs.../.../.../../S1/AS_results/50M/unified/unified.out/IRFinder_filtered.unified.out',
                           '/nfs.../.../.../../S1/AS_results/50M/unified/unified.out/SGSeq_Anno_filtered.unified.out',
                           '/nfs.../.../.../../S1/AS_results/50M/unified/unified.out/Whippet_filtered.unified.out',
                           '/nfs.../.../.../../S1/AS_results/50M/unified/unified.out/SGSeq_denovo_filtered.unified.out'],
                          '200M': ['/nfs.../.../.../../S1/AS_results/200M/unified/unified.out/SplAdder_filtered.unified.out',
                           '/nfs.../.../.../../S1/AS_results/200M/unified/unified.out/MAJIQ_filtered.unified.out',
                           '/nfs.../.../.../../S1/AS_results/200M/unified/unified.out/EventPointer_filtered.unified.out',
                           '/nfs.../.../.../../S1/AS_results/200M/unified/unified.out/Whippet.unified.out',
                           ......
    '''
    unified_paths = []
    for dirname,_,z in os.walk(dir):
        if any("unified.out" in filename for filename in z):
            for filename in z:
                if "unified.out" in filename:
                    unified_paths.append(os.path.join(dirname,filename)) 

    unified_dict = defaultdict(list)
    for path in unified_paths:
        if "unmapped" in path:
            continue
        sample = "_".join(path.split("/")[-1].split("_")[:-5])
        #maptool = path.split("/")[-1].split("_")[3]
        # if len(maptool)==0:
        #     maptool = path.split("/")[-1].split("_")[4]
        unified_dict[sample].append(path)

    return unified_dict

### for MEE, MES
def expand_coord(df1):
    dropped = df1.drop(['Start', 'End'],axis=1)
    dropped['index']=dropped.index
    exploded = pd.concat([df1['Start'].str.split(",").explode(),df1['End'].str.split(",").explode()], axis=1)
    exploded['index'] = exploded.index

    merged1 = pd.merge(dropped, exploded, how="inner", on='index')
    
    d1 = defaultdict(list)
    for ind in merged1['index'].unique():
        tmp = merged1[merged1['index']==ind]
        for i, r in tmp.iterrows():
            d1[ind].append(r)

    return d1

def create_row(x, tools):
    countdict = {}
    for t in tools:
        countdict.setdefault(t, 0)

    for xx in x:
        countdict[xx]+=1
    return list(countdict.values())    

def find_overlaps(event1, event2, combi):
    complete_event1 = pd.DataFrame(event1)
    complete_event2 = pd.DataFrame(event2)
    thisdict = {combi[0]:complete_event1.to_dict(), combi[1]:complete_event2.to_dict()}

    try:
        grs = {n: pr.from_dict(s).drop_duplicate_positions(keep=False) for n, s in thisdict.items()}
    except ValueError:
        return False

    counts = pr.count_overlaps(grs)
    countdf = counts.df

    MATCHED = False
    for ind, row in counts.df[[combi[0], combi[1]]].iterrows():
        if row.sum() == 2:
            MATCHED = True

    return MATCHED


def find_all_overlaps(mergin, merged1, merged2, combi):
    #compare tool0 event to all other event from other tool
    mergin0_res = [find_overlaps(merged1[mergin], merged2[x], combi) for x in merged2.keys()]
    #mergin1_res = [find_overlaps(merged1[x], merged2[mergin[1]]) for x in merged1.keys()]       

    return mergin0_res#, mergin1_res

def make_plots(target_paths):
    unified = defaultdict(lambda: defaultdict(list))
    tools = ["asgal", "aspli", "eventpointer", "irfinder", "majiq", "sgseq", "spladder", "whippet"]
    evs = []
    #def unified_upset():
    for file in target_paths:
        if "unified.out" in file:
            # if len(file.split("/")[-1].split(".")) > 2:
            #     tool = file.split("/")[-1].split(".")[1].split("_")[-2]
            # else:
            #     tool = file.split("/")[-1].split(".")[0].split("_")[-2]
            tool = tools[np.where(list(map(lambda x: x in file, tools)))[0][0]]
            if "_filtered" in tool:
                tool = tool.replace("_filtered","")

            try:
                tmp = pd.read_csv(file, sep="\t")
            except EmptyDataError as e:
                print(f"This file returns the following error: {file}")
                print(e)
                continue

            tmp = tmp.dropna()
            for ev in tmp.event_type.unique():
                events = tmp[tmp['event_type']==ev]
                org=events['chr'].copy(deep=False)

                events.loc[:,['chr']]=list(map(lambda x: "chr"+str(x), org))
                events.columns = ["Chromosome","gene","id","strand","event_type","count","Start","End"]
                unified[ev][tool].append(events.to_dict('list'))
                if ev not in set(evs):
                    evs.append(ev)
    
    allcomb = dict()
    for ev, X in unified.items():

        if ev=='MEE' or ev == 'MES':
            realcount = pd.DataFrame(columns=["  font-size: 1rem;Chromosome","Start","End"]+tools)
            for combi in it.combinations(X.keys(),2):
                if len(combi)<2:
                    df1 = pd.DataFrame(X[combi[0]][0]).reset_index()
                    row1 = create_row([combi[0]])
                    for _ in range(df1.shape[0]):
                        realcount.loc[realcount.shape[0]] = ["chr",0,0]+row1
                    continue

                
                df1 = pd.DataFrame(X[combi[0]][0]).reset_index()
                df1['index']=df1.index
                df2 = pd.DataFrame(X[combi[1]][0]).reset_index()
                df2['index']=df2.index

                
                merged1 = expand_coord(df1)
                merged2 = expand_coord(df2)
                
                matched_index = []
                for mergin in merged1.keys():
                    mergin0_res = find_all_overlaps(mergin, merged1, merged2, combi) 

                    if any(mergin0_res):
                        row1 = create_row([combi[0], combi[1]], tools)
                        matched_index.append(np.where(mergin0_res)[0][0]) #keep track which index of event is found overlapped, so that doesn't duplicate
                    #realcount.loc[realcount.shape[0]] = ["chr", 0, 0] +
                    else:
                        row1 = create_row([combi[0]], tools)
                    
                    
                    realcount.loc[realcount.shape[0]] = ["chr",0,0]+row1
                
                #add the events in tool 2 that doesn't have overlaps
                row2 = create_row([combi[1]], tools)
                for _ in range(len(merged2)-len(matched_index)):
                    realcount.loc[realcount.shape[0]] = ["chr",0,0]+row2
        
        else:
            grs = {n: pr.from_dict(s[0]).drop_duplicate_positions(keep=False) for n, s in X.items()}
            counts = pr.count_overlaps(grs)
            countdf = counts.df

            #check if there are tools left out
            missed = [tools[x] for x in np.where(np.isin(tools, countdf.columns)==False)[0]]
            if len(missed) > 0:
                for x in missed:
                    countdf[x] = 0
            realcount = countdf[["Chromosome","Start","End"]+tools]
        
        
        for row in realcount.itertuples():
            tmp = list(row[4:])
            tmp = [1 if x>1 else x for x in tmp]
            binkey=''.join([str(x) for x in tmp])
            if np.sum(tmp)==0:
                continue
            else:
                if binkey not in set(allcomb.keys()):
                    allcomb.setdefault(binkey, {})
                    for e in evs:
                        allcomb[binkey].setdefault(e, 0)
                

                allcomb[binkey][ev] += 1


    forplot = pd.DataFrame(columns=tools+evs)
    
    for n,j in allcomb.items():
        thisrow= [bool(int(x)) for x in n] + list(j.values())
        forplot.loc[forplot.shape[0]] = thisrow

    forplot = forplot.set_index(tools)
    return forplot



if __name__=="__main__":
    args=arg_parser().parse_args()
    print('Intersections of events identified are being plotted to: ',(args.outputdir))
    if not os.path.exists(args.outputdir):
        os.mkdir(args.outputdir)

    unified_dict = find_all_unified(args.dir)

    for sample, paths in unified_dict.items():                
      if not os.path.exists(os.path.join(args.outputdir,sample)):
          os.mkdir(os.path.join(args.outputdir,sample))



      print(f"Preparing upset plot for {sample}")
      if len(paths) < 2:
          print("Not enough unified output. Skipping this.")
          continue
      plot_this_df = make_plots(paths)
      up.plot(plot_this_df, sort_by="cardinality")
      plt.savefig(os.path.join(args.outputdir,sample,"upset_plot.png"), bbox_inches="tight")
      plt.close()
      print('Upset plots done')


    
            



    
