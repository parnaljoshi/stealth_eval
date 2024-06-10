import sys
sys.path.append("/home/rashika/CAFA4/InformationAccretion/")
from ia import *
import os


def delete_extra_go_terms(terms_df, common_terms):
    for aspect in ['BPO', 'CCO', 'MFO']:
        terms_to_delete = set(terms_df.loc[terms_df["aspect"]==aspect, "term"]).difference(common_terms[aspect])
        # Use boolean indexing to filter rows
        terms_df = terms_df[~terms_df['term'].isin(terms_to_delete)].copy()
        terms_df = terms_df.reset_index(drop=True)
    return terms_df

def keep_common_go_terms(t0, t1, t0_ont, t1_ont):
    
    # Get the three subontologies
    roots = {'BPO': 'GO:0008150', 'CCO': 'GO:0005575', 'MFO': 'GO:0003674'}
    t0_subont = {aspect: fetch_aspect(t0_ont, roots[aspect]) for aspect in roots}
    t1_subont = {aspect: fetch_aspect(t1_ont, roots[aspect]) for aspect in roots}
    
    common_terms = {aspect: set(t1_subont[aspect].nodes).intersection(set(t0_subont[aspect].nodes)) for aspect in ['BPO', 'CCO', 'MFO']}
    
    t0_eval = delete_extra_go_terms(t0, common_terms)
    t1_eval = delete_extra_go_terms(t1, common_terms)
    return t0_eval, t1_eval

def get_variable_name(value, namespace):
    for name, var in namespace.items():
        if var is value:
            return name
    return None

def get_annot_list(annot_df):
    """
    Converts the DataFrame of annotations to a nested dictionary, where
    the protein IDs are keys for 1st level, and the names MFO, BPO, CCO are the keys for 2nd level
    
    annot_list: Annotations in a nested dict format
    """
    proteins = np.unique(annot_df['EntryID'])
    annot_list = {}
    for i in range(len(proteins)):
        annot_list[proteins[i]] = {}
        
        for ont in ['BPO', 'CCO', 'MFO']:
            annot_list[proteins[i]][ont] = []
        
    for i in range(len(annot_df)):
        annot_list[annot_df.loc[i, 'EntryID']][annot_df.loc[i, "aspect"]].append(annot_df.loc[i, "term"])
        
    return annot_list
        
    
def get_annot_gained(t0_annot_list, t1_annot_list):
    """
    Returns the annotations gained by the proteins in each ontology
    
    annot_list: Annotations in a nested dict format where,
    the protein IDs are keys for 1st level, and the names MFO, BPO, CCO are the keys for 2nd level
    """
    proteins = list(set(t1_annot_list.keys()))
    
    annot_gained = {}
    for protein in proteins:
        annot_gained[protein] = {}
        for ont in ['BPO', 'CCO', 'MFO']:
            
            # If the protein is in the t0 annotations, the set difference of annotations at t1 and t0, are the annotations gained
            if protein in t0_annot_list.keys():
                annot_gained[protein][ont] = list(set(t1_annot_list[protein][ont]).difference(set(t0_annot_list[protein][ont])))
            
            # If the protein is not in the t0 annotations, all its t1 annotations are the annotations gained
            else:
                annot_gained[protein][ont] =  t1_annot_list[protein][ont]
    return annot_gained


def get_baselines(t0_annot_list, t1_annot_list, BM_path):
    """
    Returns proteins that make it to the benchmarks
    
    type1 (No Knowledge) No annotation in any ontology at t0, annotation in the evaluation ontology at t1
    type2 (Limitied Knowledge) Annotation in other ontologies but not in evaluation ontology at t0, annotations in the evaluation ontology at t1
    type3 (Partial Knowledge) Annotations in the evaluation ontology at t0, more annotations in the evaluation ontology at t0 (may or may not have annotations in other ontologies at t0 and t1
    """
    annot_gained = get_annot_gained(t0_annot_list, t1_annot_list)
    
    proteins = list(annot_gained.keys())
    
    type1 = {}
    type2 = {}
    type3 = {}
    type12 = {}
    type1_annot_gained = {}
    type2_annot_gained = {}
    type3_annot_gained = {}
    
    for ont in ['BPO', 'CCO', 'MFO']:
        type1[ont] = []
        type2[ont] = []
        type3[ont] = []
        type12[ont] = []
    
    for prot in proteins:
        
        # type 1 - No knowledge
        # If a protein does not have any annotation at t0, but has annotation in evaluation ontology at t1,
        # then it gets appended to type1 at that evaluation ontology
        if prot not in t0_annot_list.keys():
            type1_annot_gained[prot] = {}
            for ont in ['BPO', 'CCO', 'MFO']:
                if annot_gained[prot][ont]:
                    type1[ont].append(prot)
                    type1_annot_gained[prot][ont] = annot_gained[prot][ont]
                    
        # type 2 - Partial Knowledge
        # If a protein has annotations at t0, but not in the evaluation ontology,
        # and it gains annotations in the evaluation ontology
        # then it gets appended to type2 at that ontology            
        if prot in t0_annot_list.keys():
            type2_annot_gained[prot] = {}
            for ont in ['BPO', 'CCO', 'MFO']:
                if not t0_annot_list[prot][ont] and t1_annot_list[prot][ont]:
                    type2[ont].append(prot)
                    type2_annot_gained[prot][ont] = annot_gained[prot][ont]
        
        # type 3 - If a protein has annotations at t0, in the evaluation ontology,
        # and it gains annotations in the evaluation ontology
        # then it gets appended to type3 at that ontology         
        if prot in t0_annot_list.keys():
            type3_annot_gained[prot] = {}
            for ont in ['BPO', 'CCO', 'MFO']:
                if t0_annot_list[prot][ont] and annot_gained[prot][ont]:
                    type3[ont].append(prot)
                    type3_annot_gained[prot][ont] = annot_gained[prot][ont]
        
        # type 12 - union of type 1 and type 2
        #for ont in ['BPO', 'CCO', 'MFO']:
        #    type12[ont] = list(set(type1[ont] + type2[ont]))
    
    type1_df = flatten_nested_dict(type1_annot_gained)
    type2_df = flatten_nested_dict(type2_annot_gained)
    type3_df = flatten_nested_dict(type3_annot_gained)
    type12_df = pd.concat([type1_df, type2_df], axis=0)
    #type3_df.to_csv('/home/rashika/CAFA4/Type3_ann.csv', index=False, sep = "\t", header = False) 
    
    # Check if the BM_path exists, create the directory if it does not
    if not os.path.exists(BM_path):
        os.mkdir(BM_path)
    
    for type_df in [type1_df, type2_df, type3_df, type12_df]: 
        for ont in ["BPO", "CCO", "MFO"]:
            bl_file_path = BM_path + ont.lower() + "_all" + "_"+ get_variable_name(type_df, locals()).split("_")[0] +".txt"
            ont_df = type_df[type_df["aspect"] == ont].copy()
            #aspect_mapping = {'BPO': 'P', 'MFO': 'F', 'CCO':'C'}
            #ont_df['aspect'] = ont_df['aspect'].map(aspect_mapping)
            print(bl_file_path)
            ont_df.to_csv(bl_file_path, index=False, sep = "\t", header = False)  
    
    
    #return type1_df, type2_df, type3_df, type12_df
    #return type1, type2, type3, type12

def flatten_nested_dict(nested_dict):
    flat_dict = [[outer_key, inner_key, value] for outer_key, inner_dict in nested_dict.items() for inner_key, values in inner_dict.items() for value in values]
    df = pd.DataFrame(flat_dict) 
    df.columns = ['EntryID', 'aspect', 'term']
    df = df[['EntryID', 'term', 'aspect']].copy()
    return df



def process_raw_annot(ann_file_name, ont_graph, roots, remove_roots = True):
    """
    Propagates the ann_file in ann_ont and removes the roots
    """
    ann = pd.read_csv(ann_file_name,  sep = '\t', header = None)
    if(len(ann.columns)>3): # Fix the initial ann preprocessing so that this is always 3
        ann = ann.iloc[:,:3].copy()
    ann.columns = ['EntryID', 'term', 'aspect']
    
    # Map aspect
    aspect_mapping = {
    'C': 'CCO',
    'F': 'MFO',
    'P': 'BPO'}
    
    ann['aspect'] = ann['aspect'].map(aspect_mapping)
    
    # Propagate
    subontologies = {aspect: fetch_aspect(ont_graph, roots[aspect]) for aspect in roots}
    ann_prop = propagate_terms(ann, subontologies)
    
    
    if remove_roots==True:
        # Remove the roots
        ann_prop = ann_prop[~ann_prop['term'].isin(roots.values())].copy()
        
    
    return ann_prop
    
# roots is a list of the root terms
def create_bm_lists(t0_file, t1_file, t0_ont_graph, t1_ont_graph, roots, BM_path = "/data/rashika/CAFA4/eval/BM_GO/", common_path = '/data/rashika/CAFA4/common/'):
    
    #Prop t0 and t1 in their respective ontologies
    t0_prop = process_raw_annot(t0_file, t0_ont_graph, roots)
    t1_prop = process_raw_annot(t1_file, t1_ont_graph, roots)
    
    # Keep common terms
    t0_common, t1_common =  keep_common_go_terms(t0_prop, t1_prop, t0_ont_graph, t1_ont_graph)
    t0_common.to_csv(common_path + 't0.tsv', sep = '\t', header = False, index = False)
    t1_common.to_csv(common_path + 't1.tsv', sep = '\t', header = False, index = False)

    # Propagate back in the t0 ontology
    subontologies = {aspect: fetch_aspect(t0_ont_graph, roots[aspect]) for aspect in roots}
    t0_eval = propagate_terms(t0_common, subontologies)
    t1_eval = propagate_terms(t1_common, subontologies)

    # Convert the eval Dfs into annotation lists
    t0_annot_list = get_annot_list(t0_eval)
    t1_annot_list = get_annot_list(t1_eval)
    
    # Convert the eval Dfs into annotation lists
    get_baselines(t0_annot_list, t1_annot_list,  BM_path)
    

def shortlist_BM_GO(t1_truth, BM_path, BM_GO_path):
    """Shortlist the GO terms for the respective benchmark protein files, and write them"""
    dir_list = os.listdir(BM_path)
    for file in dir_list:
        print(file)
        if 'all' in file:
            x = file.split("_")[0]
            #print(x, file)
    
            # P: process-BPO, C: component-CCO,  F: function-mfo.
            tgt = pd.read_csv(BM_path + file, sep='\t', header = None)
            tgt.columns = ['EntryID']
            display(t1_truth)
            for aspect in ['bpo', 'cco', 'mfo']:
                if aspect in file:
                    tgt_GO = t1_truth.loc[(t1_truth.EntryID.isin(tgt.EntryID) & (t1_truth.aspect == aspect.upper())), :].copy()

            if 'xxo' in file:
                tgt_GO = t1_truth.loc[t1_truth.EntryID.isin(tgt.EntryID), :].copy()
            
            # Check if the BM_GO_path exists, create the directory if it does not
            if not os.path.exists(BM_GO_path):
                os.mkdir(BM_GO_path)
            
            out_file = BM_GO_path + file.split(".")[0] + '_GO.tsv'
            tgt_GO.to_csv(out_file, sep = '\t', header = False, index = False)


def calc_IA(BM_GO_path, t0_ont_file, IA_path = "/home/rashika/CAFA4/eval/IA/"):
    annot_list = os.listdir(BM_GO_path)
    
    # Check if the BM_GO_path exists, create the directory if it does not
    if not os.path.exists(IA_path):
        os.mkdir(IA_path)
        
    for file in annot_list:
        print("Calculating Information Accretion For..")
        print(file)
        file_path = BM_GO_path + file
        out_file = IA_path+'IA_'+ file.split(".")[0] + '.txt'
        cmd = 'python3 /home/rashika/CAFA4/InformationAccretion/ia.py --annot ' + file_path + ' --graph '+ t0_ont_file + ' --outfile ' + out_file + ' --prop &' 
        os.system(cmd)
        
        
def run_eval(BM_GO_path, pred_dir, ont_file, IA_file = '/data/rashika/CAFA4/eval/IA/IA.txt', result_path = '/home/rashika/CAFA4/eval/eval_results/', log_path = '/home/rashika/CAFA4/eval/log/'):
    dir_list = os.listdir(BM_GO_path) # out_path is the path to the folder containing the target GO sets

    if not os.path.exists(result_path):
        os.mkdir(result_path)
        
    if not os.path.exists(log_path):
        os.mkdir(log_path)
        
    for file in dir_list:

        print("Evaluating: " + file)
        out_dir = result_path + file.split(".")[0] + '/'
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        
        # cafaeval go-basic.obo prediction_dir test_terms.tsv -ia IA.txt -prop fill -norm cafa -th_step 0.001 -max_terms 500
        #cmd = 'cafaeval /data/yisupeng/sharing/cafa4/gene_ontology_edit.obo.2020-01-01 /data/yisupeng/sharing/cafa4/all_models/ ' + '/data/yisupeng/sharing/cafa4/t1_truth.csv' + ' -out_dir '+ out_dir + ' -prop max -th_step 0.01  -no_orphans -log_level info > '+ log_path + file.split(".")[0] + '.log'+ ' &'
        #cmd = 'cafaeval '+ ont_file + pred_dir + BL_GO_path+file +' -out_dir '+ out_dir + ' -prop max -th_step 0.01  -no_orphans -log_level info > '+ log_path + file.split(".")[0] + '.log'+ ' &'
        #With IA
        cmd = "python3 /home/rashika/CAFA4/CAFA-evaluator/src/cafaeval/__main__.py "+ ont_file +" "+ pred_dir + " " + BM_GO_path+file + " -out_dir " + out_dir + ' -ia ' + IA_file + " -prop max -th_step 0.01  -no_orphans " + " &"
        #Without IA
        #cmd = "python3 /home/rashika/CAFA4/CAFA-evaluator/src/cafaeval/__main__.py "+ ont_file +" "+ pred_dir + " " + BM_GO_path+file + " -out_dir " + out_dir + " -prop max -th_step 0.01  -no_orphans " + " &"
        
        
        print(cmd)
        os.system(cmd)


def create_plots(results_path, metric, cols,out_path='/home/rashika/CAFA4/eval/plots/', n_curves = None, names_file = None):
    dir_list = os.listdir(results_path)
    
    cumulate = True
    add_extreme_points = True
    coverage_threshold = 0.3
    axis_title_dict = {'pr': 'Precision', 'rc': 'Recall', 'f': 'F-score', 'pr_w': 'Weighted Precision', 'rc_w': 'Weighted Recall', 'f_w': 'Weighted F-score', 'mi': 'Misinformation (Unweighted)', 'ru': 'Remaining Uncertainty (Unweighted)', 'mi_w': 'Misinformation', 'ru_w': 'Remaining Uncertainty', 's': 'S-score', 'pr_micro': 'Precision (Micro)', 'rc_micro': 'Recall (Micro)', 'f_micro': 'F-score (Micro)', 'pr_micro_w': 'Weighted Precision (Micro)', 'rc_micro_w': 'Weighted Recall (Micro)', 'f_micro_w': 'Weighted F-score (Micro)'}
    ontology_dict = {'biological_process': 'BPO', 'molecular_function': 'MFO', 'cellular_component': 'CCO'}
    
    if not os.path.exists(out_path):
        os.mkdir(out_path)
    
    for file in dir_list:
        df_file = results_path + file +"/evaluation_all.tsv"
        df = pd.read_csv(df_file, sep="\t")
        out_folder = out_path + file
        if not os.path.exists(out_folder):
            os.mkdir(out_folder)
            
        
        df = pd.read_csv(df_file, sep="\t")
        
        # Set method information (optional)
        if names_file is None:
            df['group'] = df['filename']
            df['label'] = df['filename']
            df['is_baseline'] = False
        else:
            methods = pd.read_csv(names_file, sep = "\t", header=0)
            df = pd.merge(df, methods, on='filename', how='left')
            df['group'].fillna(df['filename'], inplace=True)
            df['label'].fillna(df['filename'], inplace=True)
            if 'is_baseline' not in df:
                df['is_baseline'] = False
            else:
                df['is_baseline'].fillna(False, inplace=True)
            # print(methods)
        #df = df.drop(columns='filename').set_index(['group', 'label', 'ns', 'tau'])
        df = df.set_index(['group', 'label', 'ns', 'filename','tau'])
        
        # Filter by coverage
        df = df[df['cov'] >= coverage_threshold]
        
        # Assign colors based on group
        cmap = plt.get_cmap('tab20')
        df['colors'] = df.index.get_level_values('group')
        df['colors'] = pd.factorize(df['colors'])[0]
        df['colors'] = df['colors'].apply(lambda x: cmap.colors[x % len(cmap.colors)])
        
        index_best = df.groupby(level=['group', 'ns'])[metric].idxmax() if metric in ['f', 'f_w', 'f_micro', 'f_micro_w'] else df.groupby(['group', 'ns'])[metric].idxmin()
        
        # Filter the dataframe for the best methods
        df_methods = df.reset_index('tau').loc[[ele[:-1] for ele in index_best], ['tau', 'cov', 'colors'] + cols + [metric]].sort_index()

        # Makes the curves monotonic. Cumulative max on the last column of the cols variable, e.g. "pr" --> precision
        if cumulate:
            if metric in ['f', 'f_w', 'f_micro', 'f_micro_w']:
                df_methods[cols[-1]] = df_methods.groupby(level=['label', 'ns'])[cols[-1]].cummax()
            else:
                df_methods[cols[-1]] = df_methods.groupby(level=['label', 'ns'])[cols[-1]].cummin()


        # Save to file
        df_methods.drop(columns=['colors']).to_csv('{}/fig_{}.tsv'.format(out_folder, metric), float_format="%.3f", sep="\t")
        
        # Add first last points to precision and recall curves to improve APS calculation
        def add_points(df_):
            df_ = pd.concat([df_.iloc[0:1], df_])
            df_.iloc[0, df_.columns.get_indexer(['tau', cols[0], cols[1]])] = [0, 1, 0]  # tau, rc, pr
            df_ = pd.concat([df_, df_.iloc[-1:]])
            df_.iloc[-1, df_.columns.get_indexer(['tau', cols[0], cols[1]])] = [1.1, 0, 1]
            return df_

        if metric.startswith('f') and add_extreme_points:
            df_methods = df_methods.reset_index().groupby(['group', 'label', 'ns'], as_index=False).apply(add_points).set_index(['group', 'label', 'ns'])
        
        # Filter the dataframe for the best method and threshold
        df_best = df.loc[index_best, ['cov', 'colors'] + cols + [metric]]
        
        # Calculate average precision score 
        if metric.startswith('f'):
            df_best['aps'] = df_methods.groupby(level=['group', 'label', 'ns'])[[cols[0], cols[1]]].apply(lambda x: (x[cols[0]].diff(-1).shift(1) * x[cols[1]]).sum())

        # Calculate the max coverage across all thresholds
        df_best['max_cov'] = df_methods.groupby(level=['group', 'label', 'ns'])['cov'].max()
        
        # Set a label column for the plot legend
        df_best['label'] = df_best.index.get_level_values('label')
        if 'aps' not in df_best.columns:
            df_best['label'] = df_best.agg(lambda x: f"{x['label']} ({metric.upper()}={x[metric]:.3f} C={x['max_cov']:.3f})", axis=1)
        else:
            df_best['label'] = df_best.agg(lambda x: f"{x['label']} ({metric.upper()}={x[metric]:.3f} APS={x['aps']:.3f} C={x['max_cov']:.3f})", axis=1)
        
        # Generate the figures
        plt.rcParams.update({'font.size': 22, 'legend.fontsize': 18})

        # F-score contour lines
        x = np.arange(0.01, 1, 0.01)
        y = np.arange(0.01, 1, 0.01)
        X, Y = np.meshgrid(x, y)
        Z = 2 * X * Y / (X + Y)

        
        for ns, df_g in df_best.groupby(level='ns'):
            fig, ax = plt.subplots(figsize=(15, 15))

             # Contour lines. At the moment they are provided only for the F-score
            if metric.startswith('f'):
                CS = ax.contour(X, Y, Z, np.arange(0.1, 1.0, 0.1), colors='gray')
                ax.clabel(CS, inline=True) #, fontsize=10)

            cnt = 0
            # Iterate methods
            for i, (index, row) in enumerate(df_g.sort_values(by=[metric, 'max_cov'], ascending=[False if metric.startswith('f') else True, False]).iterrows()):
                
                #data = df_methods.loc[index[:-1]]
                data = df_methods.loc[index[:-2]]
                print(row[cols[0]], row[cols[1]])

                # Precision-recall or mi-ru curves
                ax.plot(data[cols[0]], data[cols[1]], color=row['colors'], label=row['label'], lw=2, zorder=500-i)

                # F-max or S-min dots
                ax.plot(row[cols[0]], row[cols[1]], color=row['colors'], marker='o', markersize=12, mfc='none', zorder=1000-i)
                ax.plot(row[cols[0]], row[cols[1]], color=row['colors'], marker='o', markersize=6, zorder=1000-i)

                cnt+=1
                if n_curves and cnt >= n_curves:
                    break
                
            # Set axes limit
            if metric.startswith('f'):
                plt.xlim(0, 1)
                plt.ylim(0, 1)

            # plt.xlim(0, max(1, df_best.loc[:,:,ns,:][cols[0]].max()))
            # plt.ylim(0, max(1, df_best.loc[:,:,ns,:][cols[1]].max()))

            # Set titles
            ax.set_title(file)
            ax.set_xlabel(axis_title_dict[cols[0]], labelpad=20)
            ax.set_ylabel(axis_title_dict[cols[1]], labelpad=20)

            # Legend
            # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            leg = ax.legend(markerscale=6, title=file)
            for legobj in leg.get_lines():
                legobj.set_linewidth(10.0)
                
            leg.set_bbox_to_anchor((1.05, 1))  

            # Save figure on disk
            plt.savefig("{}/fig_{}_{}.png".format(out_folder, metric, ns), bbox_inches='tight', dpi=300, transparent=True)
            # plt.clf()


