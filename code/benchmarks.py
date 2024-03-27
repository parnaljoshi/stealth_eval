from InformationAccretion.ia import *
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
                annot_gained[protein][ont] = t1_annot_list[protein][ont]
    return annot_gained


def get_baselines(t0_annot_list, t1_annot_list):
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
        
        type3_df = flatten_nested_dict(type3_annot_gained[prot][ont])
            
        type3_df.to_csv('/home/rashika/CAFA4/eval/benchmarks_GO/Type3_ann.csv', index=False)
        
        # type 12 - union of type 1 and type 2
        for ont in ['BPO', 'CCO', 'MFO']:
            type12[ont] = list(set(type1[ont] + type2[ont]))
        
    return type1, type2, type3, type12

def flatten_nested_dict(nested_dict):
    flat_dict = [[outer_key, inner_key,value] for outer_key, inner_dict in nested_dict.items() for inner_key, values in inner_dict.items() for value in values]
    df = pd.DataFrame(flat_dict)
    df.columns=['EntryID','aspect','term']
    df = df.reindex(columns=['EntryID','term','aspect'])
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
    
    
def create_bm_lists(t0_file, t1_file, t0_ont_graph, t1_ont_graph, roots, BM_path = "/home/rashika/CAFA4/eval/benchmarks/"):
    
    #Prop t0 and t1 in their respective ontologies
    t0_prop = process_raw_annot(t0_file, t0_ont_graph, roots)
    t1_prop = process_raw_annot(t1_file, t1_ont_graph, roots)
    
    # Keep common terms
    t0_common, t1_common =  keep_common_go_terms(t0_prop, t1_prop, t0_ont_graph, t1_ont_graph)

    # Propagate back in the t0 ontology
    subontologies = {aspect: fetch_aspect(t0_ont_graph, roots[aspect]) for aspect in roots}
    t0_eval = propagate_terms(t0_common, subontologies)
    t1_eval = propagate_terms(t1_common, subontologies)

    # Convert the eval Dfs into annotation lists
    t0_annot_list = get_annot_list(t0_eval)
    t1_annot_list = get_annot_list(t1_eval)
    
    # Convert the eval Dfs into annotation lists
    type1, type2, type3, type1 = get_baselines(t0_annot_list, t1_annot_list)
    
    # Check if the BM_path exists, create the directory if it does not
    if not os.path.exists(BM_path):
        os.mkdir(BM_path)
    
    
    for ont in ["BPO", "CCO", "MFO"]:
        for bl_type in [type1, type2, type3, type12]:
            bl_file_path = BM_path + ont.lower() + "_all" + "_"+ get_variable_name(bl_type, locals()) +".txt"
            bl_list = bl_type[ont]
            print(bl_file_path)
            with open(bl_file_path, "w") as outfile:
                outfile.write("\n".join(bl_list))
    
    return t1_eval

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

