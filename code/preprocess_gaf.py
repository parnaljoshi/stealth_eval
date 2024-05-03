import numpy as np
import pandas as pd
import gzip
import csv
import io
import argparse

"""
Extracts columns in the extract_col_list from the goa_file and returns them as a Dataframe
Ignores the lines that start with '!'

Parameters:
    - goa_file: The input GOA file name.
    - extract_col_list: Names of the columns to be extracted. Default: ['DB Object ID', 'Qualifier','GO ID', 'Evidence Code', 'Aspect']
    
Returns Extract_ann dataframe of the extracted annotations
"""

def extract_annot(goa_file, extract_col_list = ['DB Object ID', 'Qualifier','GO ID', 'Evidence Code', 'Aspect']):
    # Uniprot gaf file description: https://geneontology.org/docs/go-annotation-file-gaf-format-2.0/
    # List of Gaf_cols in gaf-version 2.0, 2.1 and 2.2
    Gaf_cols = ['DB', 'DB Object ID', 'DB Object Symbol', 'Qualifier', 'GO ID', 'DB:Reference (|DB:Reference)', 
                       'Evidence Code', 'With (or) From', 'Aspect', 'DB Object Name', 'DB Object Synonym (|Synonym)', 'DB Object Type', 
                       'Taxon(|taxon)', 'Date', 'Assigned By', 'Annotation Extension', 'Gene Product Form ID']
    
    
    extract_col_ind = [Gaf_cols.index(col) for col in extract_col_list if col in Gaf_cols]
    print("Indices of the extracted columns are : ", extract_col_ind)
    
    rows = []
    
    with gzip.open(goa_file, 'rt') as f:
        # Skip lines starting with '!'
        filtered_lines = (line for line in f if not line.startswith('!'))
        gaf_v = f.readline().strip().split(" ")
        
        # Join the filtered lines into a single string
        joined_lines = ''.join(filtered_lines)
        
        # Use StringIO to create a file-like object
        file_like_object = io.StringIO(joined_lines)
        
        # Use csv.reader to parse TSV
        tsv_reader = csv.reader(file_like_object, delimiter='\t')
        
        # Iterate over the rows and extract specified columns
        for row in tsv_reader:
            extracted_row = [row[i] for i in extract_col_ind]
            rows.append(extracted_row)
    
    # Create a DataFrame from the extracted rows
    Extracted_ann = pd.DataFrame(rows, columns=extract_col_list)
    
    return Extracted_ann

def remove_dup_and_neg(Extracted_ann, remove_dup = True, remove_neg = True):
    N = len(Extracted_ann)
    print(N, " Rows in the input file")
    if remove_dup:
        Extracted_ann = Extracted_ann.drop_duplicates().copy()
        print(N-len(Extracted_ann), " duplicates dropped")
    if remove_neg:
        Not_qualifier = Extracted_ann["Qualifier"].apply(lambda x:"NOT" in x)
        Extracted_ann = Extracted_ann[~Not_qualifier].copy()
        print(sum(Not_qualifier), " Not Qualifiers found")
    return Extracted_ann
	    

# Documentation of evidence codes - https://geneontology.org/docs/guide-go-evidence-codes/
def filter_evidence_codes(Extracted_ann, evidence_codes = ['EXP', 'IDA', 'IPI','IMP', 'IGI', 'IEP', 'TAS', 'IC'], highTP = False):
    if highTP:
        evidence_codes += ['HTP', 'HDA', 'HMP', 'HGI', 'HEP']
        print("High Thoughput Evidence Code included")
    print("Included evidence codes are :", evidence_codes)
    evidence_True = Extracted_ann['Evidence Code'].isin(evidence_codes)
    print(sum(evidence_True))
    Filtered = Extracted_ann[evidence_True].copy() 
    return Filtered

def write_annot(Extracted_ann, out_file = "extracted.tsv", out_cols = ['DB Object ID', 'GO ID', 'Aspect']):
    print(Extracted_ann)
    print(Extracted_ann.columns)
    Extracted_ann.loc[:, out_cols].to_csv(out_file, index = False, sep = "\t")
    

def parse_args():
    parser = argparse.ArgumentParser(description='Process the arguments')
    parser.add_argument("goa_path", help="Path of the gaf.gz file")
    parser.add_argument("--extract_col_list", help="List of the columns to be extracted, default columns are - 'DB Object ID', 'Qualifier','GO ID', 'Evidence Code', 'Aspect'", nargs = '+', default = ['DB Object ID', 'Qualifier','GO ID', 'Evidence Code', 'Aspect'])
    parser.add_argument('--no_dup', action='store_true', help='Enables removal of duplicates, default = True', default = True)
    parser.add_argument('--no_neg', action='store_true', help='Enables removal of annotations with negative Qualifier, default = True', default = True) 
    parser.add_argument("--evidence_codes", help="List of the evidence codes to be included, default codes are - ['EXP', 'IDA', 'IPI','IMP', 'IGI', 'IEP', 'TAS', 'IC']", nargs = '+', default = ['EXP', 'IDA', 'IPI','IMP', 'IGI', 'IEP', 'TAS', 'IC'])
    parser.add_argument("--highTP", action='store_true', help="Include high throughput evidence codes - 'HTP', 'HDA', 'HMP', 'HGI', 'HEP'", default = False)
    parser.add_argument("--out_path", type=str, help="Path of the extracted file, default = extracted.tsv", default = "extracted.tsv")
    parser.add_argument("--only_annot", action='store_true', help="Only outputs the DB Object ID, GO ID and Aspect in the out_file, default = True", default = True)
    args = parser.parse_args()
    return args
		
def main():
    args = parse_args()
    print(args)

    # Extract the required columns from the gaf file
    Extracted_ann = extract_annot(args.goa_path, args.extract_col_list)
    
    # Filter the annotations of the required evidence codes
    Extracted_f = filter_evidence_codes(Extracted_ann, args.evidence_codes, args.highTP)
    
    # Drop duplicates and annotations with negative qualifier
    Extracted = remove_dup_and_neg(Extracted_f, args.no_dup, args.no_neg)
    
    
    # Determine whether the user want to print only the annotation fields, or all the fields passed in the extract_col_list
    if args.only_annot == True:
        out_cols = ['DB Object ID', 'GO ID', 'Aspect']
    else:
        out_cols = args.extract_col_list
        
    # Write the filtered annotations to out_file
    write_annot(Extracted, args.out_path, out_cols)
	
if __name__ == '__main__':
    main()
    print("Done")
