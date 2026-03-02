import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os as os  # For Saving to Folder

def load_counts(path_counts, coerce=True):
    """Load Count file and return Dataframe"""
    df_t = pd.read_csv(path_counts, header=None, delim_whitespace=True)
    df_t.columns = ["snp", "chr", "pos", "ref_all", "alt_all", "drop", "iid", "ref", "alt"]
    
    if coerce:
        for col in ["pos", "ref", "alt"]:
            df_t[col] = pd.to_numeric(df_t[col], errors="coerce")
            
    df_t = df_t.drop(columns="drop")
    return df_t

def load_snp_file_ISOGG(path_snps = "./data/all_snps.csv", 
                    col_pos = 'Build 37 Number', unique=True):
    """Return Dataframe in Eigenstrat Format,
    filtered for biallelic SNPs.
    unique: Whether to keep a """
    df_raw = pd.read_csv(path_snps)
    print(df_raw.columns)
    print(f"Loaded {len(df_raw)} SNPs")

    ### Process the positions
    pos = df_raw[col_pos]
    df_raw["pos"] = pd.to_numeric(pos, errors="coerce")

    idx = ~df_raw["pos"].isna()
    print(f"# Positions available: {np.sum(idx)}")
    df = df_raw[idx].reset_index(drop=True)
    df["pos"]=df["pos"].astype("int")

    idx_bi= (df["Mutation Info"].str.len()==4)
    print(f"# Biallelic SNPs: {np.sum(idx_bi)}")
    df = df[idx_bi].reset_index(drop=True)
    df["ref"] = df["Mutation Info"].str[0]
    df["alt"] = df["Mutation Info"].str[3]
    df["chrom"] = "Y"

    cols = ["Name", "chrom", "pos", "ref", "alt", 
            'Subgroup Name', 'Alternate Names', 'rs numbers']
    df = df[cols]
    df = df.replace(regex=[' ','\n'], value='_')
    ### Sort by position
    df = df.sort_values(by="pos")
    
    ### Keep only SNPs where Ref and Alt Different
    idx_same = (df["ref"]==df["alt"])
    df = df[~idx_same]
    print(f"# Ref & Alt different: {len(df)}")
    
    ### Keep only ACTG
    snps_acceptable = ["A", "C", "T", "G"]
    idx_ref = df["ref"].isin(snps_acceptable)
    idx_alt = df["alt"].isin(snps_acceptable)
    idx_both = idx_ref & idx_alt
    df = df[idx_both]
    print(f"# Ref & Alt ACTG: {len(df)}")
    
    ### Keep Unique Values
    if unique:
        idx_dup = df.duplicated(subset=["pos", "ref", "alt"], keep="first")
        df = df[~idx_dup]
        print(f"# Unique SNP positions: {len(df)}")
    
    ### Remove duplicate Names
    #idx_dup = df.duplicated(subset="Name", keep=False)
    #df = df[~idx_dup]
    #print(f"# Unique Names: {len(df)}")
    return df.copy().reset_index(drop=True)


################################################
### Calling Ys

def ref_alt_count(df_ch, bases=["A", "C", "G", "T"]):
    """Count Ref and Alt alleles in Dataframe df_ch
    with ref, alt, A, C, G, T fields and enter new columns
    ref# and alt#"""
    df_ch["ref#"]=0
    df_ch["alt#"]=0

    for p in bases:
        idx = df_ch["ref"] == p
        df_ch.loc[idx, "ref#"] = df_ch.loc[idx, p]

        idx = df_ch["alt"] == p
        df_ch.loc[idx, "alt#"] = df_ch.loc[idx, p]
    return df_ch

def pulldown_bamtable(path_bam = "", o_file = "",                   
                      bamtable = "/home/pruefer/bin/BamTable",
                      snip5=0, snip3=0,
                      path_bed = "/mnt/archgen/users/hringbauer/git/y_chrom/data/isogg_snps.bed"):
    """Pulldown a BAM at path_bam to o_file using bamtable and the bed a path_bed."""
    run_cmd = f"{bamtable} -F -A --snip5={snip5} --snip3={snip3} -f {path_bed} {path_bam} > {o_file}"
    os.system(run_cmd)

def call_y_bam(path_bam="", df=[],
               path_bed = "/mnt/archgen/users/hringbauer/git/y_chrom/data/isogg_snps.bed",
               path_temp="/mnt/archgen/users/hringbauer/git/y_chrom/temp/temp.tsv",
               snip5=0, snip3=0):
    """Return the Call Table from a .bam file"""
    ### Sanity Checks whether Input is there
    assert(os.path.exists(path_bam))
    assert(os.path.exists(path_bed))
    
    ### Create the Pulldown
    pulldown_bamtable(path_bam = path_bam,
                      path_bed = path_bed, 
                      snip5=snip5, snip3=snip3,
                      o_file = path_temp)

    ### Load and formt pulldown
    df1 = pd.read_csv(path_temp, sep="\t", header=None)
    df1.columns = ["chrom", "pos", "A", "C", "G", "T"]
    idx = df1["chrom"]=="chrY"
    if np.sum(idx)>0:
        print(f"Changing {np.sum(idx)} ChrY -> Y")
        df1.loc[idx, "chrom"] = "Y"
        
    df2 = pd.merge(df, df1, on=["chrom", "pos"])
    
    ### Coverage Statistics
    cov = df1[["A", "C", "G", "T"]].values
    cov1 = np.sum(cov, axis=1)
    print(f"Average Coverage: {np.sum(cov1)/len(df):.4f}x")
    print(f"#Sites covered: {np.sum(cov1>0)}/{len(df)}")
    
    ### Establish Ref and Alt allele
    df_ch = ref_alt_count(df2, bases=["A", "C", "G", "T"])

    ### Identify Derived    
    idx_der = df_ch["alt#"]>df_ch["ref#"]
    print(f"#Derived Loci: \n{np.sum(idx_der)} / {np.sum(cov1>0)} covered>0")
    
    df_der = df_ch[idx_der].sort_values(by="Subgroup Name").reset_index(drop=True).copy()
    
    return df_ch, df_der 

def mismatch_path(s, df):
    """Look for all mismatches in path up to s"""
    ls = [s, s+"~"]

    for i in range(1,len(s)):
        ls+= [s[:-i], s[:-i]+"~"]

    dft = df[df["Subgroup Name"].isin(ls)]
    idx = dft["ref#"]>=dft["alt#"]
    print(f"Mismatches: {np.sum(idx)} / {len(idx)}")
    return dft

###########################
### Functions for browsing Haplogroups

def find_supgroups_der(df, subgroup="", col="Subgroup Name"):
    """Output all real subgroups (but not the orginal one)."""
    dft=df[df[col].str.contains(subgroup)]
    dft2=dft[dft[col]!=subgroup]
    return dft2

def calc_error_rate(df_ch, haploclade="G"):
    """ Calculate Error Rate based on supposedly ancestral clade"""
    dft = df_ch[df_ch["Subgroup Name"].str.startswith(haploclade)]
    der = np.sum(dft["alt#"])
    ref = np.sum(dft["ref#"])
    print((ref,der))
    print(f"Error rate: {der/(ref+der)*100:.4g}%")


################################
### Functions for OY database

### Prepare Dictionary of Levels

def create_parent_dct(path_parents="/mnt/archgen/users/eric_garcia/OYdb/OYchpar.csv"):
    """"Create Dictionary of Parent Nodes"""
    chpar = {}
    
    # Create a dictionary to store all child-parent relations.
    with open(path_parents, "r") as f:
        for line in f:
            items = line.strip().split(sep=",")
            if not items:
                continue
                    
            child = items[0]
            parent = items[1]
                
            chpar[child] = parent
    
    return chpar


# Create a function to get the total number of SNPs per branch and divide between derived and ancestral.
def div_anc_der(df, chpar=[], df_exclude=[]):
    """
    Calculate total, ancestral, derived, and uncovered SNPs per branch.
    df_exclude: Dataframe of SNPs to filter
    """
    if len(df_exclude)>0:
        df = exclude_snps(df, df_exclude)
    
    if len(chpar)==0:
        chpar = create_parent_dct()
    
    # Create columns for ANC and DER SNPs.
    df = df.copy()
    df["ancestral"] = df["alt#"] < df["ref#"]
    df["derived"]   = df["alt#"] > df["ref#"]

    # Divide the df by the different haplogroups.
    grouped = df.groupby("Y-haplogroup")

    # Generate a new df with columns for the Branch, the Level and the total for SNPs.
    new_df = grouped.agg(
        Level=("Level", "first"),
        Total_SNPs=("Y-haplogroup", "size"),
        Ancestral=("ancestral", "sum"),
        Derived=("derived", "sum"),
    ).reset_index()

    # Compute uncovered as the difference of the total with ancestral and derived.
    new_df["Uncovered"] = (
        new_df["Total_SNPs"]
        - new_df["Ancestral"]
        - new_df["Derived"]
    )

    # Rename column called Y-haplogroup
    new_df = new_df.rename(columns={"Y-haplogroup": "Branch"})
    new_df = new_df.sort_values(by="Level")

    # Use previously defined function to look for all ANC and DER SNPs in parental branches. First, transform into a dictionary to speed up processing.
    anc_lookup = new_df.set_index("Branch")["Ancestral"].to_dict()
    der_lookup = new_df.set_index("Branch")["Derived"].to_dict()
    
    anc_par = []
    der_par = []

    # For each branch, find the total number of ancestral and derived SNPs in parent using ancder_par() function.
    for branch in new_df["Branch"]:
        anc_par.append(ancder_par(branch, anc_lookup, chpar=chpar))
        der_par.append(ancder_par(branch, der_lookup, chpar=chpar))

    # Insert new columns.
    new_df["#ANC in par."] = anc_par
    new_df["#DER in par."] = der_par
    
    # Save a csv file
    #new_df.to_csv(f"data/branches_statitstics_{sample}.csv")

    return new_df

def ancder_par(string, par_dict, chpar):
    """ Find the total number of ancestral and derived SNPs for all parental branches"""

    if string not in chpar:
        return 0

    else:
        parent = chpar.get(string)
        return ancder_par(parent, par_dict, chpar) + par_dict.get(parent, 0)

def get_mismatch_snps(string, chpar, df_ch):
    """Get all SNPs that are mismatches of Haplogroup String"""
    dfs_mm = [] # List of mismatching SNP dfs
    
    while True:        
        ### Find all mismatches
        dft = df_ch[df_ch["Y-haplogroup"]==string] # All SNPs in Node
        dfd =dft[dft["ref#"]>dft["alt#"]]
        dfs_mm.append(dfd)  	

        if string not in chpar: # If Root exit
            dfs_mm = pd.concat(dfs_mm)
            return dfs_mm
        
        else:
            string = chpar[string] # Get Parent Node


def exclude_snps(df_ch, df_ex, verbose=False, col="Subgroup Name"):
    """Filter SNPs from df_ch in df_ex"""

    idx= df_ch[col].isin(df_ex[col])
    df_ch2 = df_ch[~idx].copy()
    if verbose:
        print(f"Filtered to {len(df_ch2)}/{len(df_ch)}")
    return df_ch2



    