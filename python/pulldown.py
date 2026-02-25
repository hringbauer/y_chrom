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
    """Creates the Call Table from a .bam file"""
    
    ### Create the Pulldown
    pulldown_bamtable(path_bam = path_bam,
                      path_bed = path_bed, 
                      snip5=snip5, snip3=snip3,
                      o_file = path_temp)

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