import streamlit as st
from cyvcf2 import VCF
import gffutils
import pandas as pd
import pyranges as pr
# —————————————————————
# 1) Load / cache your resources
# —————————————————————

# VCF_FILE = "/Users/ryan/Desktop/variant_streamlit/files/mgp_REL2021_snps.CHROM1_sub.vcf.gz"
VCF_FILE = "/Users/ryan/Desktop/variant_streamlit/files/mgp_REL2021_snps.vcf.gz"
# VCF_FILE = "/Users/ryan/Desktop/variant_streamlit/files/mgp_REL2021_snps.vcf.gz"

# DB_FILE = "/Users/ryan/Desktop/variant_streamlit/mm39_gtf_full_test_2.db"
DB_FILE = "/Users/ryan/Desktop/variant_streamlit/mm39_gtf_filt.db"
# DB_FILE = "/Users/ryan/Desktop/variant_streamlit/files/mm39_gtf.db"
cCRE_FILE = "/Users/ryan/Desktop/variant_streamlit/files/mm39_ccres.bed"



@st.cache_resource
def load_db(db_path=DB_FILE):
    # will build or open your gffutils SQLite DB exactly once
    return gffutils.FeatureDB(db_path)

@st.cache_resource
def load_vcf(vcf_path=VCF_FILE):
    # will open and index your VCF exactly once
    return VCF(vcf_path)

@st.cache_resource
def load_ccres(bed_path="mm39_ccres.bed"):
    ccres = pr.read_bed(bed_path)
    # strip off any leading "chr"
    ccres.Chromosome = ccres.Chromosome.astype(str).str.replace(r"^chr", "", regex=True)
    return ccres



# then call them
db  = load_db()
vcf = load_vcf()
ccres = load_ccres(cCRE_FILE)

# ————————————————
# 2) SNP‐counting function
# ————————————————
def get_gene_coord(gene_symbol):
    genes = [
        g for g in db.features_of_type("gene")
        if g.attributes.get("gene_name", [None])[0] == gene_symbol
    ]
    if not genes:
        return None, None
    gene = genes[0]
    region = f"{gene.chrom}:{gene.start}-{gene.end}"
    
    return region


def count_snps(gene_symbol, strain):

    region = get_gene_coord(gene_symbol)

    # get sample index
    idx = vcf.samples.index(strain)

    hom = het = 0
    catch_variants=[]
    for var in vcf(region):
        
        variant_coord = var.CHROM + ':' + str(var.POS) + '-' + str(var.POS)
        
        
        gt = var.gt_types[idx]  # 0=hom-ref,1=het,2=hom-alt,3=./.
        if gt == 3:
            hom += 1
            catch_variants.append(variant_coord)

        elif gt == 1:
            het += 1
            catch_variants.append(variant_coord)


    return hom, het, catch_variants, region



def variant_consequences(gene_symbol, strain):
    region = get_gene_coord(gene_symbol)
    idx    = vcf.samples.index(strain)

    cons = []
    for var in vcf(region):
        csq_raw = var.INFO.get("CSQ") or ""
        entries = csq_raw.split(",") if isinstance(csq_raw, str) else csq_raw
        first   = entries[0]
        fields  = first.split("|")
        consequence = fields[1]

        gt = var.gt_types[idx]
        if gt == 3 or gt == 1:
            cons.append(consequence)

    # Now compute the counts once, after the loop
    if cons:
        # Using a pandas Series is easiest
        return pd.Series(cons).value_counts()
    else:
        # Return an empty Series so the caller can still loop or table it
        return pd.Series(dtype=int)



def count_snps_region(region: str, strain: str):
    idx = vcf.samples.index(strain)
    hom = het = 0
    
    catch_variants=[]
    
    for var in vcf(region):
        
        variant_coord = var.CHROM + ':' + str(var.POS) + '-' + str(var.POS)
        
        gt = var.gt_types[idx]
        if gt == 2:
            hom += 1
            catch_variants.append(variant_coord)
            
        elif gt == 1:
            het += 1
            catch_variants.append(variant_coord)
            
    return hom, het, catch_variants



def ccre_variants(var_list, region):
    
    chrom, coords = region.split(":")
    start, end   = map(int, coords.split("-"))
    region_pr = pr.PyRanges(
        chromosomes=[chrom],
        starts=[start],
        ends=[end]
    )
    ccres_in_region = ccres.overlap(region_pr)

    coords = [v.split(":",1)[1].split("-",1) for v in var_list]
    chroms = [v.split(":",1)[0] for v in var_list]
    starts = [int(s) for s,_ in coords]
    ends   = [int(e) for _,e in coords]

    variants_pr = pr.PyRanges(
        chromosomes=chroms,
        starts=starts,
        ends=ends
    )
    var_ccre_hits = variants_pr.join(ccres_in_region)
    
    df = var_ccre_hits.df   # a pandas.DataFrame
    if not df.empty:
        counts = (
          df.groupby("Name")
            .size()
            .reset_index(name="variant_count")
        )
        # optionally, parse cCRE type out of the Name column:
        counts["type"] = counts["Name"].str.split("/").str[0]
    else:
        counts = pd.DataFrame(columns=["Name","variant_count","type"])
        
    return counts



# Extract gene names from gtf db
genes = db.features_of_type('gene')
gene_names = [
    gene.attributes['gene_name'][0]
    for gene in db.features_of_type('gene')
    if 'gene_name' in gene.attributes
    ]
gene_names.sort()

founder_strains = [
    'CAST_EiJ','PWK_PhJ','WSB_EiJ',
    'A_J','NZO_HlLtJ','NOD_ShiLtJ','129S1_SvImJ'
]
strain_list = vcf.samples  # your full list

# Build the new list
ordered_strains = founder_strains + [
    s for s in strain_list
    if s not in founder_strains
]

# ————————————————
# 3) Streamlit UI
# ————————————————
st.title("Mouse Strain SNP Counter")

input_type = st.radio("Lookup by…", ("Gene symbol", "Genomic coordinate"))

if input_type == "Gene symbol":
    gene_in = st.selectbox("Enter gene symbol (e.g. Trp53)", gene_names)
else:
    coord_in = st.text_input("Enter region (1:10000-12000)")

strain = st.selectbox("Select strain", ordered_strains)

if st.button("Count SNPs"):
    if input_type == "Gene symbol":
        hom, het, var_list, region_coords = count_snps(gene_in, strain)
        if hom is None:
            st.error(f"Gene “{gene_in}” not found.")
            st.stop()
        ccre_table = ccre_variants(var_list, region_coords)
            

        

    else:
        hom, het, var_list = count_snps_region(coord_in, strain)
        region_coords = coord_in
        ccre_table = ccre_variants(var_list, region_coords)

        
        
        
    # 4) Display
    st.metric("Homozygous SNPs", hom)
    st.metric("Heterozygous SNPs", het)
    
    if input_type == "Gene symbol":
        var_conseq = variant_consequences(gene_in, strain)
        st.subheader("Variant consequences:")
        st.table(var_conseq)
        
        st.subheader("Variant consequence distribution")
        st.bar_chart(var_conseq)      # works directly if var_conseq is a pandas Series
        
        st.subheader("cCRE overlap with variants:")
        st.table(ccre_table)


        
    elif input_type == "Genomic coordinate":
        st.subheader("cCRE overlap with variants:")
        st.table(ccre_table)



