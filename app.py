import streamlit as st
from cyvcf2 import VCF
from config import VCF_FILE, DB_FILE, CCRE_FILE
from utils.gene_utils import load_gene_db, get_gene_coord
from utils.snp_tools import count_snps, variant_consequences
from utils.ccres import load_ccres, ccre_variants
import re
import pandas as pd

# —————————————————————
# Load / cache resources
# —————————————————————
@st.cache_resource
def get_db():
    return load_gene_db(DB_FILE)

@st.cache_resource
def get_vcf():
    return VCF(VCF_FILE)

@st.cache_resource
def get_ccres():
    return load_ccres(CCRE_FILE)

db = get_db()
vcf = get_vcf()
ccres = get_ccres()

# —————————————————————
# UI: Inputs
# —————————————————————
st.title("Mouse Strain SNP Counter")

input_type = st.radio("Lookup by…", ("Gene symbol", "Genomic coordinate"))

if input_type == "Gene symbol":
    gene_names = sorted([
        gene.attributes['gene_name'][0]
        for gene in db.features_of_type('gene')
        if 'gene_name' in gene.attributes
    ])
    gene_in = st.selectbox("Select gene", gene_names, help="Select a gene symbol from the list.")
else:
    coord_in = st.text_input("Enter region (e.g. 1:10000-12000)", help="Format: chr:start-end (e.g. 1:10000-12000)")
    # Input validation for coordinate
    coord_valid = re.match(r'^\w+:\d+-\d+$', coord_in.strip()) if coord_in else False
    if coord_in and not coord_valid:
        st.error("Invalid coordinate format. Please use chr:start-end (e.g. 1:10000-12000).")
        st.stop()

strain_list = vcf.samples
founder_strains = [
    'CAST_EiJ','PWK_PhJ','WSB_EiJ',
    'A_J','NZO_HlLtJ','NOD_ShiLtJ','129S1_SvImJ'
]
ordered_strains = founder_strains + [s for s in strain_list if s not in founder_strains]
strain = st.sidebar.selectbox("Select strain", ordered_strains, help="Choose a mouse strain.")

# —————————————————————
# UI: Actions
# —————————————————————
if st.button("Count SNPs"):
    with st.spinner("Counting SNPs..."):
        if input_type == "Gene symbol":
            region = get_gene_coord(db, gene_in)
            if region is None:
                st.error(f"Gene '{gene_in}' not found.")
                st.stop()
        else:
            region = coord_in

        idx = vcf.samples.index(strain)
        hom, het, var_list = count_snps(vcf, region, idx)
        ccre_table = ccre_variants(var_list, region, ccres)

        st.metric("Homozygous SNPs", hom)
        st.metric("Heterozygous SNPs", het)

        if input_type == "Gene symbol":
            var_conseq = variant_consequences(vcf, region, idx)
            st.subheader("Variant consequences:")
            st.table(var_conseq)
            st.subheader("Variant consequence distribution")
            st.bar_chart(var_conseq)
            # Download button for variant consequences
            st.download_button(
                "Download variant consequences as CSV",
                data=pd.DataFrame(var_conseq).to_csv(index=False),
                file_name="variant_consequences.csv",
                mime="text/csv"
            )

        st.subheader("cCRE overlap with variants:")
        st.table(ccre_table)
        # Download button for cCRE table
        st.download_button(
            "Download cCRE overlap table as CSV",
            data=pd.DataFrame(ccre_table).to_csv(index=False),
            file_name="ccre_overlap.csv",
            mime="text/csv"
        )

        