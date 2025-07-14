from cyvcf2 import VCF
import pandas as pd

def count_snps(vcf, region, strain_idx):
    hom = het = 0
    catch_variants = []
    for var in vcf(region):
        coord = f"{var.CHROM}:{var.POS}-{var.POS}"
        gt = var.gt_types[strain_idx]
        if gt == 3:
            hom += 1
            catch_variants.append(coord)
        elif gt == 1:
            het += 1
            catch_variants.append(coord)
    return hom, het, catch_variants

def variant_consequences(vcf, region, strain_idx):
    cons = []
    for var in vcf(region):
        csq_raw = var.INFO.get("CSQ") or ""
        entries = csq_raw.split(",")
        fields = entries[0].split("|") if entries else []
        if len(fields) > 1:
            consequence = fields[1]
            gt = var.gt_types[strain_idx]
            if gt == 3 or gt == 1:
                cons.append(consequence)
    return pd.Series(cons).value_counts() if cons else pd.Series(dtype=int)
