import pyranges as pr
import pandas as pd

def load_ccres(bed_path):
    ccres = pr.read_bed(bed_path)
    ccres.Chromosome = ccres.Chromosome.astype(str).str.replace(r"^chr", "", regex=True)
    return ccres

def ccre_variants(var_list, region, ccres):
    chrom, coords = region.split(":")
    start, end = map(int, coords.split("-"))
    region_pr = pr.PyRanges(chromosomes=[chrom], starts=[start], ends=[end])
    ccres_in_region = ccres.overlap(region_pr)

    coords = [v.split(":",1)[1].split("-",1) for v in var_list]
    chroms = [v.split(":",1)[0] for v in var_list]
    starts = [int(s) for s,_ in coords]
    ends = [int(e) for _,e in coords]

    variants_pr = pr.PyRanges(chromosomes=chroms, starts=starts, ends=ends)
    hits = variants_pr.join(ccres_in_region)

    df = hits.df
    if not df.empty:
        counts = df.groupby("Name").size().reset_index(name="variant_count")
        counts["type"] = counts["Name"].str.split("/").str[0]
    else:
        counts = pd.DataFrame(columns=["Name", "variant_count", "type"])
    return counts
