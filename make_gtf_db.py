import gffutils

# Constants: adjust paths as needed
GTF_FILE = "gencode.vM37.nochr.annotation.gtf"
DB_FILE = "mm39_gtf_filt.db"

gffutils.create_db(
    data=GTF_FILE,
    dbfn=DB_FILE,
    force=True,               # overwrite if it already exists
    keep_order=True,          # preserve GTF order (helps on some queries)
    merge_strategy="merge",   # combine attributes with the same key
    sort_attribute_values=True
)
