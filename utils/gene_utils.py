import gffutils
import sqlite3

def load_gene_db(db_path):
    conn = sqlite3.connect(db_path, check_same_thread=False)
    return gffutils.FeatureDB(conn, keep_order=True)

def get_gene_coord(db, gene_symbol):
    genes = [
        g for g in db.features_of_type("gene")
        if g.attributes.get("gene_name", [None])[0] == gene_symbol
    ]
    if not genes:
        return None, None
    gene = genes[0]
    return f"{gene.chrom}:{gene.start}-{gene.end}"
