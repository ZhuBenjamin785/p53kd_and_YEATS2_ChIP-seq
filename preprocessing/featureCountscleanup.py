import pandas as pd
import sys
import re

                        
if len(sys.argv) not in (3, 4):
    print("Usage: python script_name.py <input_file> <output_file> [annotation_file]")
    sys.exit(1)

                               
input_file = sys.argv[1]
output_file = sys.argv[2]
annotation_file = sys.argv[3] if len(sys.argv) == 4 else None


def read_gene_names(path):
    """Read gene_name from the project BED/GTF-like annotation."""
    if not path:
        return {}
    names = {}
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            gene_id = fields[3]
            match = re.search(r'gene_name\s+"([^"]+)"', fields[9])
            if match:
                names[gene_id] = match.group(1)
                names.setdefault(gene_id.split(".", 1)[0], match.group(1))
    return names

                     
fc = pd.read_csv(input_file, sep="\t", comment="#")

                                        
fc = fc.drop(["Chr", "Start", "End", "Strand"], axis=1, errors="ignore")

                                               
renamed_columns = []
for column in fc.columns:
    if column.startswith("BAMfiles/"):
        renamed_columns.append(column.split("/")[-1].split("_S")[0])
    elif column.endswith(".sorted.bam") or column.endswith(".bam"):
        sample = column.rsplit("/", 1)[-1]
        if sample.endswith(".sorted.bam"):
            sample = sample[:-len(".sorted.bam")]
        else:
            sample = sample[:-len(".bam")]
        renamed_columns.append(sample)
    else:
        renamed_columns.append(column)
fc.columns = renamed_columns

                                                                            
                                                                     
gene_names = read_gene_names(annotation_file)
if gene_names and "Geneid" in fc.columns:
    def display_name(gene_id):
        gene_id = str(gene_id)
        return gene_names.get(gene_id, gene_names.get(gene_id.split(".", 1)[0], gene_id))

    fc["Geneid"] = fc["Geneid"].map(display_name)
    fc = fc.rename(columns={"Geneid": "gene_name"})

                                                 
fc.to_csv(output_file, sep="\t", index=False)

        
