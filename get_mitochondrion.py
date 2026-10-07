# ex: python3 get_mitochondrion.py MK089295 Rhinolophus_ferrumequinum
from Bio import SeqIO, Entrez
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import Counter
import sys, os

# ── Input ─────────────────────────────────────────────────────────────────────
accession = sys.argv[1]
species   = sys.argv[2] if len(sys.argv) > 2 else accession
basename  = f"{species}_{accession}"

Entrez.email = "nicolas.nesi@unicaen.fr"

# ── Download or load local gb ─────────────────────────────────────────────────
if os.path.exists(f"{basename}.gb"):
    print(f"Found local {basename}.gb — skipping download")
    record = SeqIO.read(f"{basename}.gb", "genbank")
else:
    print(f"Downloading {accession} from NCBI...")
    with Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text") as handle:
        record = SeqIO.read(handle, "genbank")
    with open(f"{basename}.gb", "w") as f:
        SeqIO.write(record, f, "genbank")

chrom = record.id
print(f"Record: {chrom} — {record.description}")

# ── Map feature type to gene_biotype ─────────────────────────────────────────
biotype_map = {
    "CDS"  : "protein_coding",
    "tRNA" : "tRNA",
    "rRNA" : "rRNA"
}

# ── Pass 1: collect raw features and count name occurrences ──────────────────
raw_features = []
for feature in record.features:
    if feature.type not in ("CDS", "tRNA", "rRNA"):
        continue
    name = feature.qualifiers.get("gene",
           feature.qualifiers.get("product", ["unknown"]))[0]
    name = name.replace(" ", "_")
    name = f"mt-{name}"
    raw_features.append((name, feature))

name_totals = Counter(name for name, _ in raw_features)

# ── Pass 2: build final (possibly disambiguated) names, in encounter order ───
name_running_count = Counter()
fasta_records = []
gtf_lines     = []

for name, feature in raw_features:
    if name_totals[name] > 1:
        name_running_count[name] += 1
        final_name = f"{name}{name_running_count[name]}"   # mt-tRNA-Leu1, mt-tRNA-Leu2, ...
    else:
        final_name = name

    start   = int(feature.location.start)
    end     = int(feature.location.end)
    strand  = "+" if feature.location.strand == 1 else "-"
    biotype = biotype_map[feature.type]

    seq = feature.extract(record.seq)
    fasta_records.append(SeqRecord(seq, id=final_name, description=""))

    attr = f'gene_id "{final_name}"; transcript_id "{final_name}"; gene_biotype "{biotype}";'
    gtf_lines.append(f"{final_name}\tGenBank\t{feature.type}\t{start+1}\t{end}\t.\t{strand}\t.\t{attr}")

# ── Write outputs ─────────────────────────────────────────────────────────────
with open(f"{basename}_genes.fa", "w") as f:
    SeqIO.write(fasta_records, f, "fasta-2line")

with open(f"{basename}_genes.gtf", "w") as f:
    f.write("\n".join(gtf_lines) + "\n")

# ── Summary ───────────────────────────────────────────────────────────────────
counts = Counter(biotype_map[f.type] for f in record.features if f.type in biotype_map)
print(f"Written : {len(fasta_records)} features")
print(f"Fasta   : {basename}_genes.fa")
print(f"GTF     : {basename}_genes.gtf")
for biotype, count in counts.items():
    print(f"  {biotype:20s}: {count}")

duplicated = {n: c for n, c in name_totals.items() if c > 1}
if duplicated:
    print("Disambiguated duplicate names:")
    for n, c in duplicated.items():
        print(f"  {n} → {n}1..{n}{c}")