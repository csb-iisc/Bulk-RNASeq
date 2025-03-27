#!/bin/bash
# Genome metadata
genome_path='/mnt/nvme_sequencing/RNASeq'
thr=110

# Human genome
# Using release 109 FASTA for GRCh38 instead of release 110 FASTA -- release 110 moved from GRCh38.p13 to GRCh38.p14,
# which unmasked the pseudo-autosomal region. This causes ambiguous mappings to PAR locus genes.
# No other sequence changes were made to the primary assembly.
# org="hg38"
# genome_release="109"
# gencode_release="44"

ulimit -n 999999
# # Mouse genome
org="mm10"
genome_release="110"
gencode_release="M33"

# Set up source and build directories
mkdir -p "${genome_path}/${org}-${genome_release}"
mkdir -p "${genome_path}/${org}_STARIndexed"

# Genome and organism full names
genome=$(if [[ $org == "hg38" ]]; then echo "GRCh38"; elif [[ $org == "mm10" ]]; then echo "GRCm39"; else echo "Error: Unsupported organism '$org'" >&2; exit 1; fi)
org_fname=$(if [[ $org == "hg38" ]]; then echo "Homo_sapiens"; elif [[ $org == "mm10" ]]; then echo "Mus_musculus"; fi)
org_cname=$(if [[ $org == "hg38" ]]; then echo "human"; elif [[ $org == "mm10" ]]; then echo "mouse"; fi)

echo "Making reference genome for $org with Genome $genome-release-$genome_release and GENCODE-v$gencode_release"

fasta="${org_fname}.${genome}.dna.primary_assembly.fa"
gtf="gencode.v${gencode_release}.primary_assembly.annotation.gtf"

fasta_url="https://ftp.ensembl.org/pub/release-${genome_release}/fasta/$(echo ${org_fname,})/dna/${fasta}.gz"
gtf_url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_${org_cname}/release_${gencode_release}/${gtf}.gz"

fasta_in="${genome_path}/${org}-${genome_release}/${fasta}"
gtf_in="${genome_path}/${org}-${genome_release}/${gtf}"

# Download the source files if not present
if [ ! -f "$fasta_in" ]; then
    curl -fsSL "$fasta_url" | zcat > "$fasta_in"
fi
if [ ! -f "$gtf_in" ]; then
    curl -fsSL "$gtf_url" | zcat > "$gtf_in"
fi


# Modify sequence headers in the Ensembl FASTA to match the file
# "GRCh38.primary_assembly.genome.fa" from GENCODE. Unplaced and unlocalized
# sequences such as "KI270728.1" have the same names in both versions.
#
# Input FASTA:
#   >1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF
#
# Output FASTA:
#   >chr1 1
fasta_modified="${genome_path}/${org}-${genome_release}/${fasta}.modified"
# sed commands:
# 1. Replace metadata after space with original contig name, as in GENCODE
# 2. Add "chr" to names of autosomes and sex chromosomes
# 3. Handle the mitochrondrial chromosome
cat "$fasta_in" \
    | sed -E 's/^>(\S+).*/>\1 \1/' \
    | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' \
    | sed -E 's/^>MT />chrM /' \
    > "$fasta_modified"


# Remove version suffix from transcript, gene, and exon IDs in order to match
# previous Cell Ranger reference packages
#
# Input GTF:
#     ... gene_id "ENSG00000223972.5"; ...
# Output GTF:
#     ... gene_id "ENSG00000223972"; gene_version "5"; ...
gtf_modified="${genome_path}/${org}-${genome_release}/${gtf}.modified"
# Pattern matches Ensembl gene, transcript, and exon IDs for human or mouse:
ID="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"
cat "$gtf_in" \
    | sed -E 's/gene_id "'"$ID"'";/gene_id "\1"; gene_version "\3";/' \
    | sed -E 's/transcript_id "'"$ID"'";/transcript_id "\1"; transcript_version "\3";/' \
    | sed -E 's/exon_id "'"$ID"'";/exon_id "\1"; exon_version "\3";/' \
    > "$gtf_modified"


# Define string patterns for GTF tags
# NOTES:
# Since Ensembl 110, polymorphic pseudogenes are now just protein_coding.
# Readthrough genes are annotated with the readthrough_transcript tag.
BIOTYPE_PATTERN=\
"(protein_coding|protein_coding_LoF|lncRNA|\
IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|\
IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|\
TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|\
TR_V_pseudogene|TR_J_pseudogene)"
GENE_PATTERN="gene_type \"${BIOTYPE_PATTERN}\""
TX_PATTERN="transcript_type \"${BIOTYPE_PATTERN}\""
READTHROUGH_PATTERN="tag \"readthrough_transcript\""


# Construct the gene ID allowlist. We filter the list of all transcripts
# based on these criteria:
#   - allowable gene_type (biotype)
#   - allowable transcript_type (biotype)
#   - no "readthrough_transcript" tag
# We then collect the list of gene IDs that have at least one associated
# transcript passing the filters.
cat "$gtf_modified" \
    | awk '$3 == "transcript"' \
    | grep -E "$GENE_PATTERN" \
    | grep -E "$TX_PATTERN" \
    | grep -Ev "$READTHROUGH_PATTERN" \
    | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
    | sort \
    | uniq \
    > "${genome_path}/${org}-${genome_release}/gene_allowlist"


# NOTES:
# Since Ensembl 110, the PAR locus genes are included on chrY as copies of chrX
# Using the GRCh38.p13 assembly hard masks these regions on chrY, but removing the
# chrY PAR genes is still desirable so they do not end up as extra entries in the output.
# The awk command below excludes all PAR_Y genes, including XGY2.
# The non-coding gene XGY2 straddles the PAR1 boundary on chrY, and is homologous to XG on chrX.
# GRCh38-2024-A excludes XGY2, but includes SRY and ENSG00000286130, which are in an intron of XGY2,
# and RPS4Y1, which overlaps XGY2.

# Filter the GTF file based on the gene allowlist
gtf_filtered="${genome_path}/${org}-${genome_release}/${gtf}.filtered"
# Copy header lines beginning with "#"
grep -E "^#" "$gtf_modified" > "$gtf_filtered"
# Filter to the gene allowlist, and then remove PAR_Y genes
grep -Ff "${genome_path}/${org}-${genome_release}/gene_allowlist" "$gtf_modified" \
    | awk -F "\t" '$1 != "chrY" || $1 == "chrY" && $4 >= 2752083 && $4 < 56887903 && !/ENSG00000290840/' \
    >> "$gtf_filtered"


# # Create reference package
# # cellranger mkref --ref-version="$version" \
# #     --genome="$genome" --fasta="$fasta_modified" --genes="$gtf_filtered" \
# #     --nthreads=16

STAR --runThreadN ${thr} \
    --runMode genomeGenerate \
    --genomeDir ${genome_path}/${org}_STARIndexed \
    --genomeFastaFiles ${fasta_modified} \
    --sjdbGTFfile ${gtf_filtered} \
    --sjdbOverhang 99 

rm -f "${gtf_modified}" "${fasta_modified}" "${gtf_filtered}" "${genome_path}/${org}-${genome_release}/gene_allowlist"
echo "Reference genome for $org with Genome $genome-release-$genome_release and GENCODE-v$gencode_release has been created"
exit 0
