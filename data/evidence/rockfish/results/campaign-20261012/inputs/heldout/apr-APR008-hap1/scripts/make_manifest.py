#!/usr/bin/env python3
"""Write manifest.json for the APR008 hap1 held-out query.

Usage: python3 -I make_manifest.py DIR RAW_SHA256
DIR holds raw/, reads.fa, offsets.u64 and conversion_stats.json.
"""

import datetime
import json
import os
import struct
import sys

d, raw_sha = sys.argv[1], sys.argv[2]
st = json.load(open(os.path.join(d, "conversion_stats.json")))

off_path = os.path.join(d, "offsets.u64")
n_off = os.path.getsize(off_path) // 8
with open(off_path, "rb") as f:
    first = struct.unpack("<Q", f.read(8))[0]
    f.seek((n_off - 1) * 8)
    last = struct.unpack("<Q", f.read(8))[0]
assert n_off == st["output_records"] + 1 and first == 0 and last == st["output_bases"]

ftp = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/050/494/745/GCA_050494745.1_apr008.1_v1"
raw_name = "GCA_050494745.1_apr008.1_v1_genomic.fna.gz"
m = {
    "name": "apr-APR008-hap1",
    "description": "Arab Pangenome Reference (APR / UPR) sample APR008, haplotype 1 (Hi-C phased hifiasm), whole-genome held-out query against HPRC release 2 (Human472).",
    "source_url": f"{ftp}/{raw_name}",
    "source_dir_url": ftp + "/",
    "accession": "GCA_050494745.1",
    "accession_version": 1,
    "assembly_name": "apr008.1 v1",
    "linked_haplotype_2": "GCA_050494665.1",
    "bioproject": "PRJNA1152077 (haplotype 1; parent PRJNA1166278; reads PRJNA1108179)",
    "biosample": "SAMN41227928",
    "wgs_project": "JBHMUJ01",
    "sample_id": "pbAPR008 (GenBank isolate); APR008 (OpenHGL; assembly 210019_APR008.hap1 in human579 v4.0)",
    "sample_sex": "female",
    "sample_geo_loc": "United Arab Emirates",
    "submitter": "Mohammed Bin Rashid University of Medicine and Health Sciences (MBRU)",
    "assembly_method": "hifiasm v0.19.8 (PacBio HiFi + ONT ultralong + Hi-C)",
    "genbank_release_date": "2025-05-20",
    "publication": "Nassir et al., A draft UAE-based Arab pangenome reference, Nat Commun 16:6747 (2025), doi:10.1038/s41467-025-61645-w, PMC12290100",
    "openhgl": "OpenHGL v4.0 (human579), dataset APR v1, https://lh3.github.io/OpenHGL/ ; s3://openhgl/human/human579/human579.meta.tsv lists 210019_APR008.hap1; not in human579.misasm.tsv",
    "license": "Deposited in NCBI GenBank (INSDC), no use restrictions stated by NCBI; the paper's data availability statement says the assemblies are publicly available. The article text is CC BY-NC-ND 4.0. OpenHGL redistributes the APR assemblies under CC0 (AWS Open Data registry).",
    "publisher_checksum": {"algorithm": "md5", "file": raw_name, "value": "9217bfddb3cd4b28e0382809e8116c76", "source": f"{ftp}/md5checksums.txt", "verified": True},
    "publisher_uncompressed_checksum": {"algorithm": "md5", "file": "GCA_050494745.1_apr008.1_v1_genomic.fna", "value": "ca7975f7f0a38702fd47be559a3fadf6", "source": f"{ftp}/uncompressed_checksums.txt"},
    "raw_file": f"raw/{raw_name}",
    "raw_sha256": raw_sha,
    "reads_fa_sha256": st["reads_fa_sha256"],
    "before": {
        "records": st["input_records"],
        "total_bases": st["input_bases"],
        "acgt_bases": st["input_acgt_bases"],
        "lowercase_bases": st["input_lowercase_bases"],
        "non_acgt_by_char": st["non_acgt_by_char"],
        "record_types": "79 unplaced scaffolds (GenBank assembly_report: all unplaced-scaffold, no chromosome assignment, no mitochondrial record)",
        "ncbi_stats": {"contigs": 81, "scaffolds": 79, "contig_n50": 140723744, "scaffold_n50": 140723744, "total_sequence_length": 3038047338, "total_ungapped_length": 3038047138},
    },
    "after": {
        "records": st["output_records"],
        "total_bases": st["output_bases"],
        "longest_record": st["longest_output_record"],
        "shortest_record": st["shortest_output_record"],
        "input_records_with_non_acgt": st["input_records_with_non_acgt"],
        "input_records_fully_dropped": st["input_records_fully_dropped"],
    },
    "non_acgt_bases_dropped": st["non_acgt_bases_dropped"],
    "conversion_rule": "Keep every record of the haplotype. Uppercase all bases. Split each record at every run of characters other than A/C/G/T, dropping those characters, and drop resulting pieces shorter than 1 base. Name each piece <origname>_<start>_<end>, where origname is the first whitespace-delimited word of the input header (GenBank accession.version) and start/end are 0-based half-open coordinates within the input record. Output is two-line FASTA (header line, then the whole sequence on one line).",
    "conversion_script": "scripts/acgt_split.py",
    "offsets_u64": "offsets.u64: little-endian uint64 cumulative record start offsets, n_records+1 values, first 0, last = total bases",
    "offsets_count": n_off,
    "date": datetime.date.today().isoformat(),
    "heldout_evidence": {
        "summary": "APR008 is a UAE Arab individual sequenced by MBRU/Dubai Health for the Arab Pangenome Reference; it is a different project and sample from HPRC.",
        "hprc2_assembly_index": "https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/2a7e3a5b1ca10454762699176b1a0a1027ccbf21/data_tables/assemblies_release2_v1.0.index.csv (466 assemblies, 234 sample IDs incl. CHM13 and GRCh38): no sample_id or genbank_accession matching APR008, pbAPR008, APR, or GCA_050494745",
        "hprc2_sample_metadata": "data_tables/sample/hprc_release2_sample_metadata.csv (234 samples): no APR/pbAPR sample, BioSample SAMN41227928 absent; populations are 1000 Genomes/other codes, none UAE or Arab",
        "openhgl": "OpenHGL lists APR (v1) and HPRC (r2-v1.0.1) as separate source datasets; APR008 appears only under APR",
    },
}
with open(os.path.join(d, "manifest.json"), "w") as g:
    json.dump(m, g, indent=1)
    g.write("\n")
print(json.dumps(m, indent=1))
