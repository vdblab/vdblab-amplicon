import glob
import gzip
import os
import pandas as pd
import subprocess
import yaml

from snakemake.utils import min_version, validate


def get_pipeline_version():
    return (
        subprocess.check_output(["git", "describe", "--always"], cwd=workflow.basedir)
        .decode("utf-8")
        .strip()
    )


def make_shard_names(nshards):
    return [f"{x:03}" for x in range(1, config["nshards"] + 1)]


def load_manifest(manifest_path, manifest_schema_path):
    manifest = pd.read_csv(manifest_path, sep="\t")
    validate(manifest, manifest_schema_path)
    return manifest.set_index("sample_id")


def get_samples_from_manifest(manifest):
    return [s for s in manifest.index.tolist() if s != "orphans"]

def extract_primers(oligos_path):
    with open(config["oligos"], "r") as in_f:
        for line in in_f:
            assert line.upper().startswith("PRIMER"), "invalid oligos file, should start with `primer`"
            primers = line.strip().split("\t")
            primers.pop(0)
            return primers

def extract_barcodes(oligos_path):
    sample_ids = []
    barcodes_F = []
    barcodes_R = []

    with open(config["oligos"], "r") as in_f:
        for line in in_f:
            # this apparently can be lowercase in some instances
            if not line.upper().startswith("BARCODE"):
                continue

            pieces = line.split()

            sample_id = pieces[3].split("..")[0]
            # borrowed this regex from the old barcode cleaning script
            sample_id = re.sub(r"[_%+;: -]+", ".", sample_id)

            sample_ids.append(sample_id)
            barcodes_F.append(pieces[1])
            barcodes_R.append(pieces[2])

    assert len(sample_ids) > 0, "No samples found in oligos file"

    barcodes = pd.DataFrame({"sample_id": sample_ids, "F": barcodes_F, "R": barcodes_R})
    return barcodes[["sample_id", "F", "R"]]


# https://stackoverflow.com/questions/37874936/how-to-check-empty-gzip-file-in-python
def gz_size(fname):
    with gzip.open(fname, "rb") as f:
        return f.seek(0, whence=2)


def write_manifest_and_missing(sample_ids, fastq_template, manifest_path, missing_path):
    R1 = []
    R2 = []
    for s in sample_ids:
        R1_path = fastq_template.format(sample=s, dir=1)
        if os.path.exists(R1_path) and gz_size(R1_path) > 0:
            R1.append(R1_path)
        else:
            R1.append("")

        R2_path = fastq_template.format(sample=s, dir=2)
        if os.path.exists(R2_path) and gz_size(R2_path) > 0:
            R2.append(R2_path)
        else:
            R2.append("")

    manifest = pd.DataFrame({"sample_id": sample_ids, "R1": R1, "R2": R2})
    manifest = manifest[["sample_id", "R1", "R2"]]

    is_incomplete = (manifest["R1"] == "") | (manifest["R2"] == "")

    manifest[is_incomplete].to_csv(missing_path, sep="\t", index=False)

    manifest = manifest[~is_incomplete]
    manifest.to_csv(manifest_path, sep="\t", index=False)
