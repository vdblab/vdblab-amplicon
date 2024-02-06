#! /usr/local/env python3
import os
import pandas as pd
import re
import sys
import traceback

from io import StringIO
from zipfile import ZipFile


def extract_contam(fastqc_zip):
    outlines = []

    target = os.path.basename(fastqc_zip).replace(".zip", "")
    archive = ZipFile(fastqc_zip, "r")
    with archive.open(f"{target}/fastqc_data.txt") as content:

        in_section_of_interest = False
        for line in content:
            line = line.decode().strip()

            # Only worry about the overrepresented sequences section
            if "Overrepresented sequences" in line:
                in_section_of_interest = True
            elif ">>END_MODULE" in line and in_section_of_interest:
                break
            elif in_section_of_interest:
                outlines.append(re.sub(r"^#", "", line))

    contam = pd.read_csv(StringIO("\n".join(outlines)), sep="\t")

    # remove the "No Hit"s: for amplicon sequencing these are
    # the actual 16S sequences that are enriched, not some sequencing artifact
    contam = contam[contam["Possible Source"] != "No Hit"]

    return contam


def main(inputs, sample_id, out_f):
    zip_inputs = {"R1": inputs[0]}
    if len(inputs) == 2:
        zip_inputs["R2"] = inputs[1]
    dfs = []
    for k, v in zip_inputs.items():
        contam = extract_contam(v)
        contam["Read Direction"] = k
        dfs.append(contam)

    contam = pd.concat(dfs)
    contam["SampleID"] = sample_id

    contam.to_csv(out_f, sep="\t", index=False)
    return


if __name__ == "__main__":
    with (
        open(snakemake.log.e, "w") as ef,
        open(snakemake.log.o, "w") as of
    ):
        sys.stderr = ef
        sys.stdout = of

        try:
            main(
                snakemake.input,
                snakemake.wildcards.sample,
                snakemake.output.o,
            )
        except Exception as e:
            traceback.print_exc(file=ef)
