from contextlib import redirect_stderr
import traceback


def main(input, out_fn):

    lines = [
        "# id: 'contamination'",
        "# plot_type: 'table'",
        "# section_name: 'Contaminate Check'",
        "db\tpotential_hits",
    ]

    for key, path in input.items():
        n_hits = open(path).read().count("\n")
        lines.append(f"{key}\t{n_hits}")

    with open(out_fn, "w") as out_f:
        out_f.write("\n".join(lines))

    return

if __name__ == "__main__":
    with (
            open(snakemake.log.e, "w") as ef,
            redirect_stderr(ef)
    ):
        try:
            main(
                snakemake.input,
                snakemake.output["report"],
            )
        except Exception as e:
            traceback.print_exc(file=ef)
