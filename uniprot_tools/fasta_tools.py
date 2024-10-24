"""``.fasta`` file reading and writing."""

import os
from typing import Literal

import pandas as pd
import tqdm


def read_fasta(infile: str, sort_keys: bool = True, show_progress: bool = False) -> dict[str, str]:
    """Read a fasta file and return a mapping from description to sequence

    Parameters
    ----------
    ``infile`` :
        The path to the fasta file to read
    ``sort_keys`` :
        Whether to sort the keys in the output mapping

    Returns
    -------
        The output mapping from description to sequence
    """
    mapping = {}
    running = ""
    line = ""
    with open(infile) as fp:
        running_description = fp.readline().strip().replace(">", "")
        if show_progress:
            prog_bar = tqdm.tqdm(desc="Bytes Read", total=os.path.getsize(infile))
            prog_bar.update(len(running_description) + 2)
        else:
            prog_bar = None
        for line in fp:
            if ">" in line:
                mapping[running_description] = running
                running = ""
                running_description = line.strip().replace(">", "")
            else:
                running += line.strip()
            if prog_bar:
                prog_bar.update(len(line))
        # if prog_bar:
        #     prog_bar.update(len(line))
    mapping[running_description] = running
    if sort_keys:
        mapping = {k: v for k, v in sorted(mapping.items(), key=lambda item: item[0])}
    return mapping


def write_fasta(
    seqs: list[str],
    seq_descriptions: list[str] | None,
    outfile: str = "output.fasta",
    do_sort: bool = True,
    only_unique: bool = True,
    sort_by_desc_or_seq: Literal["desc", "seq"] = "desc",
) -> None:
    """Write a list of sequences to a fasta file, optionally with descriptions

    Parameters
    ----------
    ``seqs`` :
        List of sequences to write to output fasta file.
    ``seq_descriptions`` :
        List of descriptions for each sequence, in the same order as ``seqs``.
    ``outfile`` :
        The path to the output fasta file.
    ``do_sort`` :
        Whether to sort the sequences either by description or sequence.
    ``only_unique`` :
        Whether to only write unique sequences.
    ``sort_by_desc_or_seq`` :
        If ``do_sort`` is True, whether to sort by description or sequence. Otherwise, ignored.
    """
    if seq_descriptions is None:
        seq_descriptions = [f"Seq-{i}" for i in range(len(seqs))]

    df = pd.DataFrame({"Seq": seqs, "Description": seq_descriptions})
    if only_unique:
        df = df.drop_duplicates(subset=["Seq"], keep="first")
    if do_sort:
        match sort_by_desc_or_seq:
            case "desc":
                df = df.sort_values(by="Description")
            case "seq":
                df = df.sort_values(by="Seq")

    seqs, seq_descriptions = list(df["Seq"]), list(df["Description"])
    with open(outfile, "w") as fp:
        for i, seq in enumerate(seqs):
            if i + 1 == len(seqs):
                newline = ""
            else:
                newline = "\n"
            fp.write(f">{seq_descriptions[i]}\n{seq}{newline}")
