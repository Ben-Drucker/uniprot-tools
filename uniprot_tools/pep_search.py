import collections
import io
import json
import multiprocessing
import os
import pandas as pd
import re
import shutil
import subprocess
import tqdm
from . import this_dir
from .process_uniprot_output import process_uniprot_output
from multiprocessing.pool import ThreadPool
from shlex import split
from termcolor import colored
from textwrap import wrap


def _system_call(command: str):
    cmds = split(command)
    if not re.search("java -jar.*PeptideMatchCMD", command):
        raise ValueError(f"Invalid command: {command}")
    with open(os.devnull, "w") as devnull:
        proc = subprocess.run(cmds, stderr=subprocess.PIPE, stdout=devnull, shell=False)
    error_text = proc.stderr.decode("utf-8")
    end_code = proc.returncode
    if end_code != 0:
        err_block = "\n".join(
            wrap(
                f"« {error_text.strip()} »",
                width=80,
                initial_indent=" " * 4,
                subsequent_indent=" " * 8,
            )
        )
        command = f"` {command} `"
        raise ValueError(
            f"\nThe following command failed with code {end_code}:\n"
            f"    {colored(command, 'blue')}\n"
            "The command gave the error message:\n"
            f"{colored(err_block, 'red', attrs=['bold'])}"
        )


def create_haystacks(
    fasta_files: list[str],
    haystack_dir: str,
    num_lines: list[int] | None = None,
    seqs_per_chunk: int = 1_000_000,
    haystack_str: str = "",
) -> None:
    """Break a very large fasta file into more managerable chunks for the command line search tool

    Parameters
    ----------
    ``fasta_files`` :
        All the fasta files to be broken up
    ``haystack_dir`` :
        The directory where the "haystacks" will be stored
    ``num_lines`` :
        The number of lines in each original, large fasta file
    ``seqs_per_chunk`` :
        The number of sequences per output haystack
    ``haystack_str`` :
        Additional string to be added to the haystack file name
    """

    if num_lines is None:
        num_lines = [0] * len(fasta_files)
    for fasta_file, num_lines_i in zip(fasta_files, num_lines):
        with open(fasta_file) as fp:
            e = 0
            new = 0
            haystack_str = os.path.basename(fasta_file)
            lambda_open_haystack = lambda: open(
                f"{haystack_dir}/haystack-{int(new // seqs_per_chunk)}{'-' if haystack_str else ''}{haystack_str.split(".")[0]}.fasta",
                "w",
            )
            fp_big = lambda_open_haystack()
            for line in tqdm.tqdm(
                fp, total=num_lines_i - 1, desc="Scanning through FASTA Progress"
            ):
                e += 1
                if ">" in line:
                    if new % seqs_per_chunk == 0 and new != 0:
                        fp_big.close()
                        fp_big = lambda_open_haystack()
                    # reset
                    new += 1
                fp_big.write(line)
            fp_big.close()


def create_index(
    fasta_files: list[str],
    index_dir: str,
    num_search_procs: int = multiprocessing.cpu_count(),
) -> None:
    """Create indexes that can be used with the command line search tool

    Parameters
    ----------
    ``fasta_files`` :
        Full paths to all the fasta files whose proteins are to be indexed
    ``index_dir`` :
        The directory where the indexes will be created
    """
    output_indexes_dir_names = [
        re.sub(r"haystacks(.*)\.fasta", r"haystack-index\1", h) for h in fasta_files
    ]

    commands = sorted(
        [
            f"java -jar '{this_dir.this_dir('')}bin/PeptideMatchCMD_1.1.jar' -f -a index -d '{h}'"
            f" -i '{index_dir}/{os.path.basename(i)}'"
            for h, i in zip(fasta_files, output_indexes_dir_names)
        ]
    )

    with ThreadPool(processes=num_search_procs) as pool:
        pool.map(_system_call, commands)


class CustomTQDM(tqdm.tqdm):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.desc = "Post Processing Progress"
        self.smoothing = 0.1


def post_process_java_output(
    max_concurrency: int = 4, files: list[str] = []
) -> dict[str, set[str]]:
    """Based on outputs to java command line peptide search tool, produce a mapping from peptides \
        to sets of proteins

    Parameters
    ----------
    ``max_concurrency`` :
        The number of threads to use when searching. Each thread will handle one peptide.
    ``files`` :
        The files to post process.

    Returns
    -------
        Mapping from each peptide found in the input files to sets of proteins
    """
    args = [{"pep_search_file": x, "tid": i} for i, x in enumerate(files)]
    results = []
    with ThreadPool(processes=max_concurrency) as pool:
        with tqdm.tqdm(desc="Post Processing Progress", smoothing=0.1, total=len(args)) as pb:
            for res in pool.imap_unordered(_post_process_java_output_worker, args):
                results.append(res)
                pb.update(1)

    updated_together = collections.defaultdict(set)
    for result in results:
        for r in result:
            updated_together[r].update(result[r])
    return updated_together


def _post_process_java_output_worker(kwargs) -> dict[str, set[str]]:
    pep_search_file = kwargs["pep_search_file"]
    tid = kwargs["tid"]
    with open(pep_search_file) as fp:
        lines = fp.readlines()
    lines = lines[1:]
    lines[0] = lines[0][2:-1] + "\tX\n"
    with open(f"temp-{tid}.tsv", "wt") as fp:
        fp.write("".join(lines))
    tab = pd.read_csv(f"temp-{tid}.tsv", sep="\t")
    os.unlink(f"temp-{tid}.tsv")
    tab["Uniprot"] = [
        pd.NA if x == "No match" else re.sub(".*\\|(.*)\\|.*", r"\1", x) for x in tab["Subject"]
    ]
    tab = tab[["Query", "Uniprot"]]

    seq_to_uniprots = collections.defaultdict(set)
    for seq, uniprot in zip(tab["Query"], tab["Uniprot"]):
        if not pd.isna(uniprot):
            seq_to_uniprots[seq].add(uniprot)
        else:
            seq_to_uniprots[seq]

    return seq_to_uniprots


def pep_search(
    peptides: list[str],
    index_dir: str,
    num_search_procs: int = multiprocessing.cpu_count(),
    intermediate_dir: str = "pep_searches",
    delete_index_when_done: bool = False,
) -> None:
    """Searches for matching proteins form given peptides in a given index

    Parameters
    ----------
    ``peptides`` :
        The peptides to search for
    ``index_dir`` :
        The directory whose index to search with
    ``num_search_procs`` :
        The number of processes to use when searching. Each processess will handle one peptide.

    Returns
    -------
        Mapping peptide to list of matching protein accession IDs
    """
    indexes = [
        f"{index_dir}/{x}"
        for x in sorted(os.listdir(index_dir))
        if not re.search(r"^\.", x) and not re.search(r"tmp", x)
    ]
    if (
        os.path.exists(intermediate_dir)
        and len([x for x in os.listdir(intermediate_dir) if re.search(r"^\.", x)]) > 0
    ):
        raise FileExistsError(
            "pep_searches directory already exists containing content. Delete it and try again."
        )
    output_results_names = [
        f"{os.path.join(intermediate_dir, i.split("/")[-1])}.pepsearchres" for i in indexes
    ]

    with open(inputs := os.path.join(intermediate_dir, "pepfile.txt"), "w") as fp:
        fp.write("\n".join(peptides))

    commands = sorted(
        [
            f"java -jar '{this_dir.this_dir('')}bin/PeptideMatchCMD_1.1.jar' -f -a query -Q"
            f" {inputs} -i '{i}' -o {o} -l"
            for i, o in zip(indexes, output_results_names)
        ]
    )
    assert len(commands) > 0, f"No indexes found in {index_dir}"
    with ThreadPool(processes=num_search_procs) as pool:
        pool.map(_system_call, commands)

    if delete_index_when_done:
        for i in indexes:
            shutil.rmtree(i)


existence_codes = {
    "evidence at protein level": 1,  # best
    "evidence at transcript level": 2,
    "inferred from homology": 3,
    "predicted": 4,
    "uncertain": 5,  # worst
}


def extract_pep_to_info(combined_df: pd.DataFrame, peptide_to_accessions: dict[str, set[str]]):
    """TODO

    Parameters
    ----------
    ``combined_df`` :
        TODO
    ``peptide_to_accessions`` :
        TODO

    Returns
    -------
        TODO
    """
    acc_to_is_from_vertebrate = {
        k: ("Vertebrata" in v)
        for k, v in tqdm.tqdm(
            list(combined_df.set_index("Entry").to_dict()["Taxonomic lineage"].items()),
            desc='Mapping Accessions to "Has Vertebrate" Progress',
        )
    }
    acc_to_best_prot_evidence = {
        k: existence_codes[v.lower()]
        for k, v in combined_df.set_index("Entry").to_dict()["Protein existence"].items()
    }
    pep_to_has_vertebrate = {}
    pep_to_best_prot_evidence = {}
    key_errors = set()
    for p, v in tqdm.tqdm(
        peptide_to_accessions.items(), desc='Mapping Peptides to "Has Vertebrate" Progress'
    ):
        pep_to_has_vertebrate[p] = False
        if len(v) > 0:
            best_evidence = 6
            for vi in v:
                vi = re.sub(r"-\d+$", "", vi)
                try:
                    if acc_to_is_from_vertebrate[vi]:
                        pep_to_has_vertebrate[p] = True
                    best_evidence = min(best_evidence, acc_to_best_prot_evidence[vi])
                    pep_to_best_prot_evidence[p] = best_evidence
                except KeyError as ke:
                    key_errors.add(ke)
        else:
            pep_to_has_vertebrate[p] = None
            pep_to_best_prot_evidence[p] = None

    return (
        pep_to_has_vertebrate,
        pep_to_best_prot_evidence,
        {str(ke).replace("'", "") for ke in key_errors},
        acc_to_is_from_vertebrate,
    )


if __name__ == "__main__":
    from .get_info import accession_to_prot_info

    jsons = accession_to_prot_info(["P00519", "P06493"], format_="json", compressed_download=True)

    def jsons_to_pdb_dfs(jsons):
        read_in = process_uniprot_output(jsons)
        all_pdb_info = {
            "Uniprot ID": [],
            "PDB ID": [],
            "Method": [],
            "Resolution A": [],
            "Chain IDs": [],
            "Start Res": [],
            "End Res": [],
        }
        assert isinstance(read_in, dict)
        for results_chunk in read_in.values():
            for result in results_chunk:
                assert isinstance(result, dict)
                for pdb in [
                    d for d in result["uniProtKBCrossReferences"] if d["database"] == "PDB"
                ]:
                    all_pdb_info["PDB ID"].append(pdb["id"])
                    all_pdb_info["Uniprot ID"].append(result["primaryAccession"])
                    mini_df = pd.DataFrame.from_records(pdb["properties"])
                    mini_df = mini_df.set_index("key").T
                    res_a = mini_df.at["value", "Resolution"]
                    numberic_part = re.search(r"^(\d+\.\d+)\s*.*$", res_a)
                    if numberic_part is not None:
                        all_pdb_info["Resolution A"].append(float(numberic_part.group(1)))
                    else:
                        all_pdb_info["Resolution A"].append(None)
                    all_pdb_info["Method"].append(mini_df.at["value", "Method"])
                    id_extract = re.search(
                        r"^([A-Z0-9/]+)=(\d+)-(\d+)", mini_df.at["value", "Chains"]
                    )
                    if id_extract is not None:
                        all_pdb_info["Chain IDs"].append(id_extract.group(1).split("/"))
                        all_pdb_info["Start Res"].append(int(id_extract.group(2)))
                        all_pdb_info["End Res"].append(int(id_extract.group(3)))
                    else:
                        all_pdb_info["Chain IDs"].append(None)
                        all_pdb_info["Start Res"].append(None)
                        all_pdb_info["End Res"].append(None)

        all_pdb_df = pd.DataFrame(all_pdb_info)
        return all_pdb_df

    def get_best_res_for_range(pdb_df: pd.DataFrame, start: int, end: int):
        subset = pdb_df[(pdb_df["Start Res"] <= start) & (pdb_df["End Res"] >= end)]
        if len(subset) > 0:
            return subset.sort_values("Resolution A").iloc[0]["PDB ID"]
        else:
            raise ValueError(f"Could not find PDB entry for residue range {start}-{end}.")
