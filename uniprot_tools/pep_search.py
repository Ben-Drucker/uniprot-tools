import asyncio, collections, io, multiprocessing, os, re, shutil, ssl, subprocess
from multiprocessing.pool import ThreadPool
from shlex import quote, split

import aiohttp
import certifi
import pandas as pd
import tqdm
from requests import HTTPError
from tqdm.asyncio import tqdm_asyncio


def _system_call(command):
    cmds_sanitized = quote(command)
    cmds = split(cmds_sanitized)
    if not re.search("java -jar.*PeptideMatchCMD", command):
        raise ValueError(f"Invalid command: {command}")
    end_code = subprocess.run(cmds).returncode
    if end_code != 0:
        raise ValueError(f"Command failed with code {end_code}: {command}")


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
            f"java -jar {__file__}/../PeptideMatchCMD_1.1.jar -a index -d '{h}' -i"
            f" '{index_dir}/{os.path.basename(i)}'"
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
    print(indexes)
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
            f"java -jar {__file__}/../PeptideMatchCMD_1.1.jar -a query -Q {inputs} -i"
            f" '{i}' -o {o} -l"
            for i, o in zip(indexes, output_results_names)
        ]
    )
    assert len(commands) > 0, f"No indexes found in {index_dir}"
    with ThreadPool(processes=num_search_procs) as pool:
        pool.map(_system_call, commands)

    if delete_index_when_done:
        for i in indexes:
            shutil.rmtree(i)


async def tqdm_gather_with_concurrency(*coroutines, max_concurrency=4, **tqdm_kwargs):
    semaphore = asyncio.Semaphore(max_concurrency)

    async def semaphore_coroutine(coroutine):
        async with semaphore:
            return await coroutine

    return await tqdm_asyncio.gather(*(semaphore_coroutine(c) for c in coroutines), **tqdm_kwargs)


async def get_verts(peptide_to_accessions: dict[str, set[str]], on_chunk=0, out_of=0):
    CHUNK_SIZE = 500
    sslcontext = ssl.create_default_context(cafile=certifi.where())

    accessions = set()
    for v in peptide_to_accessions.values():
        accessions.update(v)

    accessions = sorted(accessions)

    urls = []
    fields = [
        "accession",
        "reviewed",
        "protein_name",
        "organism_name",
        "cc_alternative_products",
        "protein_existence",
        "lineage",
        "virus_hosts",
    ]

    for chunk in tqdm.tqdm(
        [accessions[i : i + CHUNK_SIZE] for i in range(0, len(accessions), CHUNK_SIZE)],
        desc="Generating URLs",
    ):
        chunk = [x for x in chunk if str(x) != "nan"]
        query_str = "query=" + "(accession:" + f")+OR+(accession:".join(chunk) + ")"
        fields_str = "fields=" + ",".join(fields)
        fmt = "format=tsv"
        url = f"https://rest.uniprot.org/uniprotkb/search?{fields_str}&{query_str}&{fmt}&size={len(chunk)}"
        urls.append(url)

    async def main(urls):
        async with aiohttp.ClientSession() as session:
            tasks = []
            for i, url in enumerate(urls):
                task = fetch(url, session, i)
                tasks.append(task)

            responses = await tqdm_asyncio.gather(
                *tasks, desc=f"Uniprot Query Progress Chunk [{on_chunk}/{out_of}]"
            )
            return responses

    return await main(urls)


existence_codes = {
    "evidence at protein level": 1,  # best
    "evidence at transcript level": 2,
    "inferred from homology": 3,
    "predicted": 4,
    "uncertain": 5,  # worst
}


def process_tsvs(string_tsvs: list[str]) -> pd.DataFrame:
    """Turn Uniprot text tsvs into one, combined pandas `DataFrame`

    Parameters
    ----------
    ``string_tsvs`` :
        List of tsv strings. Note: each tsv must have the same columns.

    Returns
    -------
        Table of combined results.
    """
    pandaized = []
    for tsv in tqdm.tqdm(string_tsvs, desc="Turning Uniprot Results into Pandas DataFrames"):
        pandaized.append(pd.read_csv(io.StringIO(tsv), sep="\t"))

    print("Concatenating Uniprot TSVs...")
    combined_df = pd.concat(pandaized, ignore_index=True)

    return combined_df


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
