import multiprocessing, os, random, re, time, warnings
from concurrent import futures
from concurrent.futures import ThreadPoolExecutor
from multiprocessing.pool import IMapIterator, ThreadPool, starmapstar
from typing import Callable, Iterable, Literal

import numpy as np
import requests
import tqdm
from pep_search import process_tsvs

# from requests.adapters import HTTPAdapter

# s = requests.Session()
# s.mount('https://rest.uniprot.org', HTTPAdapter(max_retries=5000))


def request_worker(url: str, request_kwargs: dict = {}) -> str | requests.HTTPError:
    if "timeout" not in request_kwargs:
        request_kwargs["timeout"] = 5

    def inner_main():
        with requests.get(url, **request_kwargs) as r:
            if r.status_code == 200:
                return r.text
            else:
                msg = f"Got status code {r.status_code} from {url}. Response: {r.text}"
                warnings.warn(msg)
                return requests.HTTPError(msg)

    while True:
        try:
            return inner_main()
        except requests.exceptions.ReadTimeout:
            time.sleep(10 * random.random())
            continue
        except requests.exceptions.ConnectionError:
            time.sleep(10 * random.random())
            continue


def parallel_requests(
    urls: list[str],
    num_concurrent: int | Literal["num_urls", "num_cores"] = "num_cores",
    request_kwargs: dict = {},
    prog_bar_kwargs: dict = {},
):
    results = []
    match num_concurrent:
        case "num_urls":
            num_concurrent = len(urls)
        case "num_cores":
            num_concurrent = multiprocessing.cpu_count()
        case _:
            assert isinstance(
                num_concurrent, int
            ), "num_concurrent must be an int if not 'num_urls' or 'num_cores'"

    with ThreadPoolExecutor(num_concurrent, "request-thread") as tpe:
        futures_ = [tpe.submit(request_worker, url, request_kwargs) for url in urls]
        with tqdm.tqdm(total=len(urls), **prog_bar_kwargs) as prog:
            try:
                for f in futures.as_completed(futures_):
                    prog.update(1)
                    results.append(f.result())
            except KeyboardInterrupt:
                tpe.shutdown(wait=False, cancel_futures=True)
                raise KeyboardInterrupt(
                    "Keyboard Interrupt received. Cancelling requests and shutting down..."
                )
    return results


def read_fasta(infile: str, sort_keys: bool = True) -> dict[str, str]:
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
    with open(infile) as fp:
        running_description = fp.readline().strip().replace(">", "")
        prog_bar = tqdm.tqdm(desc="Bytes Read", total=os.path.getsize(infile))
        for line in fp:
            if ">" in line:
                mapping[running_description] = running
                running = ""
                running_description = line.strip().replace(">", "")
            else:
                running += line.strip()
            prog_bar.update(len(line))
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

    if do_sort:
        match sort_by_desc_or_seq:
            case "desc":
                argsort = np.argsort(seq_descriptions)
            case "seq":
                argsort = np.argsort(seqs)
        seqs, seq_descriptions = [seqs[i] for i in argsort], [seq_descriptions[i] for i in argsort]
    if only_unique:
        unq_tuples = list(set(zip(seqs, seq_descriptions)))
        seqs, seq_descriptions = [i[0] for i in unq_tuples], [i[1] for i in unq_tuples]
    with open(outfile, "w") as fp:
        for i, seq in enumerate(seqs):
            fp.write(f">{seq_descriptions[i]}\n{seq}\n")


from multiprocessing.pool import Pool


def istarmap(p: Pool | ThreadPool, func: Callable, iterable: Iterable, chunksize=1):
    """starmap-version of imap"""
    p._check_running()
    if chunksize < 1:
        raise ValueError("Chunksize must be 1+, not {0:n}".format(chunksize))

    task_batches = Pool._get_tasks(func, iterable, chunksize)
    result = IMapIterator(p)
    p._taskqueue.put(
        (
            p._guarded_task_generation(result._job, starmapstar, task_batches),
            result._set_length,
        )
    )
    return (item for chunk in result for item in chunk)


class TqdmParallel:
    @classmethod
    def tqdm_starmap(
        cls,
        worker_fn,
        worker_args,
        num_workers=multiprocessing.cpu_count() - 1,
        cleanup_fn=lambda: None,
        processes_or_threads: Literal["processes", "threads"] = "processes",
        tqdm_class=tqdm.tqdm,
        **tqdm_kwargs,
    ):
        def end_actions(
            pool: Pool | ThreadPool | None, partial_results: list, num_expected_results: int
        ):
            if pool is not None:
                pool.terminate()
                pool.close()
                cleanup_fn()
                print(
                    f"Returning partial results ([{len(partial_results)}/{num_expected_results}])"
                )
                return partial_results  # , worker_args[len(partial_results):]

        results = []
        pool = None
        if "smoothing" not in tqdm_kwargs:
            tqdm_kwargs["smoothing"] = 0.05
        if "leave" not in tqdm_kwargs:
            tqdm_kwargs["leave"] = True
        if "desc" not in tqdm_kwargs:
            tqdm_kwargs["desc"] = "Progress: "
        if "bar_format" not in tqdm_kwargs:
            tqdm_kwargs["bar_format"] = "{desc}{percentage:3.0f}%|{bar:25}{r_bar}"
        try:
            match processes_or_threads:
                case "processes":
                    cls_ = Pool
                case "threads":
                    cls_ = ThreadPool
            with cls_(num_workers) as pool:
                for res in tqdm_class(
                    istarmap(pool, worker_fn, worker_args),
                    **tqdm_kwargs,
                    total=len(worker_args),
                ):
                    results.append(res)
        except IndexError as e:
            if str(e) != "pop from an empty deque":
                print(e)
                return end_actions(pool, results, len(worker_args))
            else:
                return end_actions(pool, results, len(worker_args))
        except KeyboardInterrupt:
            return end_actions(pool, results, len(worker_args))
        except Exception as e:
            print(e)
            return end_actions(pool, results, len(worker_args))

        cleanup_fn()
        print("Parallel tasks successfully finished!")
        return results


"""https://rest.uniprot.org/uniparc/stream?compressed=true&fields=upi,organism,accession,organism_id,protein&format=tsv&query=%28%28uniparc%3AUPI001FEE9A09%29+OR+%28uniparc%3AUPI001FEEEFDE%29+OR+%28uniparc%3AUPI001FF6F0CB%29+OR+%28uniparc%3AUPI001FFAC6BE%29%29"""


def accession_to_prot_info(
    accessions: list[str] | list[int],
    columns: list[str],
    knowledge_base: Literal["uniprotkb", "uniparc", "taxonomy"] = "uniprotkb",
    num_simultaneous_requests: int = 100,
    get_urls_only: bool = False,
) -> list[str]:
    """

    Parameters
    ----------
    ``accessions`` :
        List of uniprotkb/uniparc IDs
    ``columns`` :
        Columns to get from API. See notes for standard columns to copy and paste.
    ``knowledge_base`` :
        The Uniprot knowledge base to use, by default "uniprotkb"

    Returns
    -------
        List of TSV strings to be converted by function `process_tsvs`

    Raises
    ------
        ``ValueError`` :
            If ``knowledge_base`` is not one of ``uniprotkb`` or ``uniparc``

    Notes
    -----
    - Standard UniprotKB columns::

            [
                "accession",
                "reviewed",
                "protein_name",
                "gene_names",
                "organism_name",
                "cc_alternative_products"
            ]


    - Standard Uniparc columns::

            [
                "upi",
                "accession",
                "organism",
                "protein"
            ]

    - Standard Taxonomy columns::

            [
                "id",
                "common_name",
                "scientific_name",
                "lineage",
            ]
    """
    CHUNK_SIZE = 500
    accessions = sorted(accessions)
    match knowledge_base:
        case "uniprotkb":
            correct_id_pat = r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
            query_specifier = "accession"
        case "uniparc":
            correct_id_pat = r"UPI[A-Z0-9]{10}"
            query_specifier = "uniparc"
        case "taxonomy":
            correct_id_pat = r"\d+"
            query_specifier = "tax_id"
        case _:  # type: ignore
            raise ValueError("`database` must be either `uniprotkb` or `uniparc`")
    urls = []
    chunks = [
        [str(x) for x in accessions[i : i + CHUNK_SIZE]]
        for i in range(0, len(accessions), CHUNK_SIZE)
    ]
    for chunk in tqdm.tqdm(chunks, desc="Generating and validating URLs"):
        chunk = [x for x in chunk if str(x) != "nan"]
        for c in chunk:
            if not re.search(correct_id_pat, c):
                raise ValueError(f"{c} is not a valid accession for database `{knowledge_base}`.")

        query_str = (
            "query=" + f"({query_specifier}:" + f")+OR+({query_specifier}:".join(chunk) + ")"
        )
        fields_str = "fields=" + ",".join(columns)
        fmt = "format=tsv"
        url = (
            f"https://rest.uniprot.org/{knowledge_base}/search?"
            f"{fields_str}&{query_str}&{fmt}&size={len(chunk)}"
        )
        urls.append(url)
    try:
        if get_urls_only:
            return urls
        else:
            res = TqdmParallel.tqdm_starmap(
                worker_fn=request_worker,
                worker_args=[[url] for url in urls],
                num_workers=num_simultaneous_requests,
                processes_or_threads="threads",
            )
    except Exception:
        assert isinstance(res, list)
        return res

    return res


# if __name__ == "__main__":
#     raw_tsvs = accession_to_prot_info(
#         ids, ["upi", "accession", "organism", "protein"], knowledge_base="uniparc"
#     )
#     assert isinstance(raw_tsvs, list)
#     print(process_tsvs(raw_tsvs))