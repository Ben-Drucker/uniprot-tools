import re, requests, tqdm, warnings
from .concurrency_tools import TqdmParallel, request_worker
from typing import Literal


def accession_to_prot_info(
    accessions: list[str] | list[int],
    format_: Literal["tsv", "json"] = "tsv",
    columns: list[str] | None = None,
    knowledge_base: Literal["uniprotkb", "uniparc", "taxonomy"] = "uniprotkb",
    num_simultaneous_requests: int = 100,
    get_urls_only: bool = False,
    compressed_download: bool = False,
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
        `ValueError` :
            If ``knowledge_base`` is not one of ``uniprotkb`` or ``uniparc`` or ``taxonomy``

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
            raise ValueError("`database` must be either `uniprotkb`, `uniparc`, or `tax`.")
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
        if columns is None:
            assert format_ == "json", "`columns` must be specified if `format` is not `tsv`."
            fmt = "format=json"
            fields_str = ""
        else:
            fields_str = "fields=" + ",".join(columns) + "&"
            fmt = "format=tsv"
        url = (
            f"https://rest.uniprot.org/{knowledge_base}/search?"
            f"{fields_str}{query_str}&{fmt}&size={len(chunk)}"
            f"&compressed={str(bool(compressed_download)).lower()}"
        )
        urls.append(url)
    assert len(urls) >= len(accessions) // CHUNK_SIZE
    if get_urls_only:
        return urls
    else:
        res = TqdmParallel.tqdm_starmap(
            worker_fn=request_worker,
            worker_args=[(url,) for url in urls],
            num_workers=num_simultaneous_requests,
            processes_or_threads="threads",
        )
    assert res
    new_res = []
    for res_i in res:
        if isinstance(res_i, requests.HTTPError):
            warnings.warn(str(res_i), category=UserWarning)
        else:
            new_res.append(str(res_i))
    return new_res


# if __name__ == "__main__":
#     raw_tsvs = accession_to_prot_info(
#         ids, ["upi", "accession", "organism", "protein"], knowledge_base="uniparc"
#     )
#     assert isinstance(raw_tsvs, list)
#     print(process_tsvs(raw_tsvs))
