import io
import json
import pandas as pd
from .concurrency_tools import TqdmParallel
from typing import Any
from xml.etree import ElementTree


def _unstringify(string: str):
    try:  # xml
        return ElementTree.fromstring(string)
    except ElementTree.ParseError:
        try:  # json
            jl: dict[Any, Any] = json.loads(string)
            return jl
        except json.decoder.JSONDecodeError:
            try:  # tsv
                return pd.read_csv(io.StringIO(string), sep="\t")
            except pd.errors.ParserError:

                raise ValueError(
                    f"Could not determine string type of {string}; it could not be parsed as XML,"
                    " JSON, or TSV (the allowed types)."
                )


def process_uniprot_output(string_structures: list[str]):
    """Turn Uniprot text TSVs, JSONs, or XMLs into one, combined pandas `pd.DataFrame`, `dict`, or \
        `ElementTree`.

    Parameters
    ----------
    ``string_structures`` :
        List of TSV, JSON, or XML strings. Note: if TSVs, each TSV must have the same columns.

    Returns
    -------
        Table of combined results.
    """
    individual_results = TqdmParallel.tqdm_starmap(
        _unstringify,
        [(s,) for s in string_structures],
        processes_or_threads="threads",
    )
    assert individual_results
    representative = individual_results[0]
    match representative:
        case pd.DataFrame():
            to_concat: list[pd.DataFrame] = []
            for x in individual_results:
                assert isinstance(x, pd.DataFrame)
                to_concat.append(x)
            return pd.concat(to_concat)
        case dict():
            to_merge: list[dict] = []
            for x in individual_results:
                assert isinstance(x, dict)
                to_merge.append(x)
            for e, d in enumerate(to_merge):
                d[f"results-{e}"] = d.pop("results")
            return {k: v for d in individual_results for k, v in d.items()}
        case ElementTree.Element():
            return individual_results[0]
        case _:  # type: ignore
            raise ValueError(f'Could not determine type of the input: "{individual_results[0]}"')
