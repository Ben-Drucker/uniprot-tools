import collections
import contextlib
import gzip
import io
import json
import numpy as np
import os
import pandas as pd
import requests
import shutil
import tqdm
import unittest


def compare_files(file1, file2):
    with open(file1) as f1, open(file2) as f2:
        return f1.read() == f2.read()


class TestConcurrency(unittest.TestCase):
    def setUp(self) -> None:
        return super().setUp()

    def test_tqdm_parallel(self): ...


class TestFasta(unittest.TestCase):
    def tearDown(self) -> None:
        # for file in os.listdir("PkgContainer/tests/test_outputs"):
        #     os.remove(f"PkgContainer/tests/test_outputs/{file}")
        return super().tearDown()

    def test_read_fasta(self):
        from ..uniprot_tools.fasta_tools import read_fasta

        no_sort = read_fasta(
            "PkgContainer/tests/ground_truth/write_no_sort_yes_duplicates.fasta", sort_keys=False
        )
        yes_sort = read_fasta("PkgContainer/tests/ground_truth/write_no_sort_yes_duplicates.fasta")

        test_output = io.StringIO()
        with contextlib.redirect_stderr(test_output):
            yes_sort_prog = read_fasta(
                "PkgContainer/tests/ground_truth/write_no_sort_yes_duplicates.fasta",
                show_progress=True,
            )
        self.assertRegex(
            test_output.getvalue(),
            r"Bytes Read: 100%\|#+\| 684/684",
        )

        with open("PkgContainer/tests/ground_truth/read_fasta_test_no_sort.json") as f:
            ground_truth_no_sort = json.load(f)
        with open("PkgContainer/tests/ground_truth/read_fasta_test_yes_sort.json") as f:
            ground_truth_yes_sort = json.load(f)

        self.assertEqual(list(no_sort.keys()), list(ground_truth_no_sort.keys()))
        self.assertEqual(list(no_sort.values()), list(ground_truth_no_sort.values()))
        self.assertEqual(list(yes_sort.keys()), list(ground_truth_yes_sort.keys()))
        self.assertEqual(list(yes_sort.values()), list(ground_truth_yes_sort.values()))
        self.assertEqual(list(yes_sort_prog.keys()), list(ground_truth_yes_sort.keys()))
        self.assertEqual(list(yes_sort_prog.values()), list(ground_truth_yes_sort.values()))

    def test_write_fasta(self):
        from ..uniprot_tools.fasta_tools import write_fasta

        with open("PkgContainer/tests/test_data/test_peptides.json") as f:
            d = json.load(f)
            descriptions, peptides = list(d.keys()), list(d.values())
            file_names = [
                "PkgContainer/tests/test_outputs/write_no_sort_yes_duplicates.fasta",
                "PkgContainer/tests/test_outputs/write_yes_sort_yes_duplicates.fasta",
                "PkgContainer/tests/test_outputs/write_yes_sort_yes_duplicates_sort_by_seq.fasta",
                "PkgContainer/tests/test_outputs/write_no_sort_no_duplicates.fasta",
                "PkgContainer/tests/test_outputs/write_yes_sort_no_duplicates.fasta",
                "PkgContainer/tests/test_outputs/write_yes_sort_no_duplicates_sort_by_seq.fasta",
            ]
            write_fasta(peptides, descriptions, file_names[0], do_sort=False, only_unique=False)
            write_fasta(peptides, descriptions, file_names[1], do_sort=True, only_unique=False)
            write_fasta(
                peptides,
                descriptions,
                file_names[2],
                do_sort=True,
                only_unique=False,
                sort_by_desc_or_seq="seq",
            )
            write_fasta(peptides, descriptions, file_names[3], do_sort=False, only_unique=True)
            write_fasta(peptides, descriptions, file_names[4], do_sort=True, only_unique=True)
            write_fasta(
                peptides,
                descriptions,
                file_names[5],
                do_sort=True,
                only_unique=True,
                sort_by_desc_or_seq="seq",
            )
            ground_truth_files = [x.replace("test_outputs", "ground_truth") for x in file_names]
            for output, truth in zip(file_names, ground_truth_files):
                self.assertTrue(compare_files(output, truth))


class TestPepSearch(unittest.TestCase):
    def setUp(self) -> None:
        self._get_test_fasta()

        return super().setUp()

    # def tearDown(self) -> None:
    #     for file in os.listdir("PkgContainer/tests/test_outputs"):
    #         os.remove(f"PkgContainer/tests/test_outputs/{file}")

    #     return super().tearDown()

    def _get_test_fasta(self):
        url = (
            "https://ftp.uniprot.org/pub/databases/uniprot"
            "/current_release/knowledgebase/complete/uniprot_sprot"
            ".fasta.gz"
        )

        with requests.get(url, stream=True, timeout=120) as r:
            r.raise_for_status()
            with open("PkgContainer/tests/test_data/test_uniprot_sprot.fasta.gz", "wb") as f:
                for chunk in tqdm.tqdm(
                    r.iter_content(chunk_size=2**20),
                    total=int(np.ceil(int(r.headers["Content-Length"]) / 2**20)),
                    desc="Downloading test data from Uniprot.",
                    unit="MB",
                ):
                    f.write(chunk)

        with open(
            "PkgContainer/tests/test_data/test_uniprot_sprot-plaintext.fasta", "wt"
        ) as plain_f:
            with gzip.open(
                "PkgContainer/tests/test_data/test_uniprot_sprot.fasta.gz", "rt"
            ) as gzip_f:
                plain_f.write(gzip_f.read())

        os.remove("PkgContainer/tests/test_data/test_uniprot_sprot.fasta.gz")

    def _wrapper_create_haystacks(self):
        from ..uniprot_tools.pep_search import create_haystacks

        test_fasta = "PkgContainer/tests/test_data/test_uniprot_sprot-plaintext.fasta"
        with open(test_fasta) as f:
            num_lines = sum(1 for _ in f)

        haystack_dir = "PkgContainer/tests/test_outputs/haystacks"
        if os.path.exists(haystack_dir):
            shutil.rmtree(haystack_dir)
        os.makedirs(haystack_dir, exist_ok=True)
        create_haystacks([test_fasta], haystack_dir, num_lines=[num_lines], seqs_per_chunk=int(1e5))

        total_haystacks_size = sum(
            [os.path.getsize(f"{haystack_dir}/{x}") for x in os.listdir(haystack_dir)]
        )
        num_haystacks_files = [x for x in os.listdir(haystack_dir) if x.endswith(".fasta")]

        self.assertTrue(  # account for variations in size as uniprot updates
            250 * 1e6 < total_haystacks_size < 350 * 1e6,
            f"total_haystacks_size was {total_haystacks_size / 1e6:.2f} MB",
        )
        self.assertTrue(5 < len(num_haystacks_files) < 8)

        os.remove("PkgContainer/tests/test_data/test_uniprot_sprot-plaintext.fasta")

    # @unittest.skip(
    #     "Skipping `test_create_index` as it's run automatically by `test_post_process_java_output`"
    # )
    def _wrapper_create_index(self):
        from ..uniprot_tools.pep_search import create_index

        self._wrapper_create_haystacks()
        os.makedirs("PkgContainer/tests/test_outputs/indexes", exist_ok=True)
        create_index(
            [
                f"PkgContainer/tests/test_outputs/haystacks/{x}"
                for x in os.listdir("PkgContainer/tests/test_outputs/haystacks")
                if x.endswith(".fasta")
            ],
            "PkgContainer/tests/test_outputs/indexes",
        )

    # @unittest.skip(
    #     "Skipping `test_pep_search` as it's run automatically by `test_post_process_java_output`"
    # )
    def _wrapper_pep_search(self):
        from ..uniprot_tools.pep_search import pep_search

        self._wrapper_create_index()
        os.makedirs("PkgContainer/tests/test_outputs/pep_searches", exist_ok=True)
        with open("PkgContainer/tests/test_data/test_peptides_2.json") as f:
            pep_search(
                json.load(f),
                "PkgContainer/tests/test_outputs/indexes",
                intermediate_dir="PkgContainer/tests/test_outputs/pep_searches",
                delete_index_when_done=False,
            )

    def test_post_process_java_output(self):
        from ..uniprot_tools.pep_search import post_process_java_output

        self._wrapper_pep_search()
        mapping = post_process_java_output(
            files=[
                f"PkgContainer/tests/test_outputs/pep_searches/{x}"
                for x in os.listdir("PkgContainer/tests/test_outputs/pep_searches")
                if x.endswith(".pepsearchres")
            ],
        )

        self.assertTrue(len(mapping) > 0)
        with open("PkgContainer/tests/test_outputs/mapping.json", "w") as f:
            json.dump({k: sorted(v) for k, v in mapping.items()}, f)


class TestGetInfo(unittest.TestCase):
    def setUp(self) -> None:

        with open("PkgContainer/tests/test_data/test_ids.json") as f:
            ids = json.load(f)
            self.uniprotkb_ids = ids["uniprotkb_ids"]
            self.uniparc_ids = ids["uniparc_ids"]
            self.taxonomy_ids = ids["taxonomy_ids"]

        self.uniprotkb_cols = [
            "accession",
            "reviewed",
            "protein_name",
            "gene_names",
            "organism_name",
            "cc_alternative_products",
        ]

        self.uniparc_cols = [
            "upi",
            "accession",
            "organism",
            "protein",
        ]

        self.taxonomy_cols = [
            "id",
            "common_name",
            "scientific_name",
            "lineage",
        ]

        return super().setUp()

    def _url_test(self, knowledge_base):
        # from ..uniprot_tools.get_info import accession_to_prot_info
        from ..uniprot_tools import get_info

        accession_to_prot_info = get_info.accession_to_prot_info

        columns = self.get_info_cols(knowledge_base)

        urls = accession_to_prot_info(
            self.__dict__[f"{knowledge_base}_ids"],
            self.__dict__[f"{knowledge_base}_cols"],
            knowledge_base=knowledge_base,
            get_urls_only=True,
            columns=columns,
        )
        assert isinstance(urls, list)
        with open(f"PkgContainer/tests/ground_truth/{knowledge_base}_urls.json") as f:
            ground_truth_urls = json.load(f)

        flat_urls = "".join(urls)
        flat_ground_truth_urls = "".join(ground_truth_urls)

        urls_counter = collections.Counter(flat_urls)
        ground_truth_counter = collections.Counter(flat_ground_truth_urls)
        if set(urls_counter.keys()) != set(ground_truth_counter.keys()):
            for k in set(urls_counter.keys()) - set(ground_truth_counter.keys()):
                ground_truth_counter[k] = 0
            for k in set(ground_truth_counter.keys()) - set(urls_counter.keys()):
                urls_counter[k] = 0

        for (k1, v1), (k2, v2) in zip(
            sorted(urls_counter.items()), sorted(ground_truth_counter.items())
        ):
            self.assertEqual(k1, k2)
            self.assertTrue(
                abs(v1 - v2) < 0.1 * min(v1, v2),
                f"Proportion of {k1} was {v1} out of {min(v1, v2)} characters, which is high.",
            )


    @staticmethod
    def get_info_cols(knowledge_base):
        match knowledge_base:
            case "uniprotkb":
                columns = [
                    "accession",
                    "reviewed",
                    "protein_name",
                    "gene_names",
                    "organism_name",
                    "cc_alternative_products",
                ]
            case "uniparc":
                columns = [
                    "upi",
                    "accession",
                    "organism",
                    "protein",
                ]
            case "taxonomy":
                columns = [
                    "id",
                    "common_name",
                    "scientific_name",
                    "lineage",
                ]
            case _:
                raise ValueError(f"Unknown knowledge base: {knowledge_base}")

        return columns

    def _info_test(self, knowledge_base):
        from ..uniprot_tools.get_info import accession_to_prot_info
        from ..uniprot_tools.pep_search import process_uniprot_output
        from pandas.testing import assert_frame_equal

        columns = self.get_info_cols(knowledge_base)

        this_table = accession_to_prot_info(
            self.__dict__[f"{knowledge_base}_ids"],
            self.__dict__[f"{knowledge_base}_cols"],
            columns=columns,
            knowledge_base=knowledge_base,
            process_output_into_DataFrame=True,
        )
        # Taxon Id,Common name,Scientific name,Lineage
        # Entry,UniProtKB,Organisms,Protein names
        # Entry,Reviewed,Protein names,Gene Names,Organism,Alternative products (isoforms)

        self.assertIsInstance(this_table, pd.DataFrame)
        assert isinstance(this_table, pd.DataFrame)  # Redundant, but for the linter

        ground_truth_table = pd.read_csv(
            f"PkgContainer/tests/ground_truth/{knowledge_base}_info_table.csv"
        )

        this_table_rows = set(this_table.to_records(index=False).tolist())
        ground_truth_rows = set(ground_truth_table.to_records(index=False).tolist())
        jaccard = len(this_table_rows.intersection(ground_truth_rows)) / len(
            this_table_rows.union(ground_truth_rows)
        )
        self.assertTrue(
            jaccard > 0.80,
            f"Jaccard similarity of rows was {jaccard:.2f}, which is low. Note that a score less"
            " than 1.0 may be due to changes in Uniprot.",
        )

    def test_uniprotkb_urls(self):
        self._url_test("uniprotkb")

    def test_uniparc_urls(self):
        self._url_test("uniparc")

    def test_taxonomy_urls(self):
        self._url_test("taxonomy")

    def test_uniprotkb_info(self):
        self._info_test("uniprotkb")

    def test_uniparc_info(self):
        self._info_test("uniparc")

    def test_tax_info(self):
        self._info_test("taxonomy")
