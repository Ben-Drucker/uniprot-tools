import unittest


def compare_files(file1, file2):
    with open(file1) as f1, open(file2) as f2:
        return f1.read() == f2.read()


class TestConcurrency(unittest.TestCase):
    def setUp(self) -> None:
        from uniprot_tools.concurrency_tools import TqdmParallel

        return super().setUp()

    def test_tqdm_parallel(self): ...


class TestFasta(unittest.TestCase):

    def test_read_fasta(self): ...

    def test_write_fasta(self):
        from uniprot_tools.fasta_tools import write_fasta
        import json

        with open("tests/test_data/test_peptides.json") as f:
            d = json.load(f)
            descriptions, peptides = list(d.keys()), list(d.values())
            file_names = [
                "tests/test_outputs/write_no_sort_yes_duplicates.fasta",
                "tests/test_outputs/write_yes_sort_yes_duplicates.fasta",
                "tests/test_outputs/write_yes_sort_yes_duplicates_sort_by_seq.fasta",
                "tests/test_outputs/write_no_sort_no_duplicates.fasta",
                "tests/test_outputs/write_yes_sort_no_duplicates.fasta",
                "tests/test_outputs/write_yes_sort_no_duplicates_sort_by_seq.fasta",
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
        from uniprot_tools.pep_search import (
            create_haystacks,
            create_index,
            pep_search,
            post_process_java_output,
            process_tsvs,
        )

        return super().setUp()

    def test_pep_search(self): ...

    def test_create_haystacks(self): ...

    def test_create_index(self): ...

    def test_post_process_java_output(self): ...


class TestGetInfo(unittest.TestCase):
    def setUp(self) -> None:
        from uniprot_tools.get_info import accession_to_prot_info

        return super().setUp()

    def test_accession_to_prot_info(self): ...
