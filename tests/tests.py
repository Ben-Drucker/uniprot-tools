import os, unittest


def compare_files(file1, file2):
    with open(file1) as f1, open(file2) as f2:
        return f1.read() == f2.read()


class TestConcurrency(unittest.TestCase):
    def setUp(self) -> None:
        from uniprot_tools.concurrency_tools import TqdmParallel

        return super().setUp()

    def test_tqdm_parallel(self): ...


class TestFasta(unittest.TestCase):
    def tearDown(self) -> None:
        for file in os.listdir("tests/test_outputs"):
            os.remove(f"tests/test_outputs/{file}")
        return super().tearDown()

    def test_read_fasta(self):
        import contextlib, io, json

        from uniprot_tools.fasta_tools import read_fasta

        no_sort = read_fasta(
            "tests/ground_truth/write_no_sort_yes_duplicates.fasta", sort_keys=False
        )
        yes_sort = read_fasta("tests/ground_truth/write_no_sort_yes_duplicates.fasta")

        test_output = io.StringIO()
        with contextlib.redirect_stderr(test_output):
            yes_sort_prog = read_fasta(
                "tests/ground_truth/write_no_sort_yes_duplicates.fasta",
                show_progress=True,
            )
        self.assertRegex(
            test_output.getvalue(),
            r"Bytes Read: 100%\|#+\| 684/684",
        )

        with open("tests/ground_truth/read_fasta_test_no_sort.json") as f:
            ground_truth_no_sort = json.load(f)
        with open("tests/ground_truth/read_fasta_test_yes_sort.json") as f:
            ground_truth_yes_sort = json.load(f)

        self.assertEqual(list(no_sort.keys()), list(ground_truth_no_sort.keys()))
        self.assertEqual(list(no_sort.values()), list(ground_truth_no_sort.values()))
        self.assertEqual(list(yes_sort.keys()), list(ground_truth_yes_sort.keys()))
        self.assertEqual(list(yes_sort.values()), list(ground_truth_yes_sort.values()))
        self.assertEqual(list(yes_sort_prog.keys()), list(ground_truth_yes_sort.keys()))
        self.assertEqual(list(yes_sort_prog.values()), list(ground_truth_yes_sort.values()))

    def test_write_fasta(self):
        import json

        from uniprot_tools.fasta_tools import write_fasta

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
