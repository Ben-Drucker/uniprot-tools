#!/bin/bash

# Run from level up from the root directory of repo
export COVERAGE_PROCESS_START="/Users/druc594/Desktop/Uniprot Tools/.coveragerc"
python3.12 -m coverage run --branch --source=PkgContainer/uniprot_tools -m unittest PkgContainer.tests.tests
