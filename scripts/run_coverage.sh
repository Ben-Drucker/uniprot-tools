#!/bin/bash

# Run from level up from the root directory of repo

python3.12 -m coverage run --branch --source=PkgContainer/uniprot_tools -m unittest PkgContainer.tests.tests
