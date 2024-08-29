#!/bin/bash

# Run from level up from the root directory of repo

coverage run --branch --source=uniprot_tools -m unittest uniprot_tools.tests.tests
coverage html --skip-empty
