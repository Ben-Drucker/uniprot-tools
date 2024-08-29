#!/bin/bash

# Run from root directory of repo

coverage run --branch --source=. -m unittest tests.tests
coverage html --skip-empty
