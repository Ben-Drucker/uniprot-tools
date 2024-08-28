# Check java is installed
import os, subprocess

import requests
from setuptools import find_packages, setup

try:
    subprocess.check_output(["java", "-version"])
except subprocess.CalledProcessError:
    raise RuntimeError(
        f"Java is not installed. Please install Java and try again:"
        f" \n\thttps://www.java.com/en/download/help/download_options.html"
    )

os.makedirs("uniprot_tools/bin", exist_ok=True)
with requests.get(
    "https://research.bioinformatics.udel.edu/peptidematch/downloads/PeptideMatchCMD_1.1.jar",
    timeout=15,
) as r:
    with open("uniprot_tools/bin/PeptideMatchCMD_1.1.jar", "wb") as f:
        f.write(r.content)

setup(packages=find_packages(include=["uniprot_tools"]), package_data={"uniprot_tools": ["bin/*"]})
