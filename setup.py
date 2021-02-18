# File: setup.py
# Date: 14-Feb-2021
#
# Update:
#
import re

from setuptools import find_packages
from setuptools import setup

packages = []
thisPackage = "rcsb.utils.dictionary"

with open("rcsb/utils/dictionary/__init__.py", "r") as fd:
    version = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]', fd.read(), re.MULTILINE).group(1)

if not version:
    raise RuntimeError("Cannot find version information")

setup(
    name=thisPackage,
    version=version,
    description="RCSB Python Dictionary Utility Classes",
    long_description="See:  README.md",
    author="John Westbrook",
    author_email="john.westbrook@rcsb.org",
    url="https://github.com/rcsb/py-rcsb_utils_dictionary",
    #
    license="Apache 2.0",
    classifiers=(
        "Development Status :: 3 - Alpha",
        # 'Development Status :: 5 - Production/Stable',
        "Intended Audience :: Developers",
        "Natural Language :: English",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
    ),
    entry_points={},
    #
    install_requires=[
        "scipy",
        "numpy",
        "mmcif >= 0.57",
        "rcsb.utils.io >= 0.97",
        "rcsb.utils.config >= 0.35",
        "rcsb.utils.multiproc >= 0.17",
        "rcsb.utils.validation >= 0.22",
        "rcsb.utils.chemref >= 0.68",
        "rcsb.utils.citation >= 0.15",
        "rcsb.utils.ec >= 0.21",
        "rcsb.utils.taxonomy >= 0.32",
        "rcsb.utils.seq >= 0.43",
        "rcsb.utils.struct >= 0.26",
        "rcsb.utils.repository >= 0.11",
    ],
    packages=find_packages(exclude=["rcsb.utils.tests-dictionary", "rcsb.utils.tests-*", "tests.*"]),
    package_data={
        # If any package contains *.md or *.rst ...  files, include them:
        "": ["*.md", "*.rst", "*.txt", "*.cfg"]
    },
    #
    # These basic tests require no database services -
    test_suite="rcsb.utils.tests-dictionary",
    tests_require=["tox"],
    #
    # Not configured ...
    extras_require={"dev": ["check-manifest"], "test": ["coverage"]},
    # Added for
    command_options={"build_sphinx": {"project": ("setup.py", thisPackage), "version": ("setup.py", version), "release": ("setup.py", version)}},
    # This setting for namespace package support -
    zip_safe=False,
)
