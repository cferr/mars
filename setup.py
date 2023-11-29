from setuptools import setup, find_packages
import sys


def parse_requirements():
    reqs = []
    with open("requirements.txt", "r") as handle:
        for line in handle:
            if line.startswith("-") or line.startswith("git+"):
                continue
            idx = 0
            for i, c in enumerate(line):
                if c == "/":
                    idx = i + 1
            reqs.append(line[idx:])
    return reqs

setup(
    name="mars",
    version="0.0.2",
    description="Computation of Maximal Atomic irRedundant Sets",
    author="Corentin Ferry",
    url="https://github.com/cferr/mars.git",
    packages=find_packages(),
    include_package_data=True,
    long_description=open("README.md").read(),
    install_requires=parse_requirements()
)
