from setuptools import find_packages, setup

setup(
    name="lahuta",
    version="0.5",
    description="Lahuta is a Python package for calculating interatomic contacts based on the rules defined in CREDO.",
    author="Besian I. Sejdiu",
    project_urls={
        "Source code": "https://github.com/bisejdiu/lahuta",
        "Documentation": "https://github.com/bisejdiu/lahuta",
    },
    author_email="besian.sejdiu@stjude.org",
    # license="", # none yet
    url="https://github.com/bisejdiu/lahuta",
    install_requires=[],
    # install_requires=["MDAnalysis", "openbabel"],
    packages=find_packages(),
    python_requires=">=3.8",
    tests_require=[
        "pytest",
    ],
)
