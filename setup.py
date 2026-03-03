from setuptools import setup, find_packages

setup(
    name="straln",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "PyVCF3",
        "pandas",
        "Levenshtein"
    ],
    entry_points={
        'console_scripts': [
            'straln=stretcher_parser.__main__:main',
        ],
    },
)