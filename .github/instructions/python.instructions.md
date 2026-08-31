---
name: 'Python Standards'
description: 'Coding conventions for Python files'
applyTo: '**/*.py'
---
# Python coding standards
- Follow the PEP 8 style guide.
- Use type hints for all function signatures.
- Write docstrings for public functions.
- Use tab for indentation.
- use generators and iterators where appropriate to save memory.
- use pandas for data manipulation and analysis.
- use biopython and pysam for biological sequence analysis.
- make use of the assert statement for debugging and testing.
- use logging for error handling and debugging.
- create an argparse interface for command line scripts.
- make use of context managers for resource management (e.g., file handling).
- avoid the use of iterrows() in pandas, use vectorized operations instead.
- when writing to csv or excel files, make sure to not include the index using a pandas dataframe.