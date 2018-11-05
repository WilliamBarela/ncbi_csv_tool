# ncbi_csv_tool
Collects CSV files from NCBI search results using Entrez, the NCBI API

## How to use:
Download ncbi_cvs_tool to root directory of your project.

To get a list of lists in csv format, in your Python script include the following:

```python
from ncbi_csv_tool import query

csv_list = query.get_csv("your.email@domain.edu", "selected_ncbi_database", "search_terms")
```

## Python requirements
This library will work with any 3.xx version of Python
