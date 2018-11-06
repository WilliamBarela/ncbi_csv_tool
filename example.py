from ncbi_csv_tool import query as q

"""
    Be sure to move this file one folder up so that ncbi_csv_tool directory is on the same level.
"""

summary_csv_list, run_csv_list = q.get_csv_lists("your.email@domain.com","sra","heterodera[Organism]")
