from ncbi_csv_tool import query as q
from ncbi_csv_tool import generator as gen

"""
    Be sure to move this file one folder up so that ncbi_csv_tool directory is on the same level.
"""

summary_csv_list, run_csv_list = q.get_csv_lists("your.email@domain.com","sra","heterodera[Organism]")

summary_objects_list = gen.summary_objects(summary_csv_list)
runinfo_objects_list = gen.runinfo_objects(run_csv_list)
