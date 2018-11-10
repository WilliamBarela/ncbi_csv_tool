class Summary:
    '''
    Takes in a row from the summary_csv list of lists and creates and object
    '''

    def init(self, summary_row):
        self.experiment_accession
        self.experiment_title
        self.organism_name
        self.instrument
        self.submitter
        self.study_accession
        self.study_title
        self.sample_accession
        self.sample_title
        self.total_size_mb
        self.total_runs
        self.total_spots
        self.total_bases
        self.library_name
        self.library_strategy
        self.library_source
        self.library_selection

class Runinfo:
    '''
    Takes in a row from the runinfo_csv list of lists and creates and object
    '''

    def init(self, runinfo_row):
        self.run
        self.release_date
        self.load_date
        self.spots
        self.bases
        self.spots_with_mates
        self.avg_length
        self.size_mb
        self.assembly_name
        self.download_path
        self.experiment
        self.library_name
        self.library_strategy
        self.library_selection
        self.library_source
        self.library_layout
        self.insert_size
        self.insert_dev
        self.platform
        self.model
        self.sra_study
        self.bio_project
        self.study_pubmed_id
        self.project_id
        self.sample
        self.bio_sample
        self.sample_type
        self.tax_id
        self.scientific_name
        self.sample_name
        self.g1k_pop_code
        self.source
        self.g1k_analysis_group
        self.subject_id
        self.sex
        self.disease
        self.tumor
        self.affection_status
        self.analyte_type
        self.histological_type
        self.body_site
        self.center_name
        self.submission
        self.dbgap_study_accession
        self.consent
        self.run_hash
        self.read_hash

class Projects:
    '''
    Takes in summary and runinfo object and joins them into a project object
    '''

    def init(self, summary_object, runinfo_object):
        '''add code'''

"""

"""
