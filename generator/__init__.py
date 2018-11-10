class Summary:
    '''
    Takes in a row from the summary_csv list of lists and creates and object
    '''

    def init(self, summary_row):
        self.experiment_accession = summary_row[0]
        self.experiment_title = summary_row[1]
        self.organism_name = summary_row[2]
        self.instrument = summary_row[3]
        self.submitter = summary_row[4]
        self.study_accession = summary_row[5]
        self.study_title = summary_row[6]
        self.sample_accession = summary_row[7]
        self.sample_title = summary_row[8]
        self.total_size_mb = summary_row[9]
        self.total_runs = summary_row[10]
        self.total_spots = summary_row[11]
        self.total_bases = summary_row[12]
        self.library_name = summary_row[13]
        self.library_strategy = summary_row[14]
        self.library_source = summary_row[15]
        self.library_selection = summary_row[16]
    
    # def srx_dict(self):


class Runinfo:
    '''
    Takes in a row from the runinfo_csv list of lists and creates and object
    '''

    def init(self, runinfo_row):
        self.run = runinfo_row[0]
        self.release_date = runinfo_row[1]
        self.load_date = runinfo_row[2]
        self.spots = runinfo_row[3]
        self.bases = runinfo_row[4]
        self.spots_with_mates = runinfo_row[5]
        self.avg_length = runinfo_row[6]
        self.size_mb = runinfo_row[7]
        self.assembly_name = runinfo_row[8]
        self.download_path = runinfo_row[9]
        self.experiment = runinfo_row[10]
        self.library_name = runinfo_row[11]
        self.library_strategy = runinfo_row[12]
        self.library_selection = runinfo_row[13]
        self.library_source = runinfo_row[14]
        self.library_layout = runinfo_row[15]
        self.insert_size = runinfo_row[16]
        self.insert_dev = runinfo_row[17]
        self.platform = runinfo_row[18]
        self.model = runinfo_row[19]
        self.sra_study = runinfo_row[20]
        self.bio_project = runinfo_row[21]
        self.study_pubmed_id = runinfo_row[22]
        self.project_id = runinfo_row[23]
        self.sample = runinfo_row[24]
        self.bio_sample = runinfo_row[25]
        self.sample_type = runinfo_row[26]
        self.tax_id = runinfo_row[27]
        self.scientific_name = runinfo_row[28]
        self.sample_name = runinfo_row[29]
        self.g1k_pop_code = runinfo_row[30]
        self.source = runinfo_row[31]
        self.g1k_analysis_group = runinfo_row[32]
        self.subject_id = runinfo_row[33]
        self.sex = runinfo_row[34]
        self.disease = runinfo_row[35]
        self.tumor = runinfo_row[36]
        self.affection_status = runinfo_row[37]
        self.analyte_type = runinfo_row[38]
        self.histological_type = runinfo_row[39]
        self.body_site = runinfo_row[40]
        self.center_name = runinfo_row[41]
        self.submission = runinfo_row[42]
        self.dbgap_study_accession = runinfo_row[43]
        self.consent = runinfo_row[44]
        self.run_hash = runinfo_row[45]
        self.read_hash = runinfo_row[46]

class Projects:
    '''
    Takes in summary and runinfo object and joins them into a project object
    '''

    def init(self, summary_object, runinfo_object):
        '''add code'''
