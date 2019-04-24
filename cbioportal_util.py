"""A command line tool that retrieves cancer genomic data from CBioPortal, and summarizes the results."""
import requests
import csv
import StringIO
import sys
import subprocess
from os import system, name

class CaseSet(object):
	cases = {}
	summary = {}
	meta_data = {}
	gene_meta_data = {}
	genetic_profiles = ['gbm_tcga_gistic', 'gbm_tcga_mutations']

	def __init__(self, gene_list):
		"""Initializes the case-set object."""
		for genetic_profile_id in self.genetic_profiles:
			self.add_genetic_profile(genetic_profile_id, gene_list)

	def add_genetic_profile(self, genetic_profile_id, gene_list):
		"""Adds genetic profile data to file object."""
		meta_data_fields = ['COMMON', 'GENE_ID']
		url = self.get_request_payload(genetic_profile_id, gene_list)
		file_obj = self.fetch_API_response(url)

		i_position = 0
		for line in iter(file_obj.readline, ''):
			if line[0] == '#':
				i_position = file_obj.tell()
			else:
				break
		file_obj.seek(i_position)
		reader = csv.DictReader(file_obj, delimiter='\t')
		for row in reader:
			gene = row['COMMON']
			for field in row:
				if field in meta_data_fields:
					if not (gene in self.gene_meta_data):
						self.gene_meta_data[gene] = {}
					self.gene_meta_data[gene][field] = row[field]
				else:
					if not (field in self.cases):
						self.cases[field] = {}
					if not (gene in self.cases[field]):
						self.cases[field][gene] = {}
					self.cases[field][gene][genetic_profile_id] = row[field]
		file_obj.close()

	def get_request_payload(self, genetic_profile_id, gene_list):
		"""Gets payload body from cBioPortal's API."""
		query_payload = {}
		query_payload['gene_list'] = gene_list
		query_payload['genetic_profile_id'] = genetic_profile_id
		return 'http://www.cbioportal.org/webservice.do?cmd=getProfileData&genetic_profile_id={genetic_profile_id}&id_type=gene_symbol&gene_list={gene_list}&case_set_id=gbm_tcga_cnaseq'.format(**query_payload)

	def fetch_API_response(self, url):
		"""Gets the file object retrieved from an HTTP request."""
		response = requests.get(url)
		file_obj = StringIO.StringIO(response.content)
		return file_obj

	def get_mutation_percent(self, gene):
		return float(self.summary[gene]['mutated_case_count'])/float(self.summary[gene]['total_case_count']) # % of mutations

	def get_copy_no_altered_percent(self, gene):	
		return float(self.summary[gene]['copy_no_alteration_case_count'])/float(self.summary[gene]['total_case_count']) # % of copy-number alterations

	def get_all_altered_percent(self, gene):
		return float(self.summary[gene]['copy_no_alteration_case_count']+self.summary[gene]['mutated_case_count']-self.summary[gene]['multiple_alterations_case_count'])/float(self.summary[gene]['total_case_count']) # % of mutations OR copy-number alterations

	def is_mutated(self, case, gene):
		"""Checks if the gene is mutated."""
		genetic_profile_id = self.genetic_profiles[1]
		profile_result = case[gene][genetic_profile_id] 
		return False if (profile_result == 'NaN' or str(profile_result) == '0') else True

	def is_copy_no_altered(self, case, gene):
		"""Checks for the gene's copy number alteration."""		
		genetic_profile_id = self.genetic_profiles[0]
		profile_result = case[gene][genetic_profile_id]
		return True if (int(profile_result) in [+2, -2]) else False

	def set_summary(self):
		"""Summerizes alterations of the queried genes."""
		mutated_case_count = 0
		copy_no_alteration_case_count = 0
		multiple_alterations_case_count = 0		
		self.total_case_count = len(self.cases)

		for gene in self.gene_meta_data:
			for case in self.cases:
				if self.is_copy_no_altered(self.cases[case], gene):
					copy_no_alteration_case_count += 1
				elif self.is_mutated(self.cases[case], gene):
					mutated_case_count += 1
				elif self.is_copy_no_altered(self.cases[case], gene) and self.is_mutated(self.cases[case], gene):
					multiple_alterations_case_count += 1
			
			mutated_case_percent = float(mutated_case_count)/float(self.total_case_count)
			copy_no_alteration_percent = float(copy_no_alteration_case_count)/float(self.total_case_count)
			all_altered_percent = float(copy_no_alteration_case_count+mutated_case_count-multiple_alterations_case_count)/float(self.total_case_count)
			self.summary[gene] = {'total_case_count': self.total_case_count, 'mutated_case_count': mutated_case_count, 'copy_no_alteration_case_count': copy_no_alteration_case_count, 'multiple_alterations_case_count': multiple_alterations_case_count, 'all_altered_percent': all_altered_percent}

	def display_case_set_summary(self, display_mode):
		""""Displays the summarized analysis of queried genes."""
		display = [{'text': '{gene} is mutated in {result}% of all cases.', 'summary_function': self.get_mutation_percent, 'display_modes': ['verbose']}, {'text': '{gene} is copy number altered in {result}% of all cases.', 'summary_function': self.get_copy_no_altered_percent, 'display_modes': ['verbose']}, {'text': 'Total % of cases where {gene} is altered by either mutation or copy number alteration: {result}% of all cases.', 'summary_function': self.get_all_altered_percent, 'display_modes': ['verbose', 'concise']}]
		
		self.set_summary()
		calc = {}
		for gene in self.summary:
			calc = {'gene': gene}
			for value in display:
				if display_mode in value['display_modes']:
					result_p = value['summary_function'](gene)
					calc['result'] = '{0:.{1}f}'.format(result_p*100, 0) # percentage
					print (value['text'].format(**calc))

def main():
	#subprocess.call([sys.executable, "-m", "pip", "install", "requests"])
	if (len(sys.argv) not in [2, 3, 4]):
		#if name == 'nt':
		#	system('cls')
		#else:
		#system('clear')
		print('____________________________________________________________________________________')
		print('\nCBioPortal API Util: A Command line tool that retrieves cancer genomic data from the \ncBio Cancer Genomics Portal (http://www.cbioportal.org/), which currently supports a\nREST-based web API, and summarizes the results.\n')
		print('Usage: cbioportal_api_util.py gene1 [gene2] [gene3]')
	else:
		if (len(sys.argv) > 2):
			display_mode = 'concise'
		else:
			display_mode = 'verbose'
		gene_list = ','.join(sys.argv[1:])
		query = CaseSet(gene_list)
		query.display_case_set_summary(display_mode)

if __name__ == '__main__':
	main()