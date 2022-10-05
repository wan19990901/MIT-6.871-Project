"""
Oct. 4, 2022

BiopsyGradeExtractor
- fit
	display summary
- validate result
	user oriented evaluation of the algorithm


"""


import numpy as np
import pandas as pd
from tqdm import tqdm
import re

class BiopsyGradeExtractor:
	def __init__(self, centers = None): 
		if centers is not None:
			self.centers = centers
		else:
			self.centers = None

	def _func_det_negation(self, text_detect_neg):
		det_negation = 0
		for idx, val in enumerate(text_detect_neg):
			if 'no' == val.lower():
				bool_list = []
				for val_ in text_detect_neg[idx:]:
					if '.' not in val_:
						bool_list.append(False)
					else:
						bool_list.append(True)
				if not any(bool_list): # if there is no dot between 'no' --- 'query word', then we consider it as a negation
					det_negation = 1
		return det_negation

	def _extract_biopsy_grades(self):
		"""
		Iterate through reports and extract clinical concepts of interest
		"""
		directional_words = ['right', 'left','l','r']
		# declare lists which store abstracted biopsy features
		primary_grade_list = []; secondary_grade_list = []; overall_grade_list = []; overall_grade_merged_list = []; overall_gs_list = []; 

		# iterate through reports
		for text_oi, report_number in tqdm(zip(self.df_pathology_with_biopsy.Report_Text.values, self.df_pathology_with_biopsy.Report_Number.values), total = len(self.df_pathology_with_biopsy), desc = 'extracing info...'):
			# create lists which store abstracted biopsy features locally. "Local grade" correpsonds to grade of each tissue 
			primary_grade_list_local = []; secondary_grade_list_local = []; overall_grade_list_local = []; overall_gs_list_local = []
			
			# make all words in lowercase for consistency
			text_oi = text_oi.lower()

			# create a list of all the words in each report
			text_oi_list = text_oi.split(' ')
			found_directional_word = 0; 
			# iterate through words
			for idx, val in enumerate(text_oi_list):
				if any(dir_val in val for dir_val in directional_words):
					found_directional_word = 1

				primary_grade = None; secondary_grade = None;
				# Obtain overall grade (gleason)
				# , which appear after words like 'gleason' and directional words in each note. Note that we need to exclude false positive triggered
				# by M.D. Gleason
				if 'gleason' in val and found_directional_word and 'm.d.,' not in text_oi_list[idx - 2 : idx + 2]:
					found_primary_grade = 0; found_secondary_grade = 0; found_plus_sign = 0; found_one_dash = 0
					# get context around word "gleason"
					# heuristically choose 15. Gleason score we're after should be included in the next 15 words（might need to change this)
					raw_gleason_text = ''.join(text_oi_list[idx:idx + 15]) 

					# remove text in () or []: e.g. Score 3 (60%) + 4 (40%) = 7/10 -> Score 3 + 4 = 7/10
					raw_gleason_text = re.sub("[\(\[].*?[\)\]]", "", raw_gleason_text)
					# to capture gleason 3/5 where primary 3 and secondary 3, check if there's only one slash without any + or - signs. 
					slash_counter = 0; dash_counter = 0; plus_sign_counter = 0; only_one_slash = 0
					of_counter = 0
					# get context 
					raw_gleason_text_slash = ''.join(text_oi_list[idx:idx + 15]) 
					for char_dash in raw_gleason_text_slash:
						if '/' in char_dash:
							slash_counter += 1
						if '+' in char_dash:
							plus_sign_counter += 1
						if '-' in char_dash:
							dash_counter += 1
						if 'of' in char_dash:
							of_counter += 1
					if slash_counter == 1 and plus_sign_counter == 0 and dash_counter == 0 and of_counter == 0:
						only_one_slash = 1						

					# there are total 6 different ways to represent gleason scores. See below :
					for idx_char, char in enumerate(raw_gleason_text):
						# i) when overall grade is represented as primary + secondary
						if char == '+':
							raw_gleason_equation = raw_gleason_text[idx_char-10:idx_char+10]
							for sub_char in raw_gleason_equation:
								# checking + sign again; it covers a few corner cases 
								if sub_char == '+':
									found_plus_sign = 1
								if sub_char.isdigit():
									if not found_plus_sign:
										primary_grade_to_be = int(sub_char)
										found_primary_grade = 1
									elif found_primary_grade and found_plus_sign and not found_secondary_grade:
										primary_grade = primary_grade_to_be
										secondary_grade = int(sub_char)
										found_secondary_grade = 1
										break
							break
						elif char == '/' and not found_one_dash:
							# ii) when overall grade is represented as Gleason 3-4/5 
							raw_gleason_grade_dash = raw_gleason_text[idx_char-3:idx_char+2]
							if '-' in raw_gleason_grade_dash:
								raw_gleason_grade_dash_split = raw_gleason_grade_dash.split('/')
								if raw_gleason_grade_dash_split[1] == '5' and raw_gleason_grade_dash_split[0][0].isdigit() and raw_gleason_grade_dash_split[0][2].isdigit():
									primary_grade = int(raw_gleason_grade_dash_split[0][0])
									secondary_grade = int(raw_gleason_grade_dash_split[0][2])
									break
							# iii) when overall grade is represented as Gleason 6/10 
							raw_gleason_grade_dash_and_ten = raw_gleason_text[idx_char-1:idx_char+3]
							raw_gleason_grade_dash_and_ten_split = raw_gleason_grade_dash_and_ten.split('/')
							if raw_gleason_grade_dash_and_ten_split[0].isdigit() and raw_gleason_grade_dash_and_ten_split[1] == '10':
								# in this case primary grade = secondary grade = 3
								primary_grade = int(raw_gleason_grade_dash_and_ten_split[0])/2
								secondary_grade = primary_grade
								break

							# iv) when it is represented as Gleason grade 3/5 and 2/5 (where 3 is the primary and 2 is the secondary grade)
							raw_gleason_grade = raw_gleason_text[idx_char-1:idx_char+2]
							raw_gleason_grade_split = raw_gleason_grade.split('/')
							if len(raw_gleason_grade_split) == 2 and raw_gleason_grade_split[0].isdigit() and raw_gleason_grade_split[1] == '5':
								found_one_dash = 1
								primary_grade_to_be = int(raw_gleason_grade_split[0])								
								continue # thie ensures the next if statement is true automatically 
    					# iv - continued) encounter the second dash : most likely we have something like "gleason grade 3/5 and 2/5 confined to..." in the context
						if char == '/' and found_one_dash:
							raw_gleason_grade = raw_gleason_text[idx_char-1:idx_char+2]
							raw_gleason_grade_split = raw_gleason_grade.split('/')
							# iv - complete) when it is the form of x/5. We captured both primary and secondary grades
							if len(raw_gleason_grade_split) == 2 and raw_gleason_grade_split[0].isdigit() and raw_gleason_grade_split[1] == '5':
								primary_grade = primary_grade_to_be
								secondary_grade = int(raw_gleason_grade_split[0])
								break
						# v) when it is represented as gleason 3/5 (where primary 3 and secondary 3)
						elif found_one_dash and only_one_slash:
							primary_grade = primary_grade_to_be
							secondary_grade = primary_grade_to_be
							break
						# Another case for of
						if(char == 'of'):
							raw_gleason_grade = raw_gleason_text[idx_char-2:idx_char+1]
							raw_gleason_grade_split = raw_gleason_grade.split('of')
							if len(raw_gleason_grade_split) == 2 and (raw_gleason_grade_split[0].isdigit() or raw_gleason_grade_split[1].isdigit()):
								primary_grade = [int(s) for s in raw_gleason_grade_split[0].split() if s.isdigit()]
								secondary_grade = [int(s) for s in raw_gleason_grade_split[1].split() if s.isdigit()]
								if(len(secondary_grade) > 0 and len(primary_grade) > 0):
									if(secondary_grade[0] == 10):
										secondary_grade = None
										primary_grade = primary_grade[0]
										overall_grade_group = primary_grade
										break
									else:
										primary_grade = primary_grade[0]
										secondary_grade = secondary_grade[0]
										break
								elif(len(primary_grade) == 0):
									primary_grade = None
									secondary_grade = secondary_grade[0]
									overall_grade_group = secondary_grade
									break
								elif(len(secondary_grade) == 0):
									secondary_grade = None
									primary_grade = primary_grade[0]
									overall_grade_group = primary_grade
									break
								else:
									primary_grade = None
									secondary_grade = None
								
						# vi) when it is represented as z (x + y), where z is the total gleason score, x is the primary grade, and y is the secondary grade
						if char.isdigit() and idx_char < len(raw_gleason_text) - 5:
							if raw_gleason_text[idx_char + 1] == '(':							
								if raw_gleason_text[idx_char + 2].isdigit() and raw_gleason_text[idx_char + 3] == '+' and raw_gleason_text[idx_char + 4].isdigit():
									primary_grade = int(raw_gleason_text[idx_char + 2])
									secondary_grade = int(raw_gleason_text[idx_char + 4])
									break


									
				# get overall garde group based on abstracted gleason scores
				# https://www.pcf.org/about-prostate-cancer/diagnosis-staging-prostate-cancer/gleason-score-isup-grade/
				if primary_grade is not None and secondary_grade is not None:
					found_directional_word = 0 # set this to 0 so it can search other Gleason scores in the report
					if primary_grade + secondary_grade <= 6:
						overall_grade_group = 1
					elif primary_grade + secondary_grade == 7:
						if primary_grade == 3:
							overall_grade_group = 2
						elif primary_grade == 4:
							overall_grade_group = 3
						else: # 7/10 with no primary and seoncdary info
							overall_grade_group = 2
					elif primary_grade + secondary_grade == 8:
						overall_grade_group = 4
					elif primary_grade + secondary_grade > 8:
						overall_grade_group = 5
					else:
						overall_grade_group = None
					# update the grade group lists
					if overall_grade_group is not None:
						primary_grade_list_local.append(primary_grade)
						secondary_grade_list_local.append(secondary_grade)
						overall_grade_list_local.append(overall_grade_group)
						overall_gs_list_local.append(primary_grade + secondary_grade)
				else:
					overall_grade_group = None
			if(len(overall_grade_list_local) == 0):
				if('well-differentiated' in text_oi):
					overall_grade_list_local.append(1)
				elif('moderately differentiated' in text_oi):
					overall_grade_list_local.append(2)
				elif('moderately to poorly differentiated' in text_oi):
					overall_grade_list_local.append(3)
				elif('poorly differentiated' in text_oi):
					overall_grade_list_local.append(4)
				elif('gleason score is not considered reliable after hormonal treatment and thus not rendered' in text_oi):
					overall_grade_list_local.append(97)

			# store abstracted biopsy features
			# store local (i.e. each tissue) primary grades, secondary grades, overall grade groups
			if len(overall_grade_list_local) > 0:
				primary_grade_list.append(primary_grade_list_local)
				secondary_grade_list.append(secondary_grade_list_local)
				overall_grade_list.append(overall_grade_list_local)	
				# the maximum of local grades is the overall grade group, 
				overall_grade_merged_list.append(np.max(overall_grade_list_local))
			else:
				primary_grade_list.append(None)
				secondary_grade_list.append(None)
				overall_grade_list.append(None)	
				overall_grade_merged_list.append(None)

		# store abstracted biopsy features into the main df
		self.df_pathology_with_biopsy.loc[self.df_pathology_with_biopsy.index, 'primary_grade'] = primary_grade_list
		self.df_pathology_with_biopsy.loc[self.df_pathology_with_biopsy.index, 'secondary_grade'] = secondary_grade_list
		self.df_pathology_with_biopsy.loc[self.df_pathology_with_biopsy.index, 'overall_grade_merged'] = overall_grade_merged_list

		print('\n')
		print('Extraction complete.')
		print('\n')
		return

	def fit(self, df_pathology, display_summary = False, return_summary_df = False):
		"""
		inputs :

		"""
		# get all the reports with the following keywords : needle, core, biops, prost
		# idx_oi = [idx for idx, report in enumerate(df_pathology.Report_Text.values) if 'needle' in report.lower() or 'core' in report.lower() or 'biops' in report.lower() and 'prost' in report.lower()]
		# # exclude all the reports with the following keywords : prostatectomy, bone marrow
		# idx_to_exclude = [idx for idx, report in enumerate(df_pathology.Report_Text.values) if 'prostatectomy' in report.lower() or 'bone marrow' in report.lower(	)]# or 'nephrectomy' in report.lower() or 'cystectomy' in report.lower() or 'cystoprostatectomy' in report.lower() or 'ureter' in report.lower()] # or 'prostate, radical resection' in report.lower()]
		self.df_pathology_with_biopsy = df_pathology	
		# extract grade
		self._extract_biopsy_grades()

		if display_summary:
			print('Displaying summary...')
			primary_grade_all_list = []
			for val in self.df_pathology_with_biopsy.primary_grade.values:
				if val is not None:
					primary_grade_all_list += val
			df_count_stats_primary_grade = pd.value_counts(primary_grade_all_list, sort = True)		
			print('Counts per each primary grade: ')
			print(df_count_stats_primary_grade)
			print('\n')

			secondary_grade_all_list = []
			for val in self.df_pathology_with_biopsy.secondary_grade.values:
				if val is not None:
					secondary_grade_all_list += val
			df_count_stats_secondary_grade = pd.value_counts(secondary_grade_all_list, sort = True)		
			print('Counts per each secondary grade: ')
			print(df_count_stats_secondary_grade)
			print('\n')

			reports_mci_oi = set(self.df_pathology_with_biopsy.Report_Number.values)
			self.df_pathology_with_biopsy = self.df_pathology_with_biopsy.loc[self.df_pathology_with_biopsy.Report_Number.isin(reports_mci_oi)]

			df_count_stats_overall_grade = pd.value_counts(self.df_pathology_with_biopsy.overall_grade_merged.values, sort = True)
			print('Number of reports per each overall grade : ')
			print(df_count_stats_overall_grade)
			print('total = ', df_count_stats_overall_grade.sum())
			print('Number of unique biopsy reports : ', len(self.df_pathology_with_biopsy))
			# print('Number of unique biopsy reports per center : ')
			# print(pd.value_counts(self.df_pathology_with_biopsy.MRN_Type.values, sort = True))
			print('Number of unique patients : ', len(set(self.df_pathology_with_biopsy.EMPI.values)))

			df_pathology_with_biopsy_oi_overall_grade = self.df_pathology_with_biopsy.dropna(subset = ['overall_grade_merged'])
			print('Number of unique patients with overall_grade_merged : ', len(set(df_pathology_with_biopsy_oi_overall_grade.EMPI.values)))

		if return_summary_df:
			return self.df_pathology_with_biopsy
		else:
			return None

	def validate_results(self, feats_oi = None, report_numbers_to_query = None):
		"""
		inputs
		- feats_oi (list-type) : features of interest to validate
		- report_numbers_to_query (list-type) : specific report numbers to query. If not provided, each report is randomly sampled

		TODO : implement mechanism such that a user can obtain accuracy of the algorithm up to N reports they review
		"""
		count = 0; continue_eval = True
		if feats_oi is not None:
			feats_oi = ['Report_Number', 'Report_Text'] + list(feats_oi)
		else:
			# if feats_oi is not provided, then just evaluate overall grade group related info
			feats_oi = ['Report_Number', 'Report_Text', 'primary_grade', 'secondary_grade', 'overall_grade_merged']
		df_oi = self.df_pathology_with_biopsy[feats_oi].dropna(subset = feats_oi)
		while continue_eval:
			if report_numbers_to_query is not None:
				report_number_to_query = report_numbers_to_query[count]
				if report_number_to_query not in df_oi.Report_Number.values:
					print('Report ' + report_number_to_query + ' not found!')
					print('Randomly sampling a report...')
					sampled_entry = df_oi.sample(n = 1)
				else:
					sampled_entry = df_oi.loc[df_oi.Report_Number == report_number_to_query]
			else:
				sampled_entry = df_oi.sample(n = 1)
			print('\n')
			print('Extracted info : ')
			print(sampled_entry.iloc[0])
			print('\n')
			print('Pathology report : ')
			print(sampled_entry.Report_Text.values[0])
			print('\n')

			continue_eval = input('Would you like to continue? (Y/N) : ') == 'Y'
			count += 1
			if report_numbers_to_query is not None and count == len(report_numbers_to_query):
				print('Reached the end of list...')
				continue_eval = False
		return


		