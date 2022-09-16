"""
Oct. 17, 2021

BiopsyGradeExtractor
- fit
	display summary
- validate result
	user oriented evaluation of the algorithm

TODO :
- pT stage extractor from Radical prostetectomy
- PSA level extractor from progress notes (potentially)

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
		directional_words = ['right', 'left']
		# declare lists which store abstracted biopsy features
		primary_grade_list = []; secondary_grade_list = []; overall_grade_list = []; overall_grade_merged_list = []; overall_gs_list = []; benign_list = []
		num_pos_cores_list = []; num_pos_cores_sum_list = []; num_total_core_list = []; num_total_core_sum_list = []; max_core_involve_list = [];
		small_cell_carc = []; neuroendocrine_carc = []; adenocarcinoma = []
		stage_count = 0; margin_counter = 0; margin_counter_hard_coded = 0; pt_stage_counter_hard_coded = 0

		# create mapping between numerical numbers and numbers in string: used for number of core involvement 		
		num_words_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
		desc_num_words_list = ['one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine', 'ten']
		desc_num_to_num_dic = {desc_num:num for num, desc_num in zip(num_words_list, desc_num_words_list)}

		# iterate through reports
		for text_oi, report_number in tqdm(zip(self.df_pathology_with_biopsy.Report_Text.values, self.df_pathology_with_biopsy.Report_Number.values), total = len(self.df_pathology_with_biopsy), desc = 'extracing info...'):
			# create lists which store abstracted biopsy features locally. "Local grade" correpsonds to grade of each tissue 
			primary_grade_list_local = []; secondary_grade_list_local = []; overall_grade_list_local = []; overall_gs_list_local = []
			num_pos_core_list_local = []; num_total_core_list_local = []; max_core_involve_list_local = [];
			
			# make all words in lowercase for consistency
			text_oi = text_oi.lower()

			# check if it's benign tumor; needs more refinement
			if 'gleason score' in text_oi or 'gleason grade' in text_oi:
				benign_list.append(0)
			else:
				if ('no tumor seen' in text_oi or 'no carcinoma' in text_oi):
					benign_list.append(1)
				else:
					benign_list.append(None)

			# create a list of all the words in each report
			text_oi_list = text_oi.split(' ')
			found_directional_word = 0; 
			# declare indicator variables for histology
			# small_cell_carc_ind = 0; neuroendocrine_carc_ind = 0; adenocarcinoma_ind = 0;
			# iterate through words
			for idx, val in enumerate(text_oi_list):
				if any(dir_val in val for dir_val in directional_words):
					found_directional_word = 1
				# get histology
				# small cell carcinoma
				# if 'small' in val and 'cell' in text_oi_list[idx+1]:
				# 	if 'carcinoma' in [txt for txt in text_oi_list[idx + 1:idx + 5]]:
				# 		text_detect_neg = text_oi_list[idx - 3:idx]
				# 		det_negation = self._func_det_negation(text_detect_neg)
				# 		if not det_negation:
				# 			small_cell_carc_ind = 1
				# # neuroendocrine carcinoma
				# if 'neuroendocrine' in val and 'carcinoma' in text_oi_list[idx + 1]:
				# 	text_detect_neg = text_oi_list[idx - 3:idx]
				# 	det_negation = self._func_det_negation(text_detect_neg)
				# 	if not det_negation:
				# 		neuroendocrine_carc_ind = 1
				# # adenocarcinoma
				# if 'adenocarcinoma' in val:
				# 	text_detect_neg = text_oi_list[idx - 3:idx]
				# 	det_negation = self._func_det_negation(text_detect_neg)
				# 	if not det_negation:
				# 		adenocarcinoma_ind = 1

				primary_grade = None; secondary_grade = None; num_pos_core = None; num_total_core = None; max_core_involve = [];
				# Obtain overall grade (gleason), each biopsy grade, and number of positive cores / total cores 
				# , which appear after words like 'gleason' and directional words in each note. Note that we need to exclude false positive triggered
				# by M.D. Gleason
				if 'gleason' in val and found_directional_word and 'm.d.,' not in text_oi_list[idx - 2 : idx + 2]:
					found_primary_grade = 0; found_secondary_grade = 0; found_plus_sign = 0; found_one_dash = 0
					# get context around word "gleason"
					# heuristically choose 15. Gleason score we're after should be included in the next 15 words（might need to change this)
					raw_gleason_text = ''.join(text_oi_list[idx:idx + 15]) 
					
					# heuristically choose 25. The abstraction info we're after should be included in the next 25 words
					raw_gleason_cores_contex_joined = ' '.join(text_oi_list[idx:idx + 25])
					raw_gleason_cores_context  = text_oi_list[idx:idx + 25] 
					# Might try without the code.
					# get number of positive and total cores
					if 'of' in raw_gleason_cores_context:
						for idx_core_word, core_word in enumerate(raw_gleason_cores_context):
							if core_word == 'of':
								context_core_oi = raw_gleason_cores_context[idx_core_word-1:idx_core_word+2]
								# if successfully extracted three words
								if len(context_core_oi) == 3:
									# remove paranthesis
									pos_core_word_oi = re.sub('[()]', '', context_core_oi[0])
									total_core_word_oi = re.sub('[()]', '', context_core_oi[2])								
									# involving num of num
									if pos_core_word_oi.isdigit() and total_core_word_oi.isdigit():
										num_pos_core = int(pos_core_word_oi)
										num_total_core = int(total_core_word_oi)
									# involving num in str of num in str
									elif pos_core_word_oi in desc_num_words_list and total_core_word_oi in desc_num_words_list:
										num_pos_core = desc_num_to_num_dic[pos_core_word_oi]
										num_total_core = desc_num_to_num_dic[total_core_word_oi]
									# involving num of num in str
									elif pos_core_word_oi.isdigit() and total_core_word_oi in desc_num_words_list:
										num_pos_core = int(pos_core_word_oi)
										num_total_core = desc_num_to_num_dic[total_core_word_oi]
									# or vice versa
									elif pos_core_word_oi in desc_num_words_list and total_core_word_oi.isdigit():
										num_pos_core = desc_num_to_num_dic[pos_core_word_oi]
										num_total_core = int(total_core_word_oi)
									if type(num_total_core) == int:
										break

					# get maximum core involvement for each positive core
					if '%' in raw_gleason_cores_contex_joined and num_pos_core is not None:
						for idx_mci, char_mci in enumerate(raw_gleason_cores_contex_joined):
							if char_mci == '%' and len(max_core_involve) < num_pos_core:
								percent_word = raw_gleason_cores_contex_joined[idx_mci-3:idx_mci+1]
								# when 5-10%, get 10%
								percent_word = percent_word.split('-')[-1]
								# take only numeric parts
								percent_numeric = re.sub('\D', '', percent_word)
								if percent_numeric.isdigit():
									if percent_numeric == '00':
										max_core_involve.append(100)
									else:
										max_core_involve.append(int(percent_numeric))

					# remove text in () or []: e.g. Score 3 (60%) + 4 (40%) = 7/10 -> Score 3 + 4 = 7/10
					raw_gleason_text = re.sub("[\(\[].*?[\)\]]", "", raw_gleason_text)
					# to capture gleason 3/5 where primary 3 and secondary 3, check if there's only one slash without any + or - signs. 
					slash_counter = 0; dash_counter = 0; plus_sign_counter = 0; only_one_slash = 0
					# get context 
					raw_gleason_text_slash = ''.join(text_oi_list[idx:idx + 15]) 
					for char_dash in raw_gleason_text_slash:
						if '/' in char_dash:
							slash_counter += 1
						if '+' in char_dash:
							plus_sign_counter += 1
						if '-' in char_dash:
							dash_counter += 1

					if slash_counter == 1 and plus_sign_counter == 0 and dash_counter == 0:
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

						# vi) when it is represented as z (x + y), where z is the total gleason score, x is the primary grade, and y is the secondary grade
						if char.isdigit() and idx_char < len(raw_gleason_text) - 5:
							if raw_gleason_text[idx_char + 1] == '(':								
								if raw_gleason_text[idx_char + 2].isdigit() and raw_gleason_text[idx_char + 3] == '+' and raw_gleason_text[idx_char + 4].isdigit():
									primary_grade = int(raw_gleason_text[idx_char + 2])
									secondary_grade = int(raw_gleason_text[idx_char + 4])

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

				if num_pos_core is not None:
					num_pos_core_list_local.append(num_pos_core)
					num_total_core_list_local.append(num_total_core)

				if len(max_core_involve) > 0:# is not None:
					max_core_involve_list_local.append(max_core_involve)
				
			# store abstracted biopsy features
			# store local (i.e. each tissue) primary grades, secondary grades, overall grade groups
			if len(primary_grade_list_local) > 0:
				primary_grade_list.append(primary_grade_list_local)
				secondary_grade_list.append(secondary_grade_list_local)
				overall_grade_list.append(overall_grade_list_local)	
				# the maximum of local grades is the overall grade group, 
				if len(set(overall_grade_list_local)) > 0:
					overall_grade_merged_list.append(np.max(overall_grade_list_local))
				else:
					overall_grade_merged_list.append(None)
				overall_gs_list.append(overall_gs_list_local)
			else:
				primary_grade_list.append(None)
				secondary_grade_list.append(None)
				overall_grade_list.append(None)	
				overall_gs_list.append(None)
				overall_grade_merged_list.append(None)

			# store number of positive and total cores
			if len(num_pos_core_list_local) > 0:
				num_pos_cores_list.append(num_pos_core_list_local)
				num_pos_cores_sum_list.append(np.sum(num_pos_core_list_local))
				num_total_core_list.append(num_total_core_list_local)
				num_total_core_sum_list.append(np.sum(num_total_core_list_local))
			else:
				num_pos_cores_list.append(None)
				num_total_core_list.append(None)
				num_pos_cores_sum_list.append(None)
				num_total_core_sum_list.append(None)

			# store max core involvement
			if len(max_core_involve_list_local) > 0:
				max_core_involve_list.append(max_core_involve_list_local)
			else:
				max_core_involve_list.append(None)

			# store histology
			# if adenocarcinoma_ind:
			# 	adenocarcinoma.append(1)
			# else:
			# 	adenocarcinoma.append(0)

			# if small_cell_carc_ind and not adenocarcinoma_ind:
			# 	small_cell_carc.append(1)
			# else:
			# 	small_cell_carc.append(0)
				
			# if neuroendocrine_carc_ind and not adenocarcinoma_ind:
			# 	neuroendocrine_carc.append(1)
			# else:
			# 	neuroendocrine_carc.append(0)

		# store abstracted biopsy features into the main df
		# self.df_pathology_with_biopsy.loc[self.df_pathology_with_biopsy.index, 'num_pos_cores'] = num_pos_cores_list
		self.df_pathology_with_biopsy.loc[self.df_pathology_with_biopsy.index, 'primary_grade'] = primary_grade_list
		self.df_pathology_with_biopsy.loc[self.df_pathology_with_biopsy.index, 'secondary_grade'] = secondary_grade_list
		# self.df_pathology_with_biopsy.loc[self.df_pathology_with_biopsy.index, 'overall_grade_group'] = overall_grade_list
		# self.df_pathology_with_biopsy.loc[self.df_pathology_with_biopsy.index, 'overall_gs'] = overall_gs_list
		self.df_pathology_with_biopsy.loc[self.df_pathology_with_biopsy.index, 'overall_grade_merged'] = overall_grade_merged_list

		# self.df_pathology_with_biopsy.loc[self.df_pathology_with_biopsy.index, 'benign'] = benign_list

		# self.df_pathology_with_biopsy.loc[self.df_pathology_with_biopsy.index, 'num_pos_cores'] = num_pos_cores_list
		# self.df_pathology_with_biopsy.loc[self.df_pathology_with_biopsy.index, 'num_total_core'] = num_total_core_list
		# self.df_pathology_with_biopsy.loc[self.df_pathology_with_biopsy.index, 'num_pos_cores_sum'] = num_pos_cores_sum_list
		# self.df_pathology_with_biopsy.loc[self.df_pathology_with_biopsy.index, 'num_total_core_sum'] = num_total_core_sum_list

		self.df_pathology_with_biopsy.loc[self.df_pathology_with_biopsy.index, 'max_core_involve'] = max_core_involve_list
		# self.df_pathology_with_biopsy.loc[self.df_pathology_with_biopsy.index, 'small_cell_carc'] = small_cell_carc
		# self.df_pathology_with_biopsy.loc[self.df_pathology_with_biopsy.index, 'neuroendocrine_carc'] = neuroendocrine_carc
		# self.df_pathology_with_biopsy.loc[self.df_pathology_with_biopsy.index, 'adenocarcinoma'] = adenocarcinoma		

		print('\n')
		print('Extraction complete.')
		print('\n')
		return

	def fit(self, df_pathology, display_summary = False, return_summary_df = False):
		"""
		inputs :

		"""
		# get all the reports with the following keywords : needle, core, biops, prost
		idx_oi = [idx for idx, report in enumerate(df_pathology.Report_Text.values) if 'needle' in report.lower() or 'core' in report.lower() or 'biops' in report.lower() and 'prost' in report.lower()]
		# exclude all the reports with the following keywords : prostatectomy, bone marrow
		idx_to_exclude = [idx for idx, report in enumerate(df_pathology.Report_Text.values) if 'prostatectomy' in report.lower() or 'bone marrow' in report.lower(	)]# or 'nephrectomy' in report.lower() or 'cystectomy' in report.lower() or 'cystoprostatectomy' in report.lower() or 'ureter' in report.lower()] # or 'prostate, radical resection' in report.lower()]
		self.df_pathology_with_biopsy = df_pathology.iloc[np.sort(list(set(idx_oi) - set(idx_to_exclude)))].copy()	
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

			# Mmaximum core involvement (MCI) outlier filtering
			max_core_involve_all_list = []; outlier_report_list = []
			mci_outliers = [890, 0, 800, 320]
			for val, report_number_mci in zip(self.df_pathology_with_biopsy.max_core_involve.values, self.df_pathology_with_biopsy.Report_Number.values):
				if val is not None:
					for val_sub in val:
						if any([val_outlier in val_sub for val_outlier in mci_outliers]):
							outlier_report_list.append(report_number_mci)
						else:
							max_core_involve_all_list += val_sub
			df_count_stats_max_core_involve = pd.value_counts(max_core_involve_all_list, sort = True)		
			# print('Counts per each max core involvement : ')
			# print(df_count_stats_max_core_involve)
			# print('\n')
			# remove outliers
			reports_mci_oi = set(self.df_pathology_with_biopsy.Report_Number.values) - set(outlier_report_list)
			self.df_pathology_with_biopsy = self.df_pathology_with_biopsy.loc[self.df_pathology_with_biopsy.Report_Number.isin(reports_mci_oi)]

			# benign stats
			# df_count_stats_benign = pd.value_counts(self.df_pathology_with_biopsy.benign.values, sort = True)
			# print('Number of reports with benign indicator 1 (postiive) or 0 (negative) : ')
			# print(df_count_stats_benign)
			# print('total = ', sum(df_count_stats_benign.values))
			# print('\n')
			
			# overall grade stats
			df_count_stats_overall_grade = pd.value_counts(self.df_pathology_with_biopsy.overall_grade_merged.values, sort = True)
			# df_count_stats_num_pos_cores = pd.value_counts(self.df_pathology_with_biopsy.num_pos_cores_sum.values, sort = True)
			# outlier filtering
			# num_pos_sum_outliers = [20, 30]
			# num_pos_cores_oi = set(df_count_stats_num_pos_cores.loc[df_count_stats_num_pos_cores.values >= 9].index) | set([None]) - set(num_pos_sum_outliers)
			# self.df_pathology_with_biopsy = self.df_pathology_with_biopsy.loc[self.df_pathology_with_biopsy.num_pos_cores_sum.isin(num_pos_cores_oi)]

			# df_count_stats_num_total_cores = pd.value_counts(self.df_pathology_with_biopsy.num_total_core_sum.values, sort = True)
			# # outlier filtering
			# num_total_core_sum_outliers = [30]
			# num_total_cores_oi = set(df_count_stats_num_total_cores.loc[df_count_stats_num_total_cores.values >= 9].index) | set([None]) - set(num_total_core_sum_outliers)
			# self.df_pathology_with_biopsy = self.df_pathology_with_biopsy.loc[self.df_pathology_with_biopsy.num_total_core_sum.isin(num_total_cores_oi)]

			# print stats
			# df_count_stats_num_pos_cores = pd.value_counts(self.df_pathology_with_biopsy.num_pos_cores_sum.values, sort = True)
			# df_count_stats_num_total_cores = pd.value_counts(self.df_pathology_with_biopsy.num_total_core_sum.values, sort = True)
			# print('Number of reports per each # of positive cores : ')
			# print(df_count_stats_num_pos_cores)
			# print('Number of reports per each # of total cores : ')
			# print(df_count_stats_num_total_cores)
			print('Number of reports per each overall grade : ')
			print(df_count_stats_overall_grade)
			print('total = ', df_count_stats_overall_grade.sum())
			print('Number of unique biopsy reports : ', len(self.df_pathology_with_biopsy))
			# print('Number of unique biopsy reports per center : ')
			# print(pd.value_counts(self.df_pathology_with_biopsy.MRN_Type.values, sort = True))
			print('Number of unique patients : ', len(set(self.df_pathology_with_biopsy.EMPI.values)))

			df_pathology_with_biopsy_oi_overall_grade = self.df_pathology_with_biopsy.dropna(subset = ['overall_grade_merged'])
			print('Number of unique patients with overall_grade_merged : ', len(set(df_pathology_with_biopsy_oi_overall_grade.EMPI.values)))
			# df_pathology_with_biopsy_oi_with_num_cores_involved = df_pathology_with_biopsy_oi_overall_grade.dropna(subset = ['num_pos_cores'])
			# df_pathology_with_biopsy_oi_with_num_cores_involved_wo_grade = self.df_pathology_with_biopsy.dropna(subset = ['num_pos_cores'])
			# print('Number of unique patients with overall grade, and num. of positive and total cores : ', len(set(df_pathology_with_biopsy_oi_with_num_cores_involved.EMPI.values)))
			# print('Number of unique patients with num. of positive and total cores : ', len(set(df_pathology_with_biopsy_oi_with_num_cores_involved_wo_grade.EMPI.values)))


			# df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve = df_pathology_with_biopsy_oi_with_num_cores_involved.dropna(subset = ['max_core_involve'])
			# print('Number of unique patients with overall grade, num. of positive and total cores, and maximum core involvement : ', len(set(df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve.EMPI.values)))
			# df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve_benign = df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve.dropna(subset = ['benign'])
			# print('Number of unique patients with overall grade, num. of positive and total cores, maximum core involvement, and benign : ', len(set(df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve_benign.EMPI.values)))
			# print('\n')
			# print('\n')

		if return_summary_df:
			return self.df_pathology_with_biopsy.drop('max_core_involve',axis = 1)
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


		