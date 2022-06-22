import pandas as pd
import numpy as np
from tqdm import tqdm
import os
import re

"""
Intae Moon
Date : May 3rd, 2021

Description :
The following script extracts:

1. Overall grade group
2. Benign tumor indicator
3. Number of positive cores
4. Total number of cores
5. Maximum core involvement

from biopsy reports. 

See https://docs.google.com/document/d/1aUFjvz8bumhCnSUDJ8w4CtKoNmGBKYTGeq0cMqoOOl0/edit for more details
"""

def featurize_biopsy_df(df_biopsy):
	# get the following features 
	# 1. mci auxiiliary feature
	df_biopsy_wo_duplicate = df_biopsy.drop_duplicates(subset = 'EMPI', keep = False)
	df_biopsy_wo_duplicate.set_index('EMPI', inplace = True)
	empis_oi = np.unique(df_biopsy_wo_duplicate.index)
	counter = 0
	for empi in empis_oi:
		num_pos_cores_oi = df_biopsy_wo_duplicate.at[empi, 'num_pos_cores']
		overall_grade_group_oi = df_biopsy_wo_duplicate.at[empi, 'overall_grade_group']
		mci_oi = df_biopsy_wo_duplicate.at[empi, 'max_core_involve']
		if len(num_pos_cores_oi) == len(overall_grade_group_oi) == len(mci_oi):
			weighted_val = 0
			for num_pos_core, overall_grade, mci in zip(num_pos_cores_oi, overall_grade_group_oi, mci_oi):
				if len(mci) == num_pos_core: # mci = [a,b,c], num_pos_core = 3, result = a + b + c
					for mci_sub in mci:
						weighted_val += 1 * overall_grade * mci_sub/100 # 
				else: 
					# mci = [a,b], num_pos_core = 3, result = a + b + b or mci = [a,b,c,d], num_pos_core = 3, result = a + b + c
					# if len(mci) == 1:
					# 	weighted_val += 1 * overall_grade * mci_sub/100/num_pos_core # mci = [a], 
					# else:
					for i in range(num_pos_core):
						try:
							weighted_val += 1 * overall_grade * mci[i]/100
						except:
							weighted_val += 1 * overall_grade * mci[-1]/100
			df_biopsy_wo_duplicate.at[empi, 'auxiiliary_mci_score'] = weighted_val
		else:
			counter += 1
	df_biopsy_wo_duplicate.dropna(subset = ['auxiiliary_mci_score'], inplace = True)
	return df_biopsy_wo_duplicate

# simple function to detect negation
def func_det_negation(text_detect_neg):
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
	
def main():
	load_prcossed = input('Load the processed data? (Y/N) : ') == 'Y'
	if load_prcossed:
		df_pathology_with_biopsy = pd.read_csv('df_pathology_with_biopsy_params_extracted.tsv', sep = '\t')
		df_pathology_with_biopsy.set_index(df_pathology_with_biopsy.columns[0], inplace = True)
	else:
		# load raw pathology notes 
		# df_pathology = pd.read_csv('../data/Pat_total.csv', sep="|", low_memory=False)
		# df_pathology = df_pathology.loc[df_pathology.Report_Description == 'Surgical Pathology']
		df_pathology = pd.read_csv('../data/sample_biopsy_reports.csv')
		# df_pathology = df_pathology.iloc[[3]]
		# get all the reports with the following keywords : needle, core, biops, prost
		idx_oi = [idx for idx, report in enumerate(df_pathology.Report_Text.values) if 'needle' in report.lower() or 'core' in report.lower() or 'biops' in report.lower() and 'prost' in report.lower()]
		# exclude all the reports with the following keywords : prostatectomy, bone marrow
		idx_to_exclude = [idx for idx, report in enumerate(df_pathology.Report_Text.values) if 'prostatectomy' in report.lower()]# or 'nephrectomy' in report.lower() or 'cystectomy' in report.lower() or 'cystoprostatectomy' in report.lower() or 'ureter' in report.lower()] # or 'prostate, radical resection' in report.lower()]
		df_pathology_with_biopsy = df_pathology.iloc[np.sort(list(set(idx_oi) - set(idx_to_exclude)))].copy()
		
		# location of biopsy site within prostate
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
		for text_oi, report_number in tqdm(zip(df_pathology_with_biopsy.Report_Text.values, df_pathology_with_biopsy.Report_Number.values), total = len(df_pathology_with_biopsy), desc = 'extracing info...'):
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
			small_cell_carc_ind = 0; neuroendocrine_carc_ind = 0; adenocarcinoma_ind = 0;
			# iterate through words
			for idx, val in enumerate(text_oi_list):
				if any(dir_val in val for dir_val in directional_words):
					found_directional_word = 1
				# get histology : 
				# small cell carcinoma
				if 'small' in val and 'cell' in text_oi_list[idx+1]:
					if 'carcinoma' in [txt for txt in text_oi_list[idx + 1:idx + 5]]:
						text_detect_neg = text_oi_list[idx - 3:idx]
						det_negation = func_det_negation(text_detect_neg)
						if not det_negation:
							small_cell_carc_ind = 1
				# neuroendocrine carcinoma
				if 'neuroendocrine' in val and 'carcinoma' in text_oi_list[idx + 1]:
					text_detect_neg = text_oi_list[idx - 3:idx]
					det_negation = func_det_negation(text_detect_neg)
					if not det_negation:
						neuroendocrine_carc_ind = 1
				# adenocarcinoma
				if 'adenocarcinoma' in val:
					text_detect_neg = text_oi_list[idx - 3:idx]
					det_negation = func_det_negation(text_detect_neg)
					if not det_negation:
						adenocarcinoma_ind = 1

				primary_grade = None; secondary_grade = None; num_pos_core = None; num_total_core = None; max_core_involve = [];
				# Obtain overall grade (gleason), each biopsy grade, and number of positive cores / total cores 
				# , which appear after words like 'gleason' and directional words in each note. Note that we need to exclude false positive triggered
				# by M.D. Gleason
				if 'gleason' in val and 'm.d.,' not in text_oi_list[idx - 2 : idx + 2]:
					found_primary_grade = 0; found_secondary_grade = 0; found_plus_sign = 0; found_minus_sign = 0; found_one_dash = 0
					# get context around word "gleason"
					# heuristically choose 15. Gleason score we're after should be included in the next 15 words
					raw_gleason_text = ''.join(text_oi_list[idx:idx + 15]) 
					
					# heuristically choose 25. The abstraction info we're after should be included in the next 25 words
					raw_gleason_cores_contex_joined = ' '.join(text_oi_list[idx:idx + 25])
					raw_gleason_cores_context  = text_oi_list[idx:idx + 25]
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


					# check for groups in parens first.
					raw_gleason_paren_search = re.findall("(\([1-9]([+]|and)[1-9][^)]*\))", raw_gleason_text)
					if len(raw_gleason_paren_search) != 0:
						for i in raw_gleason_paren_search:
							raw_gleason_paren = i[0]
							#raw_gleason_paren = raw_gleason_paren_match[0]
							digits_in_paren_statement = ''.join(c for c in raw_gleason_paren if c.isdigit())
							if len(digits_in_paren_statement) == 2:
								primary_grade = int(digits_in_paren_statement[0])
								primary_grade_list_local.append(primary_grade)
								found_primary_grade = 1
								secondary_grade = int(digits_in_paren_statement[1])
								secondary_grade_list_local.append(secondary_grade)
								found_secondary_grade = 1
								if primary_grade is not None and secondary_grade is not None:
									found_directional_word = 0  # set this to 0 so it can search other Gleason scores in the report

									# had to include a separate overall grade group calculator hereq

									if primary_grade + secondary_grade <= 6:
										overall_grade_group = 1
									elif primary_grade + secondary_grade == 7:
										if found_primary_grade and found_secondary_grade:
											if primary_grade == 3:
												overall_grade_group = 2
											elif primary_grade == 4:
												overall_grade_group = 3
										else:  # 7/10 with no primary and seoncdary info
											overall_grade_group = 2.5  ## from the doc - if can't distinguish primary and secondary: divide by 2
									elif primary_grade + secondary_grade == 8:
										overall_grade_group = 4
									elif primary_grade + secondary_grade > 8:
										overall_grade_group = 5
									else:
										overall_grade_group = None
								if overall_grade_group is not None:
									overall_grade_list_local.append(overall_grade_group)
									overall_gs_list_local.append(primary_grade + secondary_grade)
								primary_grade = None


								secondary_grade = None



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
							raw_gleason_equation = raw_gleason_text[max(0,idx_char-10):idx_char+10]
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
									found_primary_grade = 1
									found_secondary_grade = 1
									break
							# iii) when overall grade is represented as Gleason 6/10 
							raw_gleason_grade_dash_and_ten = raw_gleason_text[idx_char-1:idx_char+3]
							raw_gleason_grade_dash_and_ten_split = raw_gleason_grade_dash_and_ten.split('/')
							if raw_gleason_grade_dash_and_ten_split[0].isdigit() and raw_gleason_grade_dash_and_ten_split[1] == '10':
								# in this case primary grade = secondary grade = 3
								primary_grade = int(raw_gleason_grade_dash_and_ten_split[0])/2
								secondary_grade = primary_grade
								found_primary_grade = 0
								found_secondary_grade = 0
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
								found_primary_grade = 1
								found_secondary_grade = 1
								break
						# v) when it is represented as gleason 3/5 (where primary 3 and secondary 3)
						elif found_one_dash and only_one_slash:
							primary_grade = primary_grade_to_be
							secondary_grade = primary_grade_to_be
							found_primary_grade = 1
							found_secondary_grade = 1
							break


						# vii) (mike and gha) gleason grade x-y
						if char == '-' and slash_counter == 0:
							raw_gleason_equation = raw_gleason_text[idx_char - 10:idx_char + 10]
							for sub_char in raw_gleason_equation:
								# checking - sign again; it covers a few corner cases
								if sub_char == '-':
									found_minus_sign = 1
								if sub_char.isdigit():
									if not found_minus_sign:
										primary_grade_to_be = int(sub_char)
										found_primary_grade = 1
									elif found_primary_grade and found_minus_sign and not found_secondary_grade:
										primary_grade = primary_grade_to_be
										secondary_grade = int(sub_char)
										found_secondary_grade = 1
										break
							break



					# viii) (mike and gha) gleason grade x and y; x of 5 and y of 5; x and y of 5;
					# a little redundant but can also search for /5 in addition ot "of 5"
					# should find "and" and characters around, grab any digits that aren't "of 5" or /5 and return those as primary and secondary grade
					if found_primary_grade == 0 and found_secondary_grade == 0:
						found_and = raw_gleason_text.find("and")
						found_of_ten = raw_gleason_text.find("of10")
						found_of_five = raw_gleason_text.find("of5")
						if found_and != -1:
							raw_gleason_equation = raw_gleason_text[found_and - 10:found_and + 10]
							found_core = raw_gleason_equation.find("cores")
							if found_core == -1:
								raw_gleason_equation = re.sub('of5|/5', '', raw_gleason_equation)
								digits_in_and_statement = ''.join(
									c for c in raw_gleason_equation if c.isdigit())
								if len(digits_in_and_statement) == 2:
									primary_grade = int(digits_in_and_statement[0])
									found_primary_grade = 1
									secondary_grade = int(digits_in_and_statement[1])
									found_secondary_grade = 1
								elif len(digits_in_and_statement) < 2:
									found_primary_grade = 0
									found_secondary_grade = 0
								else:
									found_primary_grade = 0
									found_secondary_grade = 0
							else:
								found_primary_grade == 0
								found_secondary_grade = 0
						if found_of_five != -1:
							if raw_gleason_text[found_of_five - 1].isdigit():
								primary_grade = int(raw_gleason_text[found_of_five - 1])
								secondary_grade = primary_grade
								found_primary_grade = 1
								found_secondary_grade = 1
						elif found_of_ten != -1:
							if raw_gleason_text[found_of_ten - 1].isdigit():
								total_score = int(raw_gleason_text[found_of_ten - 1])
								primary_grade = total_score / 2
								secondary_grade = primary_grade
								found_primary_grade = 0
								found_secondary_grade = 0





						#### can check this w Intae/test it and maybe add code to handle cases where we pick up more than two numbers here


				# get overall garde group based on abstracted gleason scores
				# https://www.pcf.org/about-prostate-cancer/diagnosis-staging-prostate-cancer/gleason-score-isup-grade/
				if primary_grade is not None and secondary_grade is not None:
					found_directional_word = 0 # set this to 0 so it can search other Gleason scores in the report
					if primary_grade + secondary_grade <= 6:
						overall_grade_group = 1
					elif primary_grade + secondary_grade == 7:
						if found_primary_grade and found_secondary_grade:
							if primary_grade == 3:
								overall_grade_group = 2
							elif primary_grade == 4:
								overall_grade_group = 3
						else: # 7/10 with no primary and seoncdary info
							overall_grade_group = 2.5 ## from the doc - if can't distinguish primary and secondary: divide by 2
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
			if adenocarcinoma_ind:
				adenocarcinoma.append(1)
			else:
				adenocarcinoma.append(0)

			if small_cell_carc_ind and not adenocarcinoma_ind:
				small_cell_carc.append(1)
			else:
				small_cell_carc.append(0)
				
			if neuroendocrine_carc_ind and not adenocarcinoma_ind:
				neuroendocrine_carc.append(1)
			else:
				neuroendocrine_carc.append(0)

		# store abstracted biopsy features into the main df
		df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'num_pos_cores'] = num_pos_cores_list
		df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'primary_grade'] = primary_grade_list
		df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'secondary_grade'] = secondary_grade_list
		df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'overall_grade_group'] = overall_grade_list
		df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'overall_gs'] = overall_gs_list
		df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'overall_grade_merged'] = overall_grade_merged_list

		df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'benign'] = benign_list

		df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'num_pos_cores'] = num_pos_cores_list
		df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'num_total_core'] = num_total_core_list
		df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'num_pos_cores_sum'] = num_pos_cores_sum_list
		df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'num_total_core_sum'] = num_total_core_sum_list

		df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'max_core_involve'] = max_core_involve_list
		df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'small_cell_carc'] = small_cell_carc
		df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'neuroendocrine_carc'] = neuroendocrine_carc
		df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'adenocarcinoma'] = adenocarcinoma

	if not load_prcossed:
		primary_grade_all_list = []
		for val in df_pathology_with_biopsy.primary_grade.values:
			if val is not None:
				primary_grade_all_list += val
		df_count_stats_primary_grade = pd.value_counts(primary_grade_all_list, sort = True)		
		print('Counts per each primary grade: ')
		print(df_count_stats_primary_grade)
		print('\n')

		secondary_grade_all_list = []
		for val in df_pathology_with_biopsy.secondary_grade.values:
			if val is not None:
				secondary_grade_all_list += val
		df_count_stats_secondary_grade = pd.value_counts(secondary_grade_all_list, sort = True)		
		print('Counts per each secondary grade: ')
		print(df_count_stats_secondary_grade)
		print('\n')

		# Mmaximum core involvement (MCI) outlier filtering
		max_core_involve_all_list = []; outlier_report_list = []
		mci_outliers = [890, 0, 800, 320]
		for val, report_number_mci in zip(df_pathology_with_biopsy.max_core_involve.values, df_pathology_with_biopsy.Report_Number.values):
			if val is not None:
				for val_sub in val:
					if any([val_outlier in val_sub for val_outlier in mci_outliers]):
						outlier_report_list.append(report_number_mci)
					else:
						max_core_involve_all_list += val_sub
		df_count_stats_max_core_involve = pd.value_counts(max_core_involve_all_list, sort = True)		
		print('Counts per each max core involvement : ')
		print(df_count_stats_max_core_involve)
		print('\n')
		# remove outliers
		reports_mci_oi = set(df_pathology_with_biopsy.Report_Number.values) - set(outlier_report_list)
		df_pathology_with_biopsy = df_pathology_with_biopsy.loc[df_pathology_with_biopsy.Report_Number.isin(reports_mci_oi)]

		# benign stats
		df_count_stats_benign = pd.value_counts(df_pathology_with_biopsy.benign.values, sort = True)
		print('Number of reports with benign indicator 1 (postiive) or 0 (negative) : ')
		print(df_count_stats_benign)
		print('total = ', sum(df_count_stats_benign.values))
		print('\n')
		
		# overall grade stats
		df_count_stats_overall_grade = pd.value_counts(df_pathology_with_biopsy.overall_grade_merged.values, sort = True)
		df_count_stats_num_pos_cores = pd.value_counts(df_pathology_with_biopsy.num_pos_cores_sum.values, sort = True)
		# outlier filtering
		# num_pos_sum_outliers = [20, 30]
		# num_pos_cores_oi = set(df_count_stats_num_pos_cores.loc[df_count_stats_num_pos_cores.values >= 9].index) | set([None]) - set(num_pos_sum_outliers)
		# df_pathology_with_biopsy = df_pathology_with_biopsy.loc[df_pathology_with_biopsy.num_pos_cores_sum.isin(num_pos_cores_oi)]

		# df_count_stats_num_total_cores = pd.value_counts(df_pathology_with_biopsy.num_total_core_sum.values, sort = True)
		# # outlier filtering
		# num_total_core_sum_outliers = [30]
		# num_total_cores_oi = set(df_count_stats_num_total_cores.loc[df_count_stats_num_total_cores.values >= 9].index) | set([None]) - set(num_total_core_sum_outliers)
		# df_pathology_with_biopsy = df_pathology_with_biopsy.loc[df_pathology_with_biopsy.num_total_core_sum.isin(num_total_cores_oi)]

		# print stats
		df_count_stats_num_pos_cores = pd.value_counts(df_pathology_with_biopsy.num_pos_cores_sum.values, sort = True)
		df_count_stats_num_total_cores = pd.value_counts(df_pathology_with_biopsy.num_total_core_sum.values, sort = True)
		print('Number of reports per each # of positive cores : ')
		print(df_count_stats_num_pos_cores)
		print('Number of reports per each # of total cores : ')
		print(df_count_stats_num_total_cores)
		print('Number of reports per each overall grade : ')
		print(df_count_stats_overall_grade)
		print('total = ', df_count_stats_overall_grade.sum())
		print('Number of unique biopsy reports : ', len(df_pathology_with_biopsy))
		print('Number of unique biopsy reports per center : ')
		print(pd.value_counts(df_pathology_with_biopsy.MRN_Type.values, sort = True))
		print('Number of unique patients : ', len(set(df_pathology_with_biopsy.EMPI.values)))

		df_pathology_with_biopsy_oi_overall_grade = df_pathology_with_biopsy.dropna(subset = ['overall_grade_group'])
		print('Number of unique patients with overall grade : ', len(set(df_pathology_with_biopsy_oi_overall_grade.EMPI.values)))
		df_pathology_with_biopsy_oi_with_num_cores_involved = df_pathology_with_biopsy_oi_overall_grade.dropna(subset = ['num_pos_cores'])
		df_pathology_with_biopsy_oi_with_num_cores_involved_wo_grade = df_pathology_with_biopsy.dropna(subset = ['num_pos_cores'])
		print('Number of unique patients with overall grade, and num. of positive and total cores : ', len(set(df_pathology_with_biopsy_oi_with_num_cores_involved.EMPI.values)))
		print('Number of unique patients with num. of positive and total cores : ', len(set(df_pathology_with_biopsy_oi_with_num_cores_involved_wo_grade.EMPI.values)))


		df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve = df_pathology_with_biopsy_oi_with_num_cores_involved.dropna(subset = ['max_core_involve'])
		print('Number of unique patients with overall grade, num. of positive and total cores, and maximum core involvement : ', len(set(df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve.EMPI.values)))
		df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve_benign = df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve.dropna(subset = ['benign'])
		print('Number of unique patients with overall grade, num. of positive and total cores, maximum core involvement, and benign : ', len(set(df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve_benign.EMPI.values)))
		print('\n')
		print('\n')
	continue_ = True; sample_rp = True; specific_report = False
	if load_prcossed:
		df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve_benign = df_pathology_with_biopsy
	
	df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve_benign.drop(columns = ['Unnamed: 0'], inplace = True)
	df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve_benign.index = np.arange(len(df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve_benign))
	df_biopsy_wo_duplicate_processed = featurize_biopsy_df(df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve_benign)
	df_biopsy_wo_duplicate_processed.to_csv('../data/df_biopsy_wo_duplicate_processed.csv')
	# breakpoint()
	# validate algorithm precision in a randomly chosen note
	while continue_:
		print('\n')
		print('Extracted info : ')
		if specific_report:
			sampled_entry = df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve_benign.loc[df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve_benign.Report_Number == report_number_to_query]
		else:
			sampled_entry = df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve_benign.sample(n = 1)
		sampled_entry_extracted_oi = sampled_entry[['Report_Number', 'overall_gs', 'overall_grade_group', 'overall_grade_merged', 'num_pos_cores', 'num_pos_cores_sum', 'num_total_core', 'num_total_core_sum', 'max_core_involve', 'benign', 'small_cell_carc', 'neuroendocrine_carc', 'adenocarcinoma']]
		print(sampled_entry_extracted_oi.iloc[0])
		print('\n')
		print('Pathology report : ')
		print(sampled_entry.Report_Text.values[0])
		print('\n')
		continue_ = input('Would you like to continue? (Y/N) : ') == 'Y'
		specific_report = input('Want to look at specific report? (Y/N) : ') == 'Y'
		if specific_report:
			report_number_to_query = input('Type report number : ')
			if report_number_to_query not in df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve_benign.Report_Number.values:
				print('Report not found!')
				print('Randomly sampling a report...')
				specific_report = False
				continue
	
	export = input('Export the pathology report datafame with the extracted info in the current directory? (Y/N) : ')
	if export not in {'Y', 'N'}:
		raise KeyError('Input Y or N only')
	if export == 'Y':
		df_pathology_with_biopsy.to_csv('df_pathology_with_biopsy_params_extracted.tsv', sep = '\t', index = True)
		df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve_benign.to_csv('df_pathology_with_biopsy_params_extracted_all_feats.tsv', sep = '\t', index = True)
		print('\n')
		print('The processed file exported to ' + os.getcwd() + '/df_pathology_with_extracted_info.tsv')
		print('\n')
	return

if __name__ == "__main__":
	explain = """
	Written by Intae Moon
	Date : May 3rd, 2021

	Description :
	The following script extracts:
	
	1. Overall grade group
	2. Benign tumor indicator
	3. Number of positive cores
	4. Total number of cores
	5. Maximum core involvement

	from biopsy reports. 

	See https://docs.google.com/document/d/1aUFjvz8bumhCnSUDJ8w4CtKoNmGBKYTGeq0cMqoOOl0/edit for more details
	"""
	print("\n")
	print(explain)
	print("\n")
	input('Type any button if you would like to continue')
	main()
