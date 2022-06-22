import pandas as pd
import numpy as np
from tqdm import tqdm
import os
import re

"""
Intae Moon
Date : April 3, 2021

Description :

The following script extracts 
1. primary grade
2. secondary grade
3. overall grade group
4. histology
5. pT stage
6. n1 status
7. surgical margin

from RP pathology notes
"""

def get_tumor_stage_info_via_context(df_pathology_with_rp):
	# hard-coded rule based pT stage abstraction:
	tumor_stage_list = []
	hard_coded_pt2_stage = ["Extraprostatic extension of tumor is not identified", "confined to the prostate", "Extraprostatic extension of tumor is not identified", "confined to the gland", "confined to prostate", "confined within the prostatic capsule", "involves the capsule but is confined to the prostate", "invades but does not transgress the capsule", "No capsular penetration demonstrated", "NOT INFILTRATING PERIPROSTATIC ADIPOSE TISSUE", "tumor does not transgress the prostatic capsule", "apparently localized", "invades into but not through the prostatic capsule", "extends to, but not through the prostatic capsule.", "the periprostatic and soft tissue resection margins are free of carcinoma", "no capsular invasion is present.", "extracapsular extension is not identified", "confined within the prostatic capsule", "TUMOR IS NOT SEEN OUTSIDE THE CAPSULE", "does not extend into periprostatic fat", "is confirmed to the prostate", "Extraprostatic extension is not present", "tumor does not transgress the prostatic capsule", "All inked margins and capsule are free of tumor", "no tumor involvement of the capsule", "Extraprostatic extension is not identified", "No extraprostatic extension of tumor"]
	hard_coded_pt2_stage = [val.lower() for val in hard_coded_pt2_stage]

	hard_coded_pt3a_stage = ["extends into the prostatic", "extends through the capsule into periprostatic", "extends slightly into periprostatic tissue", "LIMITED EXTRAPROSTATIC EXTENSION", "tumor extends into extraprostatic soft tissue", "tumor focally penetrates through the prostate capsule", "focally into the periprostatic fat"]
	hard_coded_pt3a_stage = [val.lower() for val in hard_coded_pt3a_stage]

	hard_coded_pt3b_stage = ["extends into extra prostatic soft tissues and the left seminal vesicle", "and involves the left seminal vesicle"]
	hard_coded_pt3b_stage = [val.lower() for val in hard_coded_pt3b_stage]
	
	for text_oi, report_number in tqdm(zip(df_pathology_with_rp.Report_Text.values, df_pathology_with_rp.Report_Number.values), total = len(df_pathology_with_rp), desc = 'extracing pTstage via context...'):
		# lower case all the words in each report
		text_oi_list = text_oi.lower().split(' ')
		tumor_stage_ind = 0; tumor_stage_ind_pt3a = 0; tumor_stage_ind_pt2 = 0; tumor_stage_ind_pt3b = 0

		# check for hard-coded rules (pT stage)
		for val_harcoded_pt2 in hard_coded_pt2_stage:
			if val_harcoded_pt2 in text_oi.lower():
				tumor_stage_ind = 1
				tumor_stage_list.append('pt2')
				break

		for val_harcoded_pt3a in hard_coded_pt3a_stage:
			if val_harcoded_pt3a in text_oi.lower():
				if tumor_stage_ind:
					print('(pT stage) presence of conflicting sentences in Report : ', report_number)
					print(val_harcoded_pt2)
					print(val_harcoded_pt3a)
					print('\n')
					tumor_stage_list.pop()
					tumor_stage_list.append('pt3a')
				else:
					tumor_stage_ind = 1
					tumor_stage_list.append('pt3a')
				break

		for val_harcoded_pt3b in hard_coded_pt3b_stage:
			if val_harcoded_pt3b in text_oi.lower():
				if tumor_stage_ind:
					print('(pT stage) presence of conflicting sentences in Report : ', report_number)
					print(val_harcoded_pt3a)
					print(val_harcoded_pt3b)
					print('\n')
					tumor_stage_list.pop()
					tumor_stage_list.append('pt3b')
				else:
					tumor_stage_ind = 1
					tumor_stage_list.append('pt3b')
				break


		for idx, val in enumerate(text_oi_list):
			# get tumor stage info
			if not tumor_stage_ind:
				if 'confine' in val: # check for "confined to the prostate"
					# check for negation :
					if 'not' not in text_oi_list[idx -1]:
						context_to_look = text_oi_list[idx : idx + 5]
						found_to = 0
						for context_word in context_to_look:
							if 'to' in context_word:
								found_to = 1
							if 'prostate' in context_word and found_to:
								tumor_stage_ind = 1
								tumor_stage_ind_pt2 = 1
								tumor_stage_list.append('pt2')
								break
				elif 'extend' in val:
					idx_extend = idx
					context_to_look = text_oi_list[idx : idx + 15]
					found_to = 0
					for idx, context_word in enumerate(context_to_look):
						if 'to' in context_word: # extends "to" xxx : 
							found_to = 1
						if 'soft' in context_word and found_to:
							if idx < len(context_to_look) - 1:
								if 'tissue' in context_to_look[idx+1]:
									tumor_stage_ind = 1
									tumor_stage_ind_pt3a = 1
									tumor_stage_list.append('pt3a')
									break

					found_to = 0
					for idx, context_word in enumerate(context_to_look):
						if 'to' in context_word: # extends "to" xxx : 
							found_to = 1
						if 'prostatic' in context_word and found_to:
							if idx < len(context_to_look) - 1:
								if 'capsule' in context_to_look[idx+1]:
									if not tumor_stage_ind_pt3a:
										tumor_stage_ind = 1
										tumor_stage_ind_pt3a = 1
										tumor_stage_list.append('pt3a')
									break
					if 'negative' not in context_to_look and 'no' not in context_to_look and 'not' not in context_to_look: #check for negation of seminal vesicle
						found_to = 0
						for idx, context_word in enumerate(context_to_look):
							if 'to' in context_word: # extends "to" xxx : 
								found_to = 1
							if 'seminal' in context_word and found_to:
								if idx < len(context_to_look) - 1:
									if 'vesicle' in context_to_look[idx+1] or 'vesicles' in context_to_look[idx+1]:
										if tumor_stage_ind_pt3a == 1:
											# print('duplicate assignment!')
											# print('fixing it...')
											tumor_stage_list.pop()
											tumor_stage_list.append('pt3b')
											tumor_stage_ind = 1
											tumor_stage_ind_pt3b = 1
										else:
											tumor_stage_list.append('pt3b')
											tumor_stage_ind = 1
											tumor_stage_ind_pt3b = 1
										break					
			if tumor_stage_ind_pt3a == 1: # fix pt3a case if you find "seminal" later in the report
				if 'seminal' in val:
					if idx - idx_extend > 35: # more than 35 words betweeen extend and seminal, which means seminal info is most likely irrelevant
						break
					context_to_look = text_oi_list[idx - 7 : idx + 7]
					if 'negative' not in context_to_look and 'no' not in context_to_look and 'not' not in context_to_look:
						tumor_stage_list.pop()
						tumor_stage_list.append('pt3b')
					else:
						tumor_stage_list.pop()
						tumor_stage_list.append('pt3a')
					break # break it here since you only want to take a look at seminal context right after extend
		# update tumor stage
		if not tumor_stage_ind:
			tumor_stage_list.append(None)
	df_pathology_with_rp['pT_stage_context'] = tumor_stage_list
	return df_pathology_with_rp

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
		df_pathology_with_rp = pd.read_csv('df_pathology_with_extracted_info.tsv', sep = '\t')
		df_pathology_with_rp.set_index(df_pathology_with_rp.columns[0], inplace = True)
	else:
		# hard-code rules (see : https://docs.google.com/document/d/1pRA2XcAxjbqbji8WHFfJD2bvwgo-EgryGLNRiq-gYT4/edit)
		# surgical margin
		hard_coded_margin_negative = ["margins are negative", "margins and seminal vesicles are uninvolved", "Margins uninvolved by invasive carcinoma", "does not extend to the inked peripheral resection margin", "extends to within 0.1 cm of the distal margin", "resection margins and seminal vesicles are free", "All inked margins and capsule are free of tumor", "margins, seminal vesicle and vas deferens are negative", "margins, seminal vesicles and vasa deferentia are negative for malignancy", "WITHIN 0.1 CM AT THE RIGHT ANTERIOR AND LEFT POSTERIOR RESECTION MARGINS",  "NOT PRESENT AT THE RIGHT POSTERIOR RESECTION MARGINS", "does not appear grossly to be present at the margin", "not grossly present at the margins", "tumor is not present at the evaluable soft tissue", "Margins appear negative", "margins of resection are free", "margins are free", "not present at the soft tissue resection margin", "margin is negative", "margins: uninvolved", "Margins: Negative", "resection margins, seminal vesicles and vas deferentia are free",  "within <1mm of the inked peripheral margin", "margins including urethral are negative", "within less than 0.1 cm of the right posterior resection margin", "No tumor present at the soft tissue resection margin.", "WITHIN LESS THAN 1MM FROM THE INKED MARGIN", "No tumor at the soft tissue resection margin.", "unequivocal carcinoma at the inked surface is not identified", "Surgical margins and vasa deferentia are negative.", "Negative resection margins", "within less than 1 mm of the inked surgical", "margins are free of carcinoma.", "within 0.1 cm of inked", "Margins: Uninvolved", "margins of resection and seminal vesicles are free of tumor", "<0.1 cm from the inked soft tissue resection margin"]
		hard_coded_margin_negative = [val.lower() for val in hard_coded_margin_negative]

		hard_coded_margin_positive = ["margins are positive", "is focally present at the peripheral inked margin", "carcinoma is present at the left anterior quadrant along a 1 mm", "and focally extends to the inked surface", "EXTENDING TO MULTIPLE INKED RESECTION MARGINS", "at least one, but not all margins are positive", "tumor is present at the margin", "carcinoma present at the margin", "extensively present at the inked", "extensively present at the soft tissue resection margin", "PRESENT AT THE RIGHT POSTERIOR RESECTION MARGIN", "focally present at the inked", "present at an inked resection margin", "extends through the capsule into periprostatic adipose tissue and to the inked surgical margin", "present at the distal inked margin", "PRESENT AT A RIGHT POSTERIOR INFERIOR INKED PROSTATIC SURGICAL MARGIN", "PRESENT AT THE PERIPHERIAL LEFT INFERIOR PROSTATE INKED MARGIN", "present at the inked RESECTION SURFACE.", "Except for the distal urethral margin (slide B), all other margins are free of tumor.", "INVOLVEMENT OF MARGIN"]
		hard_coded_margin_positive = [val.lower() for val in hard_coded_margin_positive]

		# tumor pT stage
		hard_coded_pt2_stage = ["Extraprostatic extension of tumor is not identified", "confined to the prostate", "Extraprostatic extension of tumor is not identified", "confined to the gland", "confined to prostate", "confined within the prostatic capsule", "involves the capsule but is confined to the prostate", "invades but does not transgress the capsule", "No capsular penetration demonstrated", "NOT INFILTRATING PERIPROSTATIC ADIPOSE TISSUE", "tumor does not transgress the prostatic capsule", "apparently localized", "invades into but not through the prostatic capsule", "extends to, but not through the prostatic capsule.", "the periprostatic and soft tissue resection margins are free of carcinoma", "no capsular invasion is present.", "extracapsular extension is not identified", "confined within the prostatic capsule", "TUMOR IS NOT SEEN OUTSIDE THE CAPSULE", "does not extend into periprostatic fat", "is confirmed to the prostate", "Extraprostatic extension is not present", "tumor does not transgress the prostatic capsule", "All inked margins and capsule are free of tumor", "no tumor involvement of the capsule", "Extraprostatic extension is not identified", "No extraprostatic extension of tumor"]
		hard_coded_pt2_stage = [val.lower() for val in hard_coded_pt2_stage]

		hard_coded_pt3a_stage = ["extends into the prostatic", "extends through the capsule into periprostatic", "extends slightly into periprostatic tissue", "LIMITED EXTRAPROSTATIC EXTENSION", "tumor extends into extraprostatic soft tissue", "tumor focally penetrates through the prostate capsule", "focally into the periprostatic fat"]
		hard_coded_pt3a_stage = [val.lower() for val in hard_coded_pt3a_stage]

		hard_coded_pt3b_stage = ["extends into extra prostatic soft tissues and the left seminal vesicle", "and involves the left seminal vesicle"]
		hard_coded_pt3b_stage = [val.lower() for val in hard_coded_pt3b_stage]

		# load the processed pathology report
		df_pathology = pd.read_csv('../data/Pat_total.csv', sep="|", low_memory=False)
		del df_pathology['Unnamed: 0']

		# inclusion/exclustion criteria based on words in each report
		idx_oi = [idx for idx, report in enumerate(df_pathology.Report_Text.values) if 'radical prostatectomy' in report.lower() or 'prostatect' in report.lower() or 'prostate radical resection' in report.lower() or 'prostate, radical resection' in report.lower() or ('prostat' in report.lower() and 'margin' in report.lower())]
		idx_to_exclude = [idx for idx, report in enumerate(df_pathology.Report_Text.values) if 'ureter' in report.lower() or 'cystect' in report.lower() or 'bone' in report.lower() or 'specimen type: prostate needle biopsy' in report.lower()] # or 'prostate, radical resection' in report.lower()]
		df_pathology_with_rp = df_pathology.iloc[np.sort(list(set(idx_oi) - set(idx_to_exclude)))].copy()

		# loading singlerp.csv which contains list of patients who received radical prostatectomy
		single_rp_patients = pd.read_csv('../data/singlerp.csv')
		df_pathology_with_rp = df_pathology_with_rp.loc[df_pathology_with_rp.EMPI.isin(single_rp_patients.EMPI.values)]

		# declare lists which store abstracted biopsy features
		primary_grade_list = []; secondary_grade_list = []; overall_grade_list = []; radical_pros = [];
		tumor_stage_list = []; margin_list = []; n1_list = []
		small_cell_carc = []; neuroendocrine_carc = []; adenocarcinoma = []
		stage_count = 0; margin_counter = 0; margin_counter_hard_coded = 0; pt_stage_counter_hard_coded = 0
		pt_stage_false_positive_terms = ['acceptor', 'transcripts']
		query_word = []

		# iterate through reports
		for text_oi, report_number in tqdm(zip(df_pathology_with_rp.Report_Text.values, df_pathology_with_rp.Report_Number.values), total = len(df_pathology_with_rp), desc = 'extracing info...'):
			radical_pros_ind = 0; gleason_reported = 0; primary_grade = None; secondary_grade = None; n1_val = None
			# all patients in the chosen cohort received rp
			radical_pros_ind = 1
			radical_pros.append(1)
			
			# make all words in lowercase for consistency
			text_oi = text_oi.lower()
			# create a list of all the words in each report
			text_oi_list = text_oi.split(' ')
			margin_ind = 0
			small_cell_carc_ind = 0; neuroendocrine_carc_ind = 0; adenocarcinoma_ind = 0; tumor_stage_ind = 0
			
			# Check for hard-coded rules (surgical margin)
			for val_margin_harcoded_pos in hard_coded_margin_positive:
				if val_margin_harcoded_pos in text_oi:
					margin_ind = 1
					margin_list.append(1)
					margin_counter_hard_coded += 1
					break

			for val_margin_harcoded_neg in hard_coded_margin_negative:
				if val_margin_harcoded_neg in text_oi:
					if margin_ind:
						margin_list.pop()
						margin_list.append(0)
					else:
						margin_ind = 1
						margin_list.append(0)
						margin_counter_hard_coded += 1
					break

			# Check for hard-coded rules (pT stage)
			for val_harcoded_pt2 in hard_coded_pt2_stage:
				if val_harcoded_pt2 in text_oi:
					tumor_stage_ind = 1
					tumor_stage_list.append('t2')
					pt_stage_counter_hard_coded += 1
					break

			for val_harcoded_pt3a in hard_coded_pt3a_stage:
				if val_harcoded_pt3a in text_oi:
					if tumor_stage_ind:
						print('(pT stage) presence of conflicting sentences in Report : ', report_number)
						print(val_harcoded_pt2)
						print(val_harcoded_pt3a)
						print('\n')
						tumor_stage_list.pop()
						tumor_stage_list.append('t3a')
					else:
						tumor_stage_ind = 1
						tumor_stage_list.append('t3a')
						pt_stage_counter_hard_coded += 1
					break

			for val_harcoded_pt3b in hard_coded_pt3b_stage:
				if val_harcoded_pt3b in text_oi:
					if tumor_stage_ind:
						print('(pT stage) presence of conflicting sentences in Report : ', report_number)
						print(val_harcoded_pt3a)
						print(val_harcoded_pt3b)
						print('\n')
						tumor_stage_list.pop()
						tumor_stage_list.append('t3b')
					else:
						tumor_stage_ind = 1
						tumor_stage_list.append('t3b')
						pt_stage_counter_hard_coded += 1
					break

			# iterate through words
			for idx, val in enumerate(text_oi_list):
				# get pT stage info
				if not tumor_stage_ind: # this ensures that the algorithm doesn't check the report again once tumor stage info has been extracted
					if 'stage' in val or 'pt' in val and 'pten' not in val and 'ptpn' not in val and val not in pt_stage_false_positive_terms:
						if 'clinical' not in text_oi_list[idx - 1]: # Do not include clinical stage for pT stage
							if 'n1' in ' '.join(text_oi_list[idx:idx+5]):
								n1_val = 1
							# if 'stage' and 'pt' are in the report: analyze at the context (i.e. next 10 words)
							for val_after in text_oi_list[idx:idx+10]:
								# remove punctuation
								val_filtered =  re.sub(r'[^\w\s]','',val_after)
								# remove any 'n' 
								val_filtered = val_filtered.replace('n', '')
								if 't' in val_filtered: # if t is in the word
									digit_detcted = 0
									for idx_char, char in enumerate(val_filtered): # check if number follows after t (e.g. t2b)
										if 't' == char:
											if idx_char < len(val_filtered) - 1:
												if val_filtered[idx_char+1].isdigit():
													digit_detcted = 1
													break #from the inner loop 
									if digit_detcted:
			                           	# remove punctuation and grab 2 chars after 't' to get stage
										tumor_stage_list.append(val_filtered[idx_char:idx_char+3])										
										tumor_stage_ind = 1
										break 
					# e.g.) get "t2a" when pT stage information doesn't appear after 'stage' or 'pT'
					elif len(val) == 3 and val[0] == 't' and val[1].isdigit() and not val[2].isdigit():
						if 'n1' in ' '.join(text_oi_list[idx:idx+5]):
							n1_val = 1
						tumor_stage_list.append(val)
						tumor_stage_ind = 1
				# get surfical margin info 
				if not margin_ind:
					if 'margin' in val:
						absence_ind_word = ['negative', 'free']
						presence_ind_word = ['positive']
						context_to_look = text_oi_list[idx : idx + 5]
						for word in presence_ind_word:
							for context_word in context_to_look:
								if word in context_word:
									margin_ind = 1
									margin_counter += 1
									margin_list.append(1)
									break
							if margin_ind:
								break
						# if not margin_ind, check if it contains negative surgical margin keywords
						margin_ind_aux = 0
						for word in absence_ind_word:
							for context_word in context_to_look:
								if word in context_word:
									if margin_ind:
										print('(Margin) presence of both positive and negative words in Report : ', report_number)
										print('\n')
										margin_list.pop()
										margin_list.append(0)
										# breakpoint()
									else:
										margin_ind = 1
										margin_counter += 1
										margin_list.append(0)
									margin_ind_aux = 1
									break
							if margin_ind_aux:
								break
	
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
				
				# get gleason score and overall grade
				if 'gleason' in val and radical_pros_ind and not gleason_reported:
					# declare indicators 
					found_primary_grade = 0; found_secondary_grade = 0; found_plus_sign = 0; found_one_dash = 0
					# get context around word "gleason"
					# heuristically choose 15. Gleason score we're after should be included in the next 15 words
					raw_gleason_text = ''.join(text_oi_list[idx:idx + 15])
					# remove text in () or [] : e.g. Score 3 (60%) + 4 (40%) = 7/10 -> Score 3 + 4 = 7/10
					raw_gleason_text = re.sub("[\(\[].*?[\)\]]", "", raw_gleason_text)
					gleason_reported = 1 # this ensures we only look at the first appearance of gleason score

					# to capture gleason 3/5 where primary 3 and secondary 3, check if there's only one slash without any + or - signs. 
					slash_counter = 0; dash_counter = 0; plus_sign_counter = 0; only_one_slash = 0
					raw_gleason_text_slash = ''.join(text_oi_list[idx:idx + 15]) # get context
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
									# in this case primary grade = secondary grade = 3
									primary_grade = int(raw_gleason_grade_dash_split[0][0])
									secondary_grade = int(raw_gleason_grade_dash_split[0][2])
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
				if primary_grade + secondary_grade <= 6:
					overall_grade_group = 1
				elif primary_grade + secondary_grade == 7:
					if primary_grade == 3:
						overall_grade_group = 2
					elif primary_grade == 4:
						overall_grade_group = 3
				elif primary_grade + secondary_grade == 8:
					overall_grade_group = 4
				elif primary_grade + secondary_grade > 8:
					overall_grade_group = 5
				else:
					overall_grade_group = None
			else:
				overall_grade_group = None
			
			# update the grade group lists
			primary_grade_list.append(primary_grade)
			secondary_grade_list.append(secondary_grade)
			overall_grade_list.append(overall_grade_group)
			
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
				
			# store tumor stage
			if not tumor_stage_ind:
				tumor_stage_list.append(None)
			# store surgical margin
			if not margin_ind:
				margin_list.append(None)

			# store n1 pT stage
			if n1_val == 1:
				n1_list.append(1)
			else:
				n1_list.append(0)

		# breakpoint()
		df_pathology_with_rp.loc[df_pathology_with_rp.index, 'n1_indicator'] = n1_list
		# df_pathology_with_rp.loc[df_pathology_with_rp.index, 'num_pos_cores'] = num_pos_cores_list
		df_pathology_with_rp.loc[df_pathology_with_rp.index, 'primary_grade'] = primary_grade_list
		df_pathology_with_rp.loc[df_pathology_with_rp.index, 'secondary_grade'] = secondary_grade_list
		df_pathology_with_rp.loc[df_pathology_with_rp.index, 'overall_grade_group'] = overall_grade_list
		# manual gleason score fix (primary, secondary) : # (3,7) => (3,4) Report MS04R67776, (3,6) => (3,3) Report S0009791E
		df_pathology_with_rp.loc[df_pathology_with_rp.Report_Number == 'MS04R67776', 'secondary_grade'] = 4
		df_pathology_with_rp.loc[df_pathology_with_rp.Report_Number == 'MS04R67776', 'overall_grade_group'] = 2

		df_pathology_with_rp.loc[df_pathology_with_rp.Report_Number == 'S0009791E', 'secondary_grade'] = 3
		df_pathology_with_rp.loc[df_pathology_with_rp.Report_Number == 'S0009791E', 'overall_grade_group'] = 1

		df_pathology_with_rp.loc[df_pathology_with_rp.index, 'radical_prostatectomy'] = radical_pros

		df_pathology_with_rp.loc[df_pathology_with_rp.index, 'small_cell_carc'] = small_cell_carc
		df_pathology_with_rp.loc[df_pathology_with_rp.index, 'neuroendocrine_carc'] = neuroendocrine_carc
		df_pathology_with_rp.loc[df_pathology_with_rp.index, 'adenocarcinoma'] = adenocarcinoma

		df_pathology_with_rp.loc[df_pathology_with_rp.index, 'margin'] = margin_list
		# manually fix margin based on Madhur's audit (see : https://docs.google.com/document/d/1aUFjvz8bumhCnSUDJ8w4CtKoNmGBKYTGeq0cMqoOOl0/edit)
		positive_margin_reports = ['S94W17361', 'S95R18902', 'MS09K38970', 'BS09R32983', 'S16-17016']
		negative_margin_reports = ['S9830284L', 'S0133510G', 'S0020803D', 'BS07T08214', 'MS09G67481', 'MS08K66009', 'S9818403C']
		df_pathology_with_rp.loc[df_pathology_with_rp.Report_Number.isin(positive_margin_reports), 'margin'] = 1
		df_pathology_with_rp.loc[df_pathology_with_rp.Report_Number.isin(negative_margin_reports), 'margin'] = 0

		df_pathology_with_rp.loc[df_pathology_with_rp.index, 'pT_stage'] = [None if stage is None else 'p' + stage for stage in tumor_stage_list]
		pt_stage_list = []
		remove_chars = ['o', 'x', 's']
		# additional curation of pT stage
		for pt_stage in df_pathology_with_rp.pT_stage.values:
			if pt_stage is None:
				pt_stage_list.append(None)
				continue
			if pt_stage == 'pt30':
				pt_stage_list.append('pt3')
				continue
			elif pt_stage == 'pt20':
				pt_stage_list.append('pt2')
				continue
			elif pt_stage == 'pt00':
				pt_stage_list.append('pt0')
				continue
			elif pt_stage == 'pt40':
				pt_stage_list.append('pt4')
				continue
			elif pt_stage == 'pt31': # manual override in accesion number : bs13j21799
				pt_stage_list.append('pt3a')
				continue
			elif pt_stage == 'pt25': 
				pt_stage_list.append('pt2')
				continue
			elif pt_stage == 'pt26': 
				pt_stage_list.append('pt2')
				continue
			pt_stage_filtered = pt_stage
			for char in remove_chars:
				pt_stage_filtered = pt_stage_filtered.replace(char, '')
			pt_stage_list.append(pt_stage_filtered)
		df_pathology_with_rp.loc[df_pathology_with_rp.index, 'pT_stage'] = pt_stage_list
		df_pathology_with_rp = get_tumor_stage_info_via_context(df_pathology_with_rp)
		# test and statistics
		df_pathology_with_rp_comparison = df_pathology_with_rp.dropna(subset = ['pT_stage_context'])
		df_pathology_with_rp_comparison = df_pathology_with_rp_comparison.dropna(subset = ['pT_stage'])
		df_compare = df_pathology_with_rp_comparison[['Report_Number', 'pT_stage', 'pT_stage_context']]
		print('\n')
		print('Discrepancies between the numeric and context-based methods among the reports which have the both information available')
		count = 0
		for entry in df_compare.values:
			compare_len = len(entry[2])
			if entry[2] != entry[1][0:compare_len]:
				print(entry)
				count = count + 1
		print('AJCC numeric-based pT stage vs. context-based pT stage')
		print('accuracy : ', (len(df_compare) - count)/len(df_compare))
		print('\n')
		print('Number of reports in each pT stage : ')
		# combining numeric based and context based pT stage 
		pt_stage_combined = []; context_based_counter = 0; numeric_based_counter = 0;
		for pt_stage_numeric, pT_stage_context in zip(df_pathology_with_rp.pT_stage.values, df_pathology_with_rp.pT_stage_context.values):
			if pt_stage_numeric is not None:
				pt_stage_combined.append(pt_stage_numeric)
				numeric_based_counter += 1
			elif pT_stage_context is not None:
				pt_stage_combined.append(pT_stage_context)
				context_based_counter += 1
			else:
				pt_stage_combined.append(None)
		df_pathology_with_rp.loc[df_pathology_with_rp.index, 'pT_stage_combined'] = pt_stage_combined

		# choose only non-N1 
		df_pathology_with_rp = df_pathology_with_rp.loc[df_pathology_with_rp.n1_indicator == 0]

		# filter out outlier ptstage
		outlier_pt_stages = ['pt32', 'pt3z', 'pt21', 'pt2.', 'pt2t', 'pt37', 'pt53', 'pt3d', 'pt1.', 'pt2d']
		pt_stages_oi = set(df_pathology_with_rp.pT_stage_combined.values) - set(outlier_pt_stages)
		df_pathology_with_rp = df_pathology_with_rp.loc[df_pathology_with_rp.pT_stage_combined.isin(pt_stages_oi)]

		# get pT stage stats
		df_count_stats = pd.value_counts(df_pathology_with_rp.pT_stage_combined.values, sort = True)
		print(df_count_stats)
		print('# of numeric-based pT stages : ', numeric_based_counter)
		print('# of context-based pT stages : ', context_based_counter)
		print('Among both methods, {0} is the number of reports with hard-coded pT stage sentences'.format(pt_stage_counter_hard_coded))
		print('total = ', sum(df_count_stats.values))
		print('\n')

		# surgical margin stats
		df_count_stats_margin = pd.value_counts(df_pathology_with_rp.margin.values, sort = True)
		print('Number of reports with margin indicator 1 (postiive) or 0 (negative) : ')
		print(df_count_stats_margin)
		print('Among those, {0} is the number of reports with hard-coded margin indicator sentences'.format(margin_counter_hard_coded))
		print('total = ', sum(df_count_stats_margin.values))
		print('\n')
	
	# overall grade stats
	df_count_stats_primary_grade = pd.value_counts(df_pathology_with_rp.primary_grade.values, sort = True)
	df_count_stats_secondary_grade = pd.value_counts(df_pathology_with_rp.secondary_grade.values, sort = True)
	df_count_stats_overall_grade = pd.value_counts(df_pathology_with_rp.overall_grade_group.values, sort = True)
	print('Number of reports per each primary grade : ')
	print(df_count_stats_primary_grade)
	print('Number of reports per each secondary grade : ')
	print(df_count_stats_secondary_grade)
	print('Number of reports per each overall grade group : ')
	print(df_count_stats_overall_grade)
	print('total = ', df_count_stats_overall_grade.sum())

	print('Number of unique reports with radical prostatectomy : ', len(df_pathology_with_rp))

	print('Number of unique reports with radical prostatectomy per center : ')
	print(pd.value_counts(df_pathology_with_rp.MRN_Type.values, sort = True))

	print('Number of unique patients with radical prostatectomy : ', len(set(df_pathology_with_rp.EMPI.values)))

	df_pathology_with_rp_oi_overall_grade = df_pathology_with_rp.dropna(subset = ['overall_grade_group'])
	print('Number of unique patients with radical prostatectomy and overall grade : ', len(set(df_pathology_with_rp_oi_overall_grade.EMPI.values)))

	df_pathology_with_rp_oi_with_pt_stage = df_pathology_with_rp_oi_overall_grade.dropna(subset = ['pT_stage_combined'])
	df_pathology_with_rp_oi_with_pt_stage_wo_grade = df_pathology_with_rp.dropna(subset = ['pT_stage_combined'])
	print('Number of unique patients with radical prostatectomy, overall grade, and pT_stage : ', len(set(df_pathology_with_rp_oi_with_pt_stage.EMPI.values)))
	print('Number of unique patients with radical prostatectomy, and pT_stage : ', len(set(df_pathology_with_rp_oi_with_pt_stage_wo_grade.EMPI.values)))


	df_pathology_with_rp_oi_with_pt_stage_margin = df_pathology_with_rp_oi_with_pt_stage.dropna(subset = ['margin'])
	print('Number of unique patients with radical prostatectomy, overall grade, pT_stage, and margin : ', len(set(df_pathology_with_rp_oi_with_pt_stage_margin.EMPI.values)))
	
	print('\n')
	print('\n')
	breakpoint()
	continue_ = True; sample_rp = True
	while continue_:
		print('\n')
		print('Extracted info : ')
		if sample_rp:
			sampled_entry = df_pathology_with_rp_oi_with_pt_stage_margin.loc[df_pathology_with_rp_oi_with_pt_stage_margin.radical_prostatectomy == 1].sample(n = 1)
		else:
			sampled_entry = df_pathology_with_rp_oi_with_pt_stage_margin.sample(n = 1)
		sampled_entry_extracted_oi = sampled_entry[['Report_Number', 'radical_prostatectomy', 'primary_grade', 'secondary_grade', 'overall_grade_group', 'pT_stage_combined', 'margin', 'small_cell_carc', 'neuroendocrine_carc', 'adenocarcinoma']]
		print(sampled_entry_extracted_oi.iloc[0])
		print('\n')
		print('Pathology report : ')
		print(sampled_entry.Report_Text.values[0])
		print('\n')
		continue_ = input('Would you like to continue? (Y/N) : ') == 'Y'
	
	export = input('Export the pathology report datafame with the extracted info in the current directory? (Y/N) : ')
	if export not in {'Y', 'N'}:
		raise KeyError('Input Y or N only')
	if export == 'Y':
		df_pathology_with_rp.to_csv('df_pathology_with_extracted_info.tsv', sep = '\t', index = True)
		df_pathology_with_rp_oi_with_pt_stage_margin.to_csv('df_pathology_with_grade_pt_stage_margin.tsv', sep = '\t', index = True)
		print('\n')
		print('The processed file exported to ' + os.getcwd() + '/df_pathology_with_extracted_info.tsv')
		print('\n')
	return

if __name__ == "__main__":
	explain = """
	Written by Intae Moon
	Date : April 26, 2021

	Description :

	The following script extracts 
	1. primary grade
	2. secondary grade
	3. overall grade group
	4. histology
	5. pT stage
	6. n1 status
	7. surgical margin

	from RP pathology notes

	See : https://docs.google.com/document/d/1pRA2XcAxjbqbji8WHFfJD2bvwgo-EgryGLNRiq-gYT4/edit?usp=sharing for more details
	"""
	print("\n")
	print(explain)
	print("\n")
	input('Type any button if you would like to continue')
	main()
