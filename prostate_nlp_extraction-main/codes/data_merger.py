import numpy as np
import pandas as pd
from tqdm import tqdm
import pickle
import re

def func_get_outcome(df_oi, empi_to_death_date_dic, empi_to_date_oi_dic, get_biopsy_based_feature = False, empi_to_last_psa_check_date_dic_biopsy_based = None, clean_mode = False, extra_clean = False):
	if get_biopsy_based_feature:
		with open('../data/empi_to_rp_date_dic.pkl', 'rb') as handle:
			empi_to_rp_date_dic = pickle.load(handle)
		df_biopsy_outcome = df_oi.copy()
		# as a heuristic
		if clean_mode:
			df_biopsy_outcome['Report_Date_Time'] = pd.to_datetime(df_biopsy_outcome['Report_Date_Time'].values)
			last_ndi_update = pd.to_datetime('2017-11-1') # latest lab result records
			# extra clean 10-years before the last death update date
			if extra_clean:
				last_ndi_update = pd.to_datetime('2007-11-1')
			df_biopsy_outcome = df_biopsy_outcome.loc[df_biopsy_outcome.Report_Date_Time < last_ndi_update]
			df_biopsy_outcome['time_to_death_in_month'] = None; df_biopsy_outcome['death_date'] = None; df_biopsy_outcome['death_ind'] = None
			for empi in df_biopsy_outcome.index:
				if empi in empi_to_death_date_dic.keys():
					death_date = empi_to_death_date_dic[empi]
					if death_date <= last_ndi_update:
						df_biopsy_outcome.at[empi, 'death_date'] = death_date
						df_biopsy_outcome.at[empi, 'death_ind'] = 1
						df_biopsy_outcome.at[empi, 'time_to_death_in_month'] = (death_date - empi_to_date_oi_dic[empi]).days/30.44
					else:
						df_biopsy_outcome.at[empi, 'death_ind'] = 0
						df_biopsy_outcome.at[empi, 'time_to_death_in_month'] = (last_ndi_update - empi_to_date_oi_dic[empi]).days/30.44
				else:
					df_biopsy_outcome.at[empi, 'death_ind'] = 0
					df_biopsy_outcome.at[empi, 'time_to_death_in_month'] = (last_ndi_update - empi_to_date_oi_dic[empi]).days/30.44
		else:
			df_biopsy_outcome['death_date'] = [empi_to_death_date_dic[empi] if empi in empi_to_death_date_dic.keys() else None for empi in df_biopsy_outcome.index]
			df_biopsy_outcome['death_ind'] = [1 if empi in empi_to_death_date_dic.keys() else 0 for empi in df_biopsy_outcome.index]
			last_ndi_update = pd.to_datetime('2018-5-1') # latest lab result records
			# breakpoint()
			time_to_death_in_month_list = []
			for empi in df_biopsy_outcome.index:
				if empi in empi_to_death_date_dic.keys():
					time_to_death_in_month_list.append((empi_to_death_date_dic[empi] - empi_to_date_oi_dic[empi]).days/30.44)
				elif empi in empi_to_last_psa_check_date_dic_biopsy_based.keys(): 
					time_to_death_in_month_list.append((empi_to_last_psa_check_date_dic_biopsy_based[empi] - empi_to_date_oi_dic[empi]).days/30.44)
				else:
					time_to_death_in_month_list.append((last_ndi_update - empi_to_date_oi_dic[empi]).days/30.44)
			df_biopsy_outcome['time_to_death_in_month'] = time_to_death_in_month_list
	else:
		df_lab_valid = df_oi.loc[df_oi.first_post_op_psa_valid == 1].copy()
		# df_lab_valid['death_date'] = None
		df_lab_valid['death_date'] = [empi_to_death_date_dic[empi] if empi in empi_to_death_date_dic.keys() else None for empi in df_lab_valid.index]
		df_lab_valid['death_ind'] = [1 if empi in empi_to_death_date_dic.keys() else 0 for empi in df_lab_valid.index]
		# as a heuristic
		max_death_date = pd.to_datetime('2020-10-06 16:10:00') # latest lab result records
		# breakpoint()
		df_lab_valid['time_to_death_in_month'] = [(empi_to_death_date_dic[empi] - empi_to_date_oi_dic[empi]).days/30.44 if empi in empi_to_death_date_dic.keys() else (max_death_date - empi_to_date_oi_dic[empi]).days/30.44 for empi in df_lab_valid.index]
		# if death date < bcr_date, censor patients at death date
		if not get_biopsy_based_feature:
			df_lab_valid['time_to_bcr_in_month'] = df_lab_valid['bcr_date_minus_rp_date_in_days'].values/30.44
			update_counter = 0
			for empi, death_date, bcr_date in zip(df_lab_valid.index, df_lab_valid.death_date.values, df_lab_valid.bcr_date.values):
				if not pd.isnull(death_date) and not pd.isnull(bcr_date):
					bcr_date_converted = pd.to_datetime(bcr_date)
					if death_date < bcr_date_converted:
						df_lab_valid.at[empi, 'time_to_bcr_in_month'] = df_lab_valid.at[empi, 'time_to_bcr_in_month'] - (bcr_date_converted - death_date).days/30.44
						update_counter += 1
	# pass
	if get_biopsy_based_feature:
		df_outcome_oi = pd.DataFrame(0.0, index = df_biopsy_outcome.index, columns = ['death_ind', 'death_ind_5_year', 'death_ind_10_year', 'time_to_death_in_month', 'biopsy_date'])
		for empi, death_ind, time_to_death, biopsy_date in zip(df_outcome_oi.index, df_biopsy_outcome.death_ind.values, df_biopsy_outcome.time_to_death_in_month.values, df_biopsy_outcome.Report_Date_Time.values):
			df_outcome_oi.at[empi, 'time_to_death_in_month'] = time_to_death
			df_outcome_oi.at[empi, 'death_ind'] = death_ind
			df_outcome_oi.at[empi, 'biopsy_date'] = pd.to_datetime(biopsy_date)
			if empi in empi_to_rp_date_dic.keys():
				rp_date = empi_to_rp_date_dic[empi]
				df_outcome_oi.at[empi, 'rp_date'] = rp_date
				df_outcome_oi.at[empi, 'rp_date_minus_biopsy_date_in_days'] = (rp_date - pd.to_datetime(biopsy_date)).days
			if death_ind == 1 and time_to_death <= (365.25*5)/30.44:
				df_outcome_oi.at[empi, 'death_ind_5_year'] = 1
				df_outcome_oi.at[empi, 'death_ind_10_year'] = 1
			elif death_ind == 1 and time_to_death <= (365.25*10)/30.44:
				df_outcome_oi.at[empi, 'death_ind_10_year'] = 1
	else:
		df_outcome_oi = pd.DataFrame(0.0, index = df_lab_valid.index, columns = ['bcr_ind', 'bcr_ind_5_year', 'bcr_ind_10_year', 'time_to_bcr_in_month', 'death_ind', 'death_ind_5_year', 'death_ind_10_year', 'time_to_death_in_month'])
		for empi, bcr_ind, time_to_bcr, death_ind, time_to_death in zip(df_outcome_oi.index, df_lab_valid.bcr_post_op.values, df_lab_valid.time_to_bcr_in_month.values, df_lab_valid.death_ind.values, df_lab_valid.time_to_death_in_month.values):
			df_outcome_oi.at[empi, 'time_to_bcr_in_month'] = time_to_bcr
			df_outcome_oi.at[empi, 'time_to_death_in_month'] = time_to_death
			df_outcome_oi.at[empi, 'bcr_ind'] = bcr_ind
			df_outcome_oi.at[empi, 'death_ind'] = death_ind

			if bcr_ind == 1 and time_to_bcr <= (365.25*5)/30.44:
				df_outcome_oi.at[empi, 'bcr_ind_5_year'] = 1
				df_outcome_oi.at[empi, 'bcr_ind_10_year'] = 1
			elif bcr_ind == 1 and time_to_bcr <= (365.25*10)/30.44:
				df_outcome_oi.at[empi, 'bcr_ind_10_year'] = 1

			if death_ind == 1 and time_to_death <= (365.25*5)/30.44:
				df_outcome_oi.at[empi, 'death_ind_5_year'] = 1
				df_outcome_oi.at[empi, 'death_ind_10_year'] = 1
			elif death_ind == 1 and time_to_death <= (365.25*10)/30.44:
				df_outcome_oi.at[empi, 'death_ind_10_year'] = 1

	# breakpoint()
	return df_outcome_oi

def get_psa_prior_to_rp(df_lab_total, empi_to_date_oi_dic, empis_oi, export_dic = False, mode = 'rp'):#, empis_oi, get_biopsy_based_feature = False):
	# breakpoint()
	if mode in ['rp', 'radiation']:
		pass
	else:
		raise KeyError('Only supports rp or radiation atm...')

	df_lab_results_oi = df_lab_total.loc[df_lab_total.EMPI.isin(empis_oi)].copy()
	df_lab_results_oi['Seq_Date_Time'] = pd.to_datetime(df_lab_results_oi.Seq_Date_Time)

	unique_empi_vals = np.unique(df_lab_results_oi.EMPI.values)
	units_oi = ['ng/mL', 'ng/ml', 'NG/ML', 'MCG/L']
	tests_oi = ['PSA (Total and Screening)', 'PSA Monitoring', 'PSA (ultrasensitive)']
	values_to_ignore = ['<0.7', '<0.5']
	values_to_consider = ['FREE PSA = 0.47', 'ASSAY RANGE, 0.1', '[0.1]', '[1.16]', '[0.47]', 'AS.01', ',ASSAY RANGE, 0.1', '[.01]']
	offset_days = 3650
	lab_val_thresh = 0.2
	empi_recurrence_psa_date_dic = {}; lab_val_abnomral_list = []
	for empi in tqdm(unique_empi_vals, desc = 'getting prior PSA...'):
		# if empi == 100003779:
		# 	breakpoint()
		df_lab_results_oi_per_empi = df_lab_results_oi.loc[df_lab_results_oi.EMPI == empi].copy()
		rp_date = empi_to_date_oi_dic[empi]
		df_lab_results_oi_per_empi_filtered = df_lab_results_oi_per_empi.loc[df_lab_results_oi_per_empi.Seq_Date_Time >= rp_date - pd.DateOffset(days = offset_days)].copy()
		df_lab_results_oi_per_empi_filtered.sort_values(by = 'Seq_Date_Time', ascending = True, inplace = True)
		lab_val_list = []
		for lab_val, test_date, unit, test_type in zip(df_lab_results_oi_per_empi_filtered.Result.values, df_lab_results_oi_per_empi_filtered.Seq_Date_Time.values, df_lab_results_oi_per_empi_filtered.Reference_Units.values, df_lab_results_oi_per_empi_filtered.Group_Id.values):
			if not pd.isnull(lab_val) and not pd.isnull(unit) and not pd.isnull(test_type):
				if unit in units_oi and test_type in tests_oi: # ng/ml only
					if '<' in lab_val and lab_val not in values_to_ignore:
						# if lab_val == '<0.5':
						# 	lab_val = '0.5777'
						# if empi == 100003779:
						# 	breakpoint()
						try: 
							lab_val_filtered = re.findall('\d*\.?\d+', lab_val)[0] # only recover decimal point 
							lab_val_float = float(lab_val_filtered)
						except:
							lab_val_abnomral_list.append(lab_val)
							continue
						# try: # convert
						if lab_val_float > lab_val_thresh:
#                             days_oi = (test_date - rp_date).days
							lab_val_list.append((lab_val_float, pd.to_datetime(test_date), (pd.to_datetime(test_date) - rp_date).days, 1, 'abnormal', 'after rp' if test_date >= rp_date else 'before rp'))
						else:
							lab_val_list.append((lab_val_float, pd.to_datetime(test_date), (pd.to_datetime(test_date) - rp_date).days, 1, 'normal', 'after rp' if test_date >= rp_date else 'before rp'))
						
					else:
						# handles some corner cases
						if 'ng/mL' in lab_val or lab_val in values_to_consider:
							# print(lab_val)
							# breakpoint()
							lab_val = re.findall('\d*\.?\d+', lab_val)[0]
						try:
							lab_val_float = float(lab_val)
						except:
							lab_val_abnomral_list.append(lab_val)
							continue
						if lab_val_float > lab_val_thresh:
#                             breakpoint()
							lab_val_list.append((lab_val_float, pd.to_datetime(test_date), (pd.to_datetime(test_date) - rp_date).days, 0, 'abnormal', 'after rp' if test_date >= rp_date else 'before rp'))
						else:
							lab_val_list.append((lab_val_float, pd.to_datetime(test_date), (pd.to_datetime(test_date) - rp_date).days, 0, 'normal', 'after rp' if test_date >= rp_date else 'before rp'))
						# except:
						# 	lab_val_abnomral_list.append(lab_val)
						# 	continue
		if len(lab_val_list) > 0:
			empi_recurrence_psa_date_dic[empi] = lab_val_list

	# print('lab_val_abnomral_list: ', pd.value_counts(lab_val_abnomral_list,sort = True))
	# breakpoint()
	# # obtain the dataframe 
	# df_outcome = pd.DataFrame([], index = empi_recurrence_psa_date_dic.keys(), columns = ['first normal psa level after RP', 'first normal test date after RP', 'flag for < (normal)',
	# 																					  'first normal test date - RP date (months)', 'first abnormal psa level after RP', 
	# 																					  'first abnormal test date after RP', 'flag for < (abnormal)', 'first abnormal test date - RP date (months)', 
	# 																					  'RP date', 'PSA levels before RP (10 year)', 'PSA levels before RP dates', 'PSA levels after RP', 'PSA levels after RP dates', 'Last PSA test date'])
	if mode == 'rp':
		df_outcome = pd.DataFrame([], index = empi_recurrence_psa_date_dic.keys(), columns = ['rp_date', 'first_post_op_psa', 'first_post_op_psa_date', 'first_post_op_psa_date_minus_rp_date_in_days', 
																								'first_post_op_psa_valid', 'post_op_psa_trajectory', 'post_op_psa_trajectory_dates', 'hormone_post_op', 'bcr_post_op',
																								'bcr_date', 'bcr_date_minus_rp_date_in_days', 'pre_op_psa_trajectory_10_year_window', 'pre_op_psa_trajectory_dates', 'pre_op_psa_trajectory_dates_relative']) # bcr_rp is set to NaN if they first_post_op_psa_date_minus_rp_date_in_days > 6 month																		  
	else:
		df_outcome = pd.DataFrame([], index = empi_recurrence_psa_date_dic.keys(), columns = ['radiation_date', 'first_post_op_psa', 'first_post_op_psa_date', 'first_post_op_psa_date_minus_radiation_date_in_days', 
																								'first_post_op_psa_valid', 'post_op_psa_trajectory', 'post_op_psa_trajectory_dates', 'nadir_psa', 'nadir_psa_date', 'nadir_psa_date_minus_radiation_date_in_days', 'min_post_op_psa', 'min_post_op_psa_date', 'min_post_op_psa_date_minus_radiation_date_in_days', 'nadir_psa_v2', 'nadir_psa_date_v2', 'nadir_psa_date_minus_radiation_date_in_days_v2' , 'bcr_post_rad', # failure cases will have NaN nadir PSA
																								'bcr_date', 'bcr_date_minus_radiation_date_in_days', 'censor_date', 'censor_date_minus_radiation_date_in_days', 'pre_op_psa_trajectory_10_year_window', 'pre_op_psa_trajectory_dates', 'pre_op_psa_trajectory_dates_relative']) # bcr_rp is set to NaN if they first_post_op_psa_date_minus_rp_date_in_days > 6 month																		  
	empi_to_last_psa_check_date_dic = {}
	# breakpoint()
	# if get_biopsy_based_feature:
	# 	df_outcome['rp_failure'] = None
	lab_val_thresh_validity = 0.2
	for empi, list_oi in empi_recurrence_psa_date_dic.items():
		first_post_op_recorded = 0; record_bcr_event = 0
		RP_date = empi_to_date_oi_dic[empi]    
		psa_levels_before_rp = []; psa_levels_before_rp_dates = []; psa_levels_before_rp_dates_relative = []
		psa_levels_after_rp = []; psa_levels_after_rp_dates = []
		for tuple_oi in list_oi:
			if tuple_oi[5] == 'before rp':
				psa_levels_before_rp.append(tuple_oi[0])
				psa_levels_before_rp_dates.append(tuple_oi[1])
				# breakpoint()
				psa_levels_before_rp_dates_relative.append(-1 * np.round((RP_date - tuple_oi[1]).days/30.44))
			else: 
				if not first_post_op_recorded:
					# if (tuple_oi[1] - RP_date).days <= 6*30.44:
					df_outcome.at[empi, 'first_post_op_psa'] = tuple_oi[0]
					df_outcome.at[empi, 'first_post_op_psa_date'] = tuple_oi[1]
					if mode == 'rp':
						df_outcome.at[empi, 'first_post_op_psa_date_minus_rp_date_in_days'] = (tuple_oi[1] - RP_date).days
					else:
						df_outcome.at[empi, 'first_post_op_psa_date_minus_radiation_date_in_days'] = (tuple_oi[1] - RP_date).days
					if (tuple_oi[1] - RP_date).days <= 6*30.44:
						if tuple_oi[0] <= lab_val_thresh_validity:
							df_outcome.at[empi, 'first_post_op_psa_valid'] = 1 # rp initially worked for them
							# df_outcome.at[empi, 'rp_failure'] = 0
						else:
							# if get_biopsy_based_feature:
							# 	df_outcome.at[empi, 'first_post_op_psa_valid'] = 1 # it is still valid. but need to be flagged for rp failure 
							# 	df_outcome.at[empi, 'rp_failure'] = 1
							# else:
							df_outcome.at[empi, 'first_post_op_psa_valid'] = 0 # rp didnt' work for them, most likely patients were put on other therapy
					else:
						df_outcome.at[empi, 'first_post_op_psa_valid'] = 0 # don't know if RP worked for them or not
					psa_levels_after_rp.append(tuple_oi[0])
					psa_levels_after_rp_dates.append(tuple_oi[1])
					first_post_op_recorded = 1
				else:
					psa_levels_after_rp.append(tuple_oi[0])
					psa_levels_after_rp_dates.append(tuple_oi[1])

				# if tuple_oi[4] == 'normal' and not first_normal_psa_recorded:
				# 	df_outcome.at[empi, 'first normal psa level after RP'] = tuple_oi[0]
				# 	df_outcome.at[empi, 'first normal test date after RP'] = tuple_oi[1]
				# 	df_outcome.at[empi, 'first normal test date - RP date (months)'] = (tuple_oi[1] - RP_date).days/29.53
				# 	df_outcome.at[empi, 'flag for < (normal)'] = tuple_oi[3]
				# 	first_normal_psa_recorded = 1
				# elif tuple_oi[4] == 'abnormal' and not first_abnormal_psa_recorded: #abnormal
				# 	df_outcome.at[empi, 'first abnormal psa level after RP'] = tuple_oi[0]
				# 	df_outcome.at[empi, 'first abnormal test date after RP'] = tuple_oi[1]
				# 	df_outcome.at[empi, 'first abnormal test date - RP date (months)'] = (tuple_oi[1] - RP_date).days/29.53
				# 	df_outcome.at[empi, 'flag for < (abnormal)'] = tuple_oi[3]
				# 	first_abnormal_psa_recorded = 1
		df_outcome.at[empi, 'pre_op_psa_trajectory_10_year_window'] = psa_levels_before_rp
		df_outcome.at[empi, 'pre_op_psa_trajectory_dates'] = psa_levels_before_rp_dates
		df_outcome.at[empi, 'post_op_psa_trajectory'] = psa_levels_after_rp
		df_outcome.at[empi, 'post_op_psa_trajectory_dates'] = psa_levels_after_rp_dates
		df_outcome.at[empi, 'pre_op_psa_trajectory_dates_relative'] = psa_levels_before_rp_dates_relative
		if mode == 'rp':
			df_outcome.at[empi, 'rp_date'] = RP_date
		else:
			df_outcome.at[empi, 'radiation_date'] = RP_date
			if len(psa_levels_after_rp) > 0:
				min_post_rad_psa = np.min(psa_levels_after_rp)
				min_post_rad_psa_date = psa_levels_after_rp_dates[np.argmin(psa_levels_after_rp)]
				found_min = 0; found_nadir = 0
				for post_psa, post_psa_date in zip(psa_levels_after_rp, psa_levels_after_rp_dates):
					if post_psa <= lab_val_thresh_validity and not found_nadir:
						df_outcome.at[empi, 'nadir_psa'] = post_psa
						df_outcome.at[empi, 'nadir_psa_date'] = post_psa_date
						df_outcome.at[empi, 'nadir_psa_date_minus_radiation_date_in_days'] = (post_psa_date - RP_date).days
						found_nadir = 1
						if found_min: # only if you already have found min, break
							# if df_outcome.at[empi, 'min_post_op_psa'] <= lab_val_thresh_validity:		
							# if empi == 100000339:
							# 	breakpoint()						
							break
					if post_psa == min_post_rad_psa:
						found_min = 1
						df_outcome.at[empi, 'min_post_op_psa'] = min_post_rad_psa
						df_outcome.at[empi, 'min_post_op_psa_date'] = min_post_rad_psa_date
						df_outcome.at[empi, 'min_post_op_psa_date_minus_radiation_date_in_days'] = (min_post_rad_psa_date - RP_date).days
						if found_nadir:
							df_outcome.at[empi, 'nadir_psa_v2'] = min_post_rad_psa
							df_outcome.at[empi, 'nadir_psa_date_v2'] = min_post_rad_psa_date
							df_outcome.at[empi, 'nadir_psa_date_minus_radiation_date_in_days_v2'] = (min_post_rad_psa_date - RP_date).days
							# if empi == 100000339:
							# 	breakpoint()	
							break
		if len(psa_levels_after_rp_dates) > 0:
			empi_to_last_psa_check_date_dic[empi] = psa_levels_after_rp_dates[-1]
		# choose bcr date 
		"""
		Another thing I am noticing is that some patients do not meet the definition of BCR because they are treated before a second confirmatory test i.e. 
		patient has undetectable PSAs and then it goes to 0.6 and then the next value is undetectable. I've cut receipt of hormonal therapy (in processed data, under meds) 
		and the first date of that can be used as an event indicator (hormonepost_rp=1). Considering that some patients are being treated before BCR, 
		we should consider rephrasing our outcome to "failure-free survival" with events being first of BCR (PSA>0.2 on two occasions, 
		receipt of hormonal therapy, or radiation). (I am still cutting the radiation treatment).

		bcr_date is reported min(date of detactable PSA between undetectable and detectable tests, latter date of two consecutive PSA dates)
		"""
		if mode == 'rp':
			if df_outcome.at[empi, 'first_post_op_psa_valid'] == 1:
				for idx, psa_level, psa_date in zip(np.arange(len(psa_levels_after_rp)), psa_levels_after_rp, psa_levels_after_rp_dates):
					if not record_bcr_event:
						if idx < len(psa_levels_after_rp) -2:
							if psa_level <= lab_val_thresh and psa_levels_after_rp[idx + 1] > lab_val_thresh and psa_levels_after_rp[idx + 2] <= lab_val_thresh:
								df_outcome.at[empi, 'hormone_post_op'] = 1
								df_outcome.at[empi, 'bcr_date'] = psa_levels_after_rp_dates[idx + 1]
								df_outcome.at[empi, 'bcr_post_op'] = 1
								df_outcome.at[empi, 'bcr_date_minus_rp_date_in_days'] = (psa_levels_after_rp_dates[idx + 1] - RP_date).days
								record_bcr_event = 1
								break
							elif psa_level > lab_val_thresh and psa_levels_after_rp[idx + 1] > lab_val_thresh:
								df_outcome.at[empi, 'bcr_date'] = psa_levels_after_rp_dates[idx + 1] # get confirmatory date
								df_outcome.at[empi, 'bcr_post_op'] = 1
								df_outcome.at[empi, 'bcr_date_minus_rp_date_in_days'] = (psa_levels_after_rp_dates[idx + 1] - RP_date).days
								record_bcr_event = 1
								break
							elif psa_level > lab_val_thresh + 0.3: # if psa level > 0.5, then no need for a subsequent confirmatory test
								df_outcome.at[empi, 'bcr_date'] = psa_date
								df_outcome.at[empi, 'bcr_post_op'] = 1
								df_outcome.at[empi, 'bcr_date_minus_rp_date_in_days'] = (psa_date - RP_date).days
								record_bcr_event = 1
								break
						elif idx < len(psa_levels_after_rp) -1:
							if psa_level > lab_val_thresh and psa_levels_after_rp[idx + 1] > lab_val_thresh:
								df_outcome.at[empi, 'bcr_date'] = psa_levels_after_rp_dates[idx + 1] # get confirmatory date
								df_outcome.at[empi, 'bcr_post_op'] = 1
								df_outcome.at[empi, 'bcr_date_minus_rp_date_in_days'] = (psa_levels_after_rp_dates[idx + 1] - RP_date).days
								record_bcr_event = 1
								break
							elif psa_level > lab_val_thresh + 0.3: # if psa level > 0.5, then no need for a subsequent confirmatory test
								df_outcome.at[empi, 'bcr_date'] = psa_date
								df_outcome.at[empi, 'bcr_post_op'] = 1
								df_outcome.at[empi, 'bcr_date_minus_rp_date_in_days'] = (psa_date - RP_date).days
								record_bcr_event = 1
								break
						# else: # edge-case at the end of the trajectory
						# 	if psa_level > 0.2 and psa_levels_after_rp[idx -1] > 0.2:
						# 		df_outcome.at[empi, 'bcr_date'] = psa_date
						# 		df_outcome.at[empi, 'bcr_post_op'] = 1
						# 		df_outcome.at[empi, 'bcr_date_minus_rp_date_in_days'] = (psa_date - RP_date).days
						# 		record_bcr_event = 1
						# 		break

				if not record_bcr_event:
					df_outcome.at[empi, 'bcr_post_op'] = 0		
					df_outcome.at[empi, 'bcr_date']	= psa_levels_after_rp_dates[-1]	
					df_outcome.at[empi, 'bcr_date_minus_rp_date_in_days'] = (psa_levels_after_rp_dates[-1] - RP_date).days
		else: # if a patient has the valid nadir PSA
			if not pd.isnull(df_outcome.at[empi, 'nadir_psa_v2']):
				for idx, psa_level, psa_date in zip(np.arange(len(psa_levels_after_rp)), psa_levels_after_rp, psa_levels_after_rp_dates):
					if psa_date > df_outcome.at[empi, 'nadir_psa_date_v2']:
						if not record_bcr_event:
							if psa_level >= df_outcome.at[empi, 'nadir_psa_v2'] + 2:
								df_outcome.at[empi, 'bcr_date'] = psa_date
								df_outcome.at[empi, 'bcr_post_rad'] = 1
								df_outcome.at[empi, 'bcr_date_minus_radiation_date_in_days'] = (psa_date - RP_date).days
								record_bcr_event = 1
								break
							# if idx < len(psa_levels_after_rp) -2:
							# 	if psa_level <= lab_val_thresh and psa_levels_after_rp[idx + 1] > lab_val_thresh and psa_levels_after_rp[idx + 2] <= lab_val_thresh:
							# 		df_outcome.at[empi, 'hormone_post_op'] = 1
							# 		df_outcome.at[empi, 'bcr_date'] = psa_levels_after_rp_dates[idx + 1]
							# 		df_outcome.at[empi, 'bcr_post_rad'] = 1
							# 		df_outcome.at[empi, 'bcr_date_minus_rp_date_in_days'] = (psa_levels_after_rp_dates[idx + 1] - RP_date).days
							# 		record_bcr_event = 1
							# 		break
							# 	elif psa_level > lab_val_thresh and psa_levels_after_rp[idx + 1] > lab_val_thresh:
							# 		df_outcome.at[empi, 'bcr_date'] = psa_levels_after_rp_dates[idx + 1] # get confirmatory date
							# 		df_outcome.at[empi, 'bcr_post_op'] = 1
							# 		df_outcome.at[empi, 'bcr_date_minus_rp_date_in_days'] = (psa_levels_after_rp_dates[idx + 1] - RP_date).days
							# 		record_bcr_event = 1
							# 		break
							# elif idx < len(psa_levels_after_rp) -1:
							# 	if psa_level > lab_val_thresh and psa_levels_after_rp[idx + 1] > lab_val_thresh:
							# 		df_outcome.at[empi, 'bcr_date'] = psa_levels_after_rp_dates[idx + 1] # get confirmatory date
							# 		df_outcome.at[empi, 'bcr_post_op'] = 1
							# 		df_outcome.at[empi, 'bcr_date_minus_rp_date_in_days'] = (psa_levels_after_rp_dates[idx + 1] - RP_date).days
							# 		record_bcr_event = 1
							# 		break
						# else: # edge-case at the end of the trajectory
						# 	if psa_level > 0.2 and psa_levels_after_rp[idx -1] > 0.2:
						# 		df_outcome.at[empi, 'bcr_date'] = psa_date
						# 		df_outcome.at[empi, 'bcr_post_op'] = 1
						# 		df_outcome.at[empi, 'bcr_date_minus_rp_date_in_days'] = (psa_date - RP_date).days
						# 		record_bcr_event = 1
						# 		break

				if not record_bcr_event:
					df_outcome.at[empi, 'bcr_post_rad'] = 0		
					df_outcome.at[empi, 'censor_date']	= psa_levels_after_rp_dates[-1]	
					df_outcome.at[empi, 'censor_date_minus_radiation_date_in_days'] = (psa_levels_after_rp_dates[-1] - RP_date).days
					# possibly 

	if export_dic:
		f = open("../data/empi_to_last_psa_check_date_dic.pkl", "wb")
		pickle.dump(empi_to_last_psa_check_date_dic,f)
		f.close()
	# breakpoint()
	return df_outcome

def merge_lab_data():
	# df_lab_results_1 = pd.read_csv('../data/KS185_20200918_114153_Lab.txt', sep = '|', error_bad_lines = False)
	df_lab_results_1 = pd.read_csv('../data/lab_1.csv')
	df_lab_results_2 = pd.read_csv('../data/lab_2.csv')
	df_lab_results_3 = pd.read_csv('../data/lab_3.csv')
	df_lab_results_4 = pd.read_csv('../data/lab_4.csv')

	cols_oi = ['EMPI', 'Seq_Date_Time', 'Group_Id', 'Test_Id', 'Result', 'Result_Text', 'Reference_Units', 'Reference_Range']
	cols_oi_others = ['ID_MERGE', 'time_lab_result', 'lab_group', 'lab_testID', 'lab_result', 'lab_result_txt', 'lab_result_unit', 'lab_result_range']
	
	df_lab_results_1_cols = df_lab_results_1[cols_oi_others]
	df_lab_results_1_cols.columns = cols_oi

	df_lab_results_2_cols = df_lab_results_2[cols_oi_others]
	df_lab_results_2_cols.columns = cols_oi

	df_lab_results_3_cols = df_lab_results_3[cols_oi_others]
	df_lab_results_3_cols.columns = cols_oi

	df_lab_results_4_cols = df_lab_results_4[cols_oi_others]
	df_lab_results_4_cols.columns = cols_oi

	df_lab_total = pd.concat([df_lab_results_1_cols, df_lab_results_2_cols, df_lab_results_3_cols, df_lab_results_4_cols])
	df_lab_total.Seq_Date_Time = pd.to_datetime(df_lab_total.Seq_Date_Time.values)

	df_lab_total.index = np.arange(len(df_lab_total))
	return df_lab_total

def process_bmi(df_bmi, empi_to_date_oi_dic, empis_to_use, days_after_rp_margin = 5):
	unique_empis_bmi = np.unique(df_bmi.EMPI.values)
	df_bmi_oi = df_bmi.loc[df_bmi.EMPI.isin(empis_to_use)]
	# empis_oi = empis_to_use & set(df_bmi.EMPI.values)
	df_bmi_total = pd.DataFrame([])
	for empi in tqdm(np.unique(df_bmi_oi.EMPI.values), desc = 'getting the most recent bmi before RP...'):
		df_bmi_oi_per_empi = df_bmi_oi.loc[df_bmi_oi.EMPI == empi]
		df_bmi_oi_per_empi.sort_values(by = 'diadate_parsed', ascending = False, inplace = True)
		if empi in empi_to_date_oi_dic.keys():
			rp_date = empi_to_date_oi_dic[empi]
			df_bmi_oi_per_empi.diadate_parsed = pd.to_datetime(df_bmi_oi_per_empi.diadate_parsed.values)
			df_bmi_oi_per_empi_filtered = df_bmi_oi_per_empi.loc[df_bmi_oi_per_empi.diadate_parsed <= rp_date + pd.DateOffset(days = days_after_rp_margin)]
			if len(df_bmi_oi_per_empi_filtered) > 0:
				df_bmi_total = pd.concat([df_bmi_total, df_bmi_oi_per_empi_filtered.iloc[0:1]])
	return df_bmi_total 

def process_comorbidity(df_comorb, empi_to_date_oi_dic, empis_to_use):
	unique_empis_comorb = np.unique(df_comorb.EMPI.values)
	cols_oi = ['EMPI', 'Date'] + [col for col in df_comorb.columns if '_agg' in col]
	
	# df_comorb_oi = df_comorb[cols_oi]
	# empis_oi = empis_to_use & set(df_comorb.EMPI.values)
	df_comorb_oi = df_comorb[cols_oi]
	df_comorb_oi = df_comorb_oi.loc[df_comorb_oi.EMPI.isin(empis_to_use)]
	df_comorb_total = pd.DataFrame([])
	for empi in tqdm(np.unique(df_comorb_oi.EMPI.values), desc = 'getting the most recent comorbidity before RP...'):
		df_comorb_oi_per_empi = df_comorb_oi.loc[df_comorb_oi.EMPI == empi].copy()
		df_comorb_oi_per_empi.sort_values(by = 'Date', ascending = False, inplace = True)
		if empi in empi_to_date_oi_dic.keys():
			date_oi = empi_to_date_oi_dic[empi]
			df_comorb_oi_per_empi.Date = pd.to_datetime(df_comorb_oi_per_empi.Date.values)
			df_comorb_oi_per_empi_filtered = df_comorb_oi_per_empi.loc[df_comorb_oi_per_empi.Date <= date_oi]
			# breakpoint()
			if len(df_comorb_oi_per_empi_filtered) > 0:
				df_comorb_total = pd.concat([df_comorb_total, df_comorb_oi_per_empi_filtered.iloc[0:1]])
		# else:
		# 	df_comorb_total = pd.concat([df_comorb_total, df_comorb_oi_per_empi.iloc[0:1]])
	return df_comorb_total

def merge_dfs(dfs_list_dic, df_type_to_cols_oi_dic, empis_common):
	# df_merged = pd.DataFrame([], index = empis_common)
	counter = 0
	for df_type, df in dfs_list_dic.items():
		df_oi = df.loc[empis_common][df_type_to_cols_oi_dic[df_type]]
		if counter == 0:
			df_merged = df_oi
		else:
			df_merged = pd.concat([df_merged, df_oi], axis = 1)
		counter += 1
	return df_merged

def main(data_comb, load_processed = False, get_biopsy_based_feature = False):
	# load processed files
	print('Loading and processing...')
	# df_height = pd.read_csv('../data/processed_data/height.csv')
	# df_weight = pd.read_csv('../data/processed_data/processed_weights.csv')
	df_race = pd.read_csv('../data/processed_data/race.csv')
	# convert race into one-hot-encoding
	df_race_ont_hot = pd.DataFrame([], index = df_race.EMPI.values, columns = ['White', 'Black', 'Asian', 'Hispanic', 'Unknown/other'])
	for empi, race in zip(df_race.EMPI.values, df_race.race_cat.values):
		if race == 'White':
			df_race_ont_hot.at[empi, 'White'] = 1
			df_race_ont_hot.at[empi, 'Black'] = 0
			df_race_ont_hot.at[empi, 'Asian'] = 0
			df_race_ont_hot.at[empi, 'Hispanic'] = 0
			df_race_ont_hot.at[empi, 'Unknown/other'] = 0
		elif race == 'Black':
			df_race_ont_hot.at[empi, 'White'] = 0
			df_race_ont_hot.at[empi, 'Black'] = 1
			df_race_ont_hot.at[empi, 'Asian'] = 0
			df_race_ont_hot.at[empi, 'Hispanic'] = 0
			df_race_ont_hot.at[empi, 'Unknown/other'] = 0
		elif race == 'Asian':
			df_race_ont_hot.at[empi, 'White'] = 0
			df_race_ont_hot.at[empi, 'Black'] = 0
			df_race_ont_hot.at[empi, 'Asian'] = 1
			df_race_ont_hot.at[empi, 'Hispanic'] = 0
			df_race_ont_hot.at[empi, 'Unknown/other'] = 0
		elif race == 'Hispanic':
			df_race_ont_hot.at[empi, 'White'] = 0
			df_race_ont_hot.at[empi, 'Black'] = 0
			df_race_ont_hot.at[empi, 'Asian'] = 0
			df_race_ont_hot.at[empi, 'Hispanic'] = 1
			df_race_ont_hot.at[empi, 'Unknown/other'] = 0
		else:
			df_race_ont_hot.at[empi, 'White'] = 0
			df_race_ont_hot.at[empi, 'Black'] = 0
			df_race_ont_hot.at[empi, 'Asian'] = 0
			df_race_ont_hot.at[empi, 'Hispanic'] = 0
			df_race_ont_hot.at[empi, 'Unknown/other'] = 1


	# single radical prostatectomy patients from procedure data
	single_rp_patients = pd.read_csv('../data/singlerp.csv')
	# empi_to_rp_date_dic = {empi:pd.to_datetime(date) for empi, date in zip(single_rp_patients.EMPI.values, single_rp_patients.prdate_parsed.values)}

	# f = open("../data/empi_to_rp_date_dic.pkl", "wb")
	# pickle.dump(empi_to_rp_date_dic,f)
	# f.close()

	with open('../data/empi_to_rp_date_dic.pkl', 'rb') as handle:
		empi_to_rp_date_dic = pickle.load(handle)

	# get empi to RP date dic
	# empi_to_date_oi_dic = {empi:pd.to_datetime(date) for empi, date in zip(single_rp_patients.EMPI.values, single_rp_patients.prdate_parsed.values)}
	# empi_to_biopsy_report_date_dic = {empi:pd.to_datetime(date) for empi, date in zip(df_pathology_biopsy_final.EMPI.values, df_pathology_biopsy_final.Report_Date_Time.values)}
	# df_lab_results = pd.read_csv('../data/KS185_20200918_114153_Lab.txt', sep = '|', error_bad_lines = False)

	# radical prostatectomy pathology report 

	
	# df_pathology_with_biopsy_all_feats.Report_Date_Time = pd.to_datetime(df_pathology_with_biopsy_all_feats.Report_Date_Time.values)
	# df_pathology_with_biopsy_all_feats.drop(columns = ['Unnamed: 0', 'Unnamed: 0.1'], inplace = True)
	# breakpoint()

	# empis of interest from rp pathology
	if load_processed:
		if get_biopsy_based_feature:
			df_pathology_biopsy_final = pd.read_csv('../data/df_pathology_biopsy_final.csv')
			# if load_processed:
			df_pathology_biopsy_final.set_index('EMPI', inplace = True)
			# empis_oi = set(df_pathology_biopsy_final.EMPI.values)
			df_labs_processed = pd.read_csv('../data/df_labs_processed_biopsy_based.csv')
			df_labs_processed.rename(columns = {'Unnamed: 0': 'EMPI'}, inplace = True)
			df_labs_processed.set_index('EMPI', inplace = True)

			with open('../data/empi_to_prior_psa_levels_dic_biopsy_based.pkl', 'rb') as handle:
				empi_to_prior_psa_levels_dic = pickle.load(handle)
			df_psa_prior_to_rp = pd.DataFrame([], index = empi_to_prior_psa_levels_dic.keys(), columns = ['psa_trajectory', 'max_psa', 'min_psa', 'mean_psa'])
			df_psa_prior_to_rp.index.name = 'EMPI'
			for empi, psa_list in empi_to_prior_psa_levels_dic.items():
				df_psa_prior_to_rp.at[empi, 'psa_trajectory'] = psa_list
				df_psa_prior_to_rp.at[empi, 'max_psa'] = np.max(psa_list)
				df_psa_prior_to_rp.at[empi, 'min_psa'] = np.min(psa_list)
				df_psa_prior_to_rp.at[empi, 'mean_psa'] = np.mean(psa_list)
				df_psa_prior_to_rp.at[empi, 'psa_prior_to_rp'] = psa_list[-1]

			df_comorb_processed = pd.read_csv('../data/df_comorb_processed_biopsy_based.csv')
			df_comorb_processed.drop(columns = ['Unnamed: 0'], inplace = True)

			df_bmi_processed = pd.read_csv('../data/df_bmi_processed_biopsy_based.csv')
			df_bmi_processed.drop(columns = ['Unnamed: 0'], inplace = True)

			df_age_at_rp = pd.read_csv('../data/df_age_at_rp_biopsy_based.csv')
			df_age_at_rp.set_index('EMPI', inplace = True)
		else:
			df_pathology_with_rp_all_feats = pd.read_csv('df_pathology_with_grade_pt_stage_margin.tsv', sep = '\t')
			# empis_oi = set(df_pathology_with_rp_all_feats.EMPI.values)
			df_labs_processed = pd.read_csv('../data/df_labs_processed.csv')
			df_labs_processed.rename(columns = {'Unnamed: 0': 'EMPI'}, inplace = True)
			df_labs_processed.set_index('EMPI', inplace = True)

			with open('../data/empi_to_prior_psa_levels_dic.pkl', 'rb') as handle:
				empi_to_prior_psa_levels_dic = pickle.load(handle)
			df_psa_prior_to_rp = pd.DataFrame([], index = empi_to_prior_psa_levels_dic.keys(), columns = ['psa_trajectory', 'max_psa', 'min_psa', 'mean_psa'])
			df_psa_prior_to_rp.index.name = 'EMPI'
			for empi, psa_list in empi_to_prior_psa_levels_dic.items():
				df_psa_prior_to_rp.at[empi, 'psa_trajectory'] = psa_list
				df_psa_prior_to_rp.at[empi, 'max_psa'] = np.max(psa_list)
				df_psa_prior_to_rp.at[empi, 'min_psa'] = np.min(psa_list)
				df_psa_prior_to_rp.at[empi, 'mean_psa'] = np.mean(psa_list)
				df_psa_prior_to_rp.at[empi, 'psa_prior_to_rp'] = psa_list[-1]

			df_comorb_processed = pd.read_csv('../data/df_comorb_processed.csv')
			df_comorb_processed.drop(columns = ['Unnamed: 0'], inplace = True)

			df_bmi_processed = pd.read_csv('../data/df_bmi_processed.csv')
			df_bmi_processed.drop(columns = ['Unnamed: 0'], inplace = True)

			df_age_at_rp = pd.read_csv('../data/df_age_at_rp.csv')
			df_age_at_rp.set_index('EMPI', inplace = True)
	else:
		# biopsy paramters for both positive and negative RP patients
		if get_biopsy_based_feature:
			df_pathology_with_biopsy_all_feats = pd.read_csv('../data/df_biopsy_wo_duplicate_processed.csv')
			empis_rp_in_biopsy_reports = set(df_pathology_with_biopsy_all_feats.EMPI.values) & set(single_rp_patients.EMPI.values)
			df_pathology_with_biopsy_all_feats.set_index('EMPI', inplace = True)
			df_pathology_with_biopsy_all_feats['rp_indicator'] = 0
			df_pathology_with_biopsy_all_feats.loc[empis_rp_in_biopsy_reports, 'rp_indicator'] = 1
			df_pathology_with_biopsy_all_feats['rp_date'] = None
			df_pathology_with_biopsy_all_feats['rp_date_minus_report_date_in_days'] = None
			for empi in df_pathology_with_biopsy_all_feats.index:
				if empi in empis_rp_in_biopsy_reports:
					df_pathology_with_biopsy_all_feats.at[empi, 'rp_date'] = empi_to_rp_date_dic[empi]
					df_pathology_with_biopsy_all_feats.at[empi, 'rp_date_minus_report_date_in_days'] = (empi_to_rp_date_dic[empi] - pd.to_datetime(df_pathology_with_biopsy_all_feats.at[empi, 'Report_Date_Time'])).days
			df_pathology_with_biopsy_all_feats_valid_rp_positive = df_pathology_with_biopsy_all_feats.loc[df_pathology_with_biopsy_all_feats.rp_date_minus_report_date_in_days >= 0]
			df_pathology_with_biopsy_all_feats_valid_rp_negative = df_pathology_with_biopsy_all_feats.loc[df_pathology_with_biopsy_all_feats.rp_indicator == 0]
			df_pathology_biopsy_final = pd.concat([df_pathology_with_biopsy_all_feats_valid_rp_positive, df_pathology_with_biopsy_all_feats_valid_rp_negative])
			df_pathology_biopsy_final.to_csv('../data/df_pathology_biopsy_final.csv')
			# df_pathology_biopsy_final.set_index('EPIC_PMRN')
			empis_oi = set(df_pathology_biopsy_final.index)

			# get mapping from empi to biopsy report date dic
			empi_to_date_oi_dic = {empi:pd.to_datetime(date) for empi, date in zip(df_pathology_biopsy_final.index, df_pathology_biopsy_final.Report_Date_Time.values)}
			f = open("../data/empi_to_date_oi_dic.pkl", "wb")
			pickle.dump(empi_to_date_oi_dic,f)
			f.close()

			# breakpoint()

			# get lab values
			df_labs = pd.read_csv('../data/labs_merged.csv')
			cols_oi = ['EMPI', 'Seq_Date_Time', 'Group_Id', 'Test_Id', 'Result', 'Result_Text', 'Reference_Units', 'Reference_Range']
			cols_oi_others = ['ID_MERGE', 'time_lab_result', 'lab_group', 'lab_testID', 'lab_result', 'lab_result_txt', 'lab_result_unit', 'lab_result_range']
			df_labs = df_labs[cols_oi_others]
			df_labs.columns = cols_oi
			# df_labs = pd.read_csv('../data/lab_total.csv')
			# df_labs.drop(columns = ['Unnamed: 0'], inplace = True)

			df_labs_processed = get_psa_prior_to_rp(df_labs, empi_to_date_oi_dic, empis_oi, export_dic = True)#, get_biopsy_based_feature = get_biopsy_based_feature)
			# breakpoint()
			# post-fix : fix rp date to biopsy date and add actual rp date
			df_labs_processed.rename(columns = {'rp_date':'biopsy_report_date'}, inplace = True)
			df_labs_processed['rp_date'] = None
			df_labs_processed['rp_date_minus_biopsy_date_in_days'] = None
			for empi in df_labs_processed.index:
				if empi in empi_to_rp_date_dic.keys():
					df_labs_processed.at[empi, 'rp_date'] = empi_to_rp_date_dic[empi]
					df_labs_processed.at[empi, 'rp_date_minus_biopsy_date_in_days'] = (empi_to_rp_date_dic[empi] - df_labs_processed.at[empi, 'biopsy_report_date']).days
			df_labs_processed.to_csv('../data/df_labs_processed_biopsy_based.csv')

			empi_to_prior_psa_levels_dic = {empi : [float(val) for val in prior_psa_list] for empi, prior_psa_list in zip(df_labs_processed.index, df_labs_processed['pre_op_psa_trajectory_10_year_window'].values) if len(prior_psa_list) > 0}
			
			f = open("../data/empi_to_prior_psa_levels_dic_biopsy_based.pkl", "wb")
			pickle.dump(empi_to_prior_psa_levels_dic,f)
			f.close()

			# breakpoint()
			df_psa_prior_to_rp = pd.DataFrame([], index = empi_to_prior_psa_levels_dic.keys(), columns = ['psa_trajectory', 'max_psa', 'min_psa', 'mean_psa'])
			df_psa_prior_to_rp.index.name = 'EMPI'
			for empi, psa_list in empi_to_prior_psa_levels_dic.items():
				df_psa_prior_to_rp.at[empi, 'psa_trajectory'] = psa_list
				df_psa_prior_to_rp.at[empi, 'max_psa'] = np.max(psa_list)
				df_psa_prior_to_rp.at[empi, 'min_psa'] = np.min(psa_list)
				df_psa_prior_to_rp.at[empi, 'mean_psa'] = np.mean(psa_list)
			# breakpoint()

			# get comorbidity values prior to RP
			df_comorb = pd.read_csv('../data/processed_data/dia_comorbDL_total.csv')
			df_comorb_processed = process_comorbidity(df_comorb, empi_to_date_oi_dic, empis_oi)
			df_comorb_processed.to_csv('../data/df_comorb_processed_biopsy_based.csv')
			# breakpoint()

			# get bmi values prior to RP
			df_bmi = pd.read_csv('../data/processed_data/bmi.csv')
			df_bmi_processed = process_bmi(df_bmi, empi_to_date_oi_dic, empis_oi)
			df_bmi_processed.to_csv('../data/df_bmi_processed_biopsy_based.csv')

			# get age at RP
			df_birth_death_date = pd.read_csv('../data/processed_data/birthdeath.csv')
			empi_to_death_date_dic = {empi:pd.to_datetime(death_date) for empi, death_date in zip(df_birth_death.EMPI.values, df_birth_death.dth_date_parsed.values) if not pd.isnull(death_date)}
			f = open("../data/empi_to_death_date_dic.pkl", "wb")
			pickle.dump(empi_to_death_date_dic,f)
			f.close()
			empi_to_birth_date_dic = {empi:pd.to_datetime(date) for empi, date in zip(df_birth_death_date.EMPI.values, df_birth_death_date.birth_date_parsed.values)}

			empis_oi_age = set(empi_to_date_oi_dic.keys()) & set(empi_to_birth_date_dic.keys())
			df_age_at_rp = pd.DataFrame(0.0, index = empis_oi_age, columns = ['Age at RP'])
			for empi in df_age_at_rp.index:
				df_age_at_rp.at[empi, 'Age at RP'] = int((empi_to_date_oi_dic[empi] - empi_to_birth_date_dic[empi]).days/365.2425)
			df_age_at_rp.index.name = 'EMPI'
			df_age_at_rp.to_csv('../data/df_age_at_rp_biopsy_based.csv')
		else:
			df_pathology_with_rp_all_feats = pd.read_csv('df_pathology_with_grade_pt_stage_margin.tsv', sep = '\t')
			empis_oi = set(df_pathology_with_rp_all_feats.EMPI.values)# | set(df_pathology_biopsy_final.EMPI.values)
			# rp datea and empi mapping
			empi_to_date_oi_dic = {empi:pd.to_datetime(date) for empi, date in zip(single_rp_patients.EMPI.values, single_rp_patients.prdate_parsed.values)}
			# get PSA values prior to RP
			# df_labs = merge_lab_data()
			# df_labs.to_csv('../data/lab_total.csv')
			df_labs = pd.read_csv('../data/labs_merged.csv')
			cols_oi = ['EMPI', 'Seq_Date_Time', 'Group_Id', 'Test_Id', 'Result', 'Result_Text', 'Reference_Units', 'Reference_Range']
			cols_oi_others = ['ID_MERGE', 'time_lab_result', 'lab_group', 'lab_testID', 'lab_result', 'lab_result_txt', 'lab_result_unit', 'lab_result_range']
			df_labs = df_labs[cols_oi_others]
			df_labs.columns = cols_oi
			# df_labs = pd.read_csv('../data/lab_total.csv')
			# df_labs.drop(columns = ['Unnamed: 0'], inplace = True)

			df_labs_processed = get_psa_prior_to_rp(df_labs, empi_to_date_oi_dic, empis_oi)
			df_labs_processed.to_csv('../data/df_labs_processed.csv')
			# breakpoint()

			empi_to_prior_psa_levels_dic = {empi : [float(val) for val in prior_psa_list] for empi, prior_psa_list in zip(df_labs_processed.index, df_labs_processed['pre_op_psa_trajectory_10_year_window'].values) if len(prior_psa_list) > 0}
			
			f = open("../data/empi_to_prior_psa_levels_dic.pkl", "wb")
			pickle.dump(empi_to_prior_psa_levels_dic,f)
			f.close()
			# breakpoint()
			df_psa_prior_to_rp = pd.DataFrame([], index = empi_to_prior_psa_levels_dic.keys(), columns = ['psa_trajectory', 'max_psa', 'min_psa', 'mean_psa'])
			df_psa_prior_to_rp.index.name = 'EMPI'
			for empi, psa_list in empi_to_prior_psa_levels_dic.items():
				df_psa_prior_to_rp.at[empi, 'psa_trajectory'] = psa_list
				df_psa_prior_to_rp.at[empi, 'max_psa'] = np.max(psa_list)
				df_psa_prior_to_rp.at[empi, 'min_psa'] = np.min(psa_list)
				df_psa_prior_to_rp.at[empi, 'mean_psa'] = np.mean(psa_list)

			# get comorbidity values prior to RP
			df_comorb = pd.read_csv('../data/processed_data/dia_comorbDL.csv')
			df_comorb_processed = process_comorbidity(df_comorb, empi_to_date_oi_dic, empis_oi)
			df_comorb_processed.to_csv('../data/df_comorb_processed.csv')

			# get bmi values prior to RP
			df_bmi = pd.read_csv('../data/processed_data/bmi.csv')
			df_bmi_processed = process_bmi(df_bmi, empi_to_date_oi_dic, empis_oi)
			df_bmi_processed.to_csv('../data/df_bmi_processed.csv')

			# get age at RP
			df_birth_death_date = pd.read_csv('../data/processed_data/birthdeath.csv')
			empi_to_birth_date_dic = {empi:pd.to_datetime(date) for empi, date in zip(df_birth_death_date.EMPI.values, df_birth_death_date.birth_date_parsed.values)}
			empi_to_death_date_dic = {empi:pd.to_datetime(death_date) for empi, death_date in zip(df_birth_death.EMPI.values, df_birth_death.dth_date_parsed.values) if not pd.isnull(death_date)}
			f = open("../data/empi_to_death_date_dic.pkl", "wb")
			pickle.dump(empi_to_death_date_dic,f)
			f.close()

			empis_oi_age = set(empi_to_date_oi_dic.keys()) & set(empi_to_birth_date_dic.keys())
			df_age_at_rp = pd.DataFrame(0.0, index = empis_oi_age, columns = ['Age at RP'])
			for empi in df_age_at_rp.index:
				df_age_at_rp.at[empi, 'Age at RP'] = int((empi_to_date_oi_dic[empi] - empi_to_birth_date_dic[empi]).days/365.2425)
			df_age_at_rp.index.name = 'EMPI'
			df_age_at_rp.to_csv('../data/df_age_at_rp.csv')

	print('Processing/Loading complete!')
	# breakpoint()
	
	# empis_oi = set(df_pathology_with_rp_all_feats.EMPI.values) & set(df_lab_total.EMPI.values)
	# empis_oi_biopsies = set(df_pathology_with_biopsy_all_feats.EMPI.values) & set(df_lab_total.EMPI.values)
	

	# merge bmi, race, age, comorb
	if get_biopsy_based_feature:
		df_pathology_oi = df_pathology_biopsy_final[['Report_Date_Time', 'overall_grade_merged', 'benign', 'num_pos_cores_sum', 'auxiiliary_mci_score', 'rp_indicator']]
		df_type_to_cols_oi_dic = {'pat': ['overall_grade_merged', 'benign', 'num_pos_cores_sum', 'auxiiliary_mci_score', 'rp_indicator'], 'comorb': ['ami_agg', 'chf_agg', 'pvd_agg', 'cevd_agg', 'copd_agg', 'rheumd_agg', 'pud_agg', 'mld_agg', 'diab_agg', 'diabwc_agg', 'hp_agg', 'rend_agg', 'metacanc_agg', 'aids_agg', 'wscore_agg'], 'age' : ['Age at RP'], 'race' : ['White', 'Black', 'Asian', 'Hispanic', 'Unknown/other'], 'bmi' : ['BMI'], 'psa' : ['max_psa', 'min_psa', 'mean_psa', 'psa_prior_to_rp']}
	else:
		df_pathology_oi = df_pathology_with_rp_all_feats[['EMPI', 'Report_Date_Time', 'overall_grade_group', 'pT_stage_combined', 'margin']]
		# drop patients with more than one pathology report
		df_pathology_oi.drop_duplicates(subset = 'EMPI', keep = False, inplace = True)
		df_pathology_oi.set_index('EMPI', inplace = True)
		df_type_to_cols_oi_dic = {'pat': ['overall_grade_group', 'pT_stage_combined', 'margin'], 'comorb': ['ami_agg', 'chf_agg', 'pvd_agg', 'cevd_agg', 'copd_agg', 'rheumd_agg', 'pud_agg', 'mld_agg', 'diab_agg', 'diabwc_agg', 'hp_agg', 'rend_agg', 'metacanc_agg', 'aids_agg', 'wscore_agg'], 'age' : ['Age at RP'], 'race' : ['White', 'Black', 'Asian', 'Hispanic', 'Unknown/other'], 'bmi' : ['BMI'], 'psa' : ['max_psa', 'min_psa', 'mean_psa', 'psa_prior_to_rp']}
	
	# set EMPI as index
	df_comorb_processed.set_index('EMPI', inplace = True)
	df_bmi_processed.set_index('EMPI', inplace = True)

	dfs_list_dic = {'pat' : df_pathology_oi, 'comorb' : df_comorb_processed, 'age': df_age_at_rp, 'race':df_race_ont_hot, 'bmi':df_bmi_processed, 'psa':df_psa_prior_to_rp}
	empis_common = set(df_pathology_oi.index) & set(df_comorb_processed.index) & set(df_age_at_rp.index) & set(df_race_ont_hot.index) & set(df_bmi_processed.index) & set(df_psa_prior_to_rp.index) # need to imput BMI
	# A : Pathology & Age & Race & Comorbidity & BMI & PSA prior (n = 1454)
	# B : Pathology & Age & Race (n = 9778)
	# C : Pathology & Age & Race & Comorbidity (n = 7032)
	# D : Pathology & Age & Race & BMI (n = 5000)
	# E : Pathology & Age & Race & PSA prior (n = 3676)
	print('Total number of unique patients in (A) : ', len(empis_common))

	if data_comb == 'B':
		remove_list = ['comorb', 'bmi', 'psa']
		[dfs_list_dic.pop(key) for key in remove_list]
		[df_type_to_cols_oi_dic.pop(key) for key in remove_list]
		empis_common = set(df_pathology_oi.index) & set(df_age_at_rp.index) & set(df_race_ont_hot.index)
		print('Total number of unique patients in (B) : ', len(empis_common))
	elif data_comb == 'C':
		remove_list = ['bmi', 'psa']
		[dfs_list_dic.pop(key) for key in remove_list]
		[df_type_to_cols_oi_dic.pop(key) for key in remove_list]
		empis_common = set(df_pathology_oi.index) & set(df_comorb_processed.index) & set(df_age_at_rp.index) & set(df_race_ont_hot.index)
		print('Total number of unique patients in (C) : ', len(empis_common))
	elif data_comb == 'D':
		remove_list = ['comorb', 'psa']
		[dfs_list_dic.pop(key) for key in remove_list]
		[df_type_to_cols_oi_dic.pop(key) for key in remove_list]
		empis_common = set(df_pathology_oi.index) & set(df_bmi_processed.index) & set(df_age_at_rp.index) & set(df_race_ont_hot.index) 
		print('Total number of unique patients in (D) : ', len(empis_common))
	elif data_comb == 'E':
		remove_list = ['comorb', 'bmi']
		[dfs_list_dic.pop(key) for key in remove_list]
		[df_type_to_cols_oi_dic.pop(key) for key in remove_list]
		empis_common = set(df_pathology_oi.index) & set(df_psa_prior_to_rp.index) & set(df_age_at_rp.index) & set(df_race_ont_hot.index)
		print('Total number of unique patients in (E) : ', len(empis_common))

	
	# df_pathology_oi.set_index('EMPI', inplace = True)

	# check if there's any duplicate
	print('Ensuring there is no duplicate...')
	print('Length of df_pathology_oi : ', len(df_pathology_oi) , len(df_pathology_oi) == len(set(df_pathology_oi.index)))
	print('Length of df_comorb_processed : ', len(df_comorb_processed), len(df_comorb_processed) == len(set(df_comorb_processed.index)))
	print('Length of df_age_at_rp : ', len(df_age_at_rp), len(df_age_at_rp) == len(set(df_age_at_rp.index)))
	print('Length of df_bmi_processed : ', len(df_bmi_processed), len(df_bmi_processed) == len(set(df_bmi_processed.index)))
	print('Length of df_psa_prior_to_rp : ', len(df_psa_prior_to_rp), len(df_psa_prior_to_rp) == len(set(df_psa_prior_to_rp.index)))

	df_merged = merge_dfs(dfs_list_dic, df_type_to_cols_oi_dic, empis_common)
	if get_biopsy_based_feature:
		df_merged.to_csv('../data/merged/df_merged_biopsy_based_' + data_comb + '.csv')
	else:
		df_merged.to_csv('../data/merged/df_merged_rp_positive_based_' + data_comb + '.csv')
	get_outcome = True
	if get_outcome:
		with open('../data/empi_to_death_date_dic.pkl', 'rb') as handle:
			empi_to_death_date_dic = pickle.load(handle)
		# if get_biopsy_based_feature:
		# 	func_get_outcome(df_labs_processed_biopsy_based, empi_to_death_date_dic)
		# else:
		if get_biopsy_based_feature:
			with open('../data/empi_to_date_oi_dic.pkl', 'rb') as handle:
				empi_to_date_oi_dic = pickle.load(handle)
			with open('../data/empi_to_last_psa_check_date_dic_biopsy_based.pkl', 'rb') as handle:
				empi_to_last_psa_check_date_dic_biopsy_based = pickle.load(handle)
			clean_mode = input('Only use patients data before last NDI update (Nov-1-2017) (Y/N)? : ') == 'Y'
			extra_clean = True
			df_outcome_oi = func_get_outcome(df_pathology_biopsy_final, empi_to_death_date_dic, empi_to_date_oi_dic, get_biopsy_based_feature = get_biopsy_based_feature, empi_to_last_psa_check_date_dic_biopsy_based = empi_to_last_psa_check_date_dic_biopsy_based, clean_mode = clean_mode, extra_clean = extra_clean)
			if clean_mode:
				if extra_clean:
					df_outcome_oi.to_csv('../data/df_outcome_final_biopsy_based_extra_clean.csv')
				else:
					df_outcome_oi.to_csv('../data/df_outcome_final_biopsy_based_clean.csv')
			else:
				df_outcome_oi.to_csv('../data/df_outcome_final_biopsy_based.csv')
		else:
			df_outcome_oi = func_get_outcome(df_labs_processed, empi_to_death_date_dic, empi_to_rp_date_dic, get_biopsy_based_feature = get_biopsy_based_feature)
			# df_outcome_oi.to_csv()
	breakpoint()
	return


if __name__ == "__main__":
	explain = """
	Written by Intae Moon
	Date : May 6th, 2021

	Description :
	The following script merges data from different sources including pathology, comorbidity, age, race, bmi, and prior PSA values
	Type the letter for the combination of the data you want :
	
	Pathology data (rp pathology based : n = 9,781, biopsy based : n = 18,350)
	Age data (n = 18,345)
	Race data (n = 86,071)
	Comorbidity data (rp pathology based : n = 7,072, biopsy based : n = 10,542)
	BMI data  (rp pathology based : n = 4,625, biopsy based : n = 5,627)
	PSA prior (rp pathology based : n = 5,158, biopsy based : n = 6,042)

	A : Pathology & Age & Race & Comorbidity & BMI & PSA prior (rp pathology based : n = 2,327, biopsy based : n = 3,266)
	B : Pathology & Age & Race (rp pathology based : n = 9,778, biopsy based : n = 18,345)
	C : Pathology & Age & Race & Comorbidity (rp pathology based : n = 7,032, biopsy based : n = 10,539)
	D : Pathology & Age & Race & BMI (rp pathology based : n = 4,585, biopsy based : n = 5,620)
	E : Pathology & Age & Race & PSA prior (rp pathology based : n = 5,144, biopsy based : n = 6,027)

	"""
	print("\n")
	print(explain)
	print("\n")
	data_comb = input('Type the combination of datasets : ')
	if data_comb not in ['A', 'B', 'C', 'D', 'E']:
		raise KeyError('Choose it among A, B, C, D, and E')
	load_processed = input('Load the processed data? (Y/N) : ') == 'Y'
	get_biopsy_based_feature = input('Get biopsy-based data? (Y/N) : ') == 'Y'
	main(data_comb, load_processed = load_processed, get_biopsy_based_feature = get_biopsy_based_feature)