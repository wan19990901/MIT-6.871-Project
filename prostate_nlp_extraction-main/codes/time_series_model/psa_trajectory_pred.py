import pandas as pd
import numpy as np
from tqdm import tqdm
import re
from data_merger import get_psa_prior_to_rp
import pickle

def get_bcr_info(df_oi, dic_oi, mode = 'radiation'):
	if mode == 'radiation':
		for empi, bcr, days_oi_bcr, days_oi_censor, rad_date in zip(df_oi.index, df_oi.bcr_post_rad.values, 
														  			df_oi.bcr_date_minus_radiation_date_in_days.values, 
														  			df_oi.censor_date_minus_radiation_date_in_days.values,
														  			df_oi.radiation_date.values):
			if not pd.isnull(bcr):
				if bcr == 0:
					dic_oi[empi] = (bcr, days_oi_censor, rad_date, 'rad')
				elif bcr == 1:
					dic_oi[empi] = (bcr, days_oi_bcr, rad_date, 'rad')
	elif mode == 'rp':
		for empi, bcr, days_oi_bcr, rp_date in zip(df_oi.index, df_oi.bcr_post_op.values, 
												   df_oi.bcr_date_minus_rp_date_in_days.values, 
												   df_oi.rp_date.values):
			if not pd.isnull(bcr):
				dic_oi[empi] = (bcr, days_oi_bcr, rp_date, 'rp')
	else:
		raise KeyError('Only supports radiation or rp atm...')
	return dic_oi

def main():
	# load relevant data
	# patients with pathology biopsy with overall grade extracted
	df_pathology_with_biopsy_overall_grade = pd.read_csv('../data/df_pathology_with_biopsy_oi_overall_grade.csv')

	empis_to_bcr_info_dic = {}
	"""
	i) patients who underwent radiation therapy
	"""
	first_date_radiationt = pd.read_csv('../data/first_date_radiationtx.csv', index_col = 'EMPI')

	# labs data
	df_labs = pd.read_csv('../data/labs_merged.csv')
	cols_oi = ['EMPI', 'Seq_Date_Time', 'Group_Id', 'Test_Id', 'Result', 'Result_Text', 'Reference_Units', 'Reference_Range']
	cols_oi_others = ['ID_MERGE', 'time_lab_result', 'lab_group', 'lab_testID', 'lab_result', 'lab_result_txt', 'lab_result_unit', 'lab_result_range']
	df_labs = df_labs[cols_oi_others]
	df_labs.columns = cols_oi

	empi_to_date_oi_dic = {empi:pd.to_datetime(date) for empi, date in zip(first_date_radiationt.index, first_date_radiationt.rdtdate_first.values)}
	empis_oi = set(df_labs.EMPI.values) & set(first_date_radiationt.index) & set(df_pathology_with_biopsy_overall_grade.EMPI.values)
	# len(empis_oi)

	# process patients who received radiation
	print('Processing patients who received radiation therapy...')
	df_labs_processed_rad = get_psa_prior_to_rp(df_labs, empi_to_date_oi_dic, empis_oi, export_dic = False, mode = 'radiation')
	df_count_stats = pd.value_counts(df_labs_processed_rad.bcr_post_rad.values, sort = True)
	print(df_count_stats)
	print('BCR rates', np.round(df_count_stats[1]/df_count_stats.values.sum(), 3))
	# breakpoint()
	empis_to_bcr_info_dic = get_bcr_info(df_labs_processed_rad, empis_to_bcr_info_dic, mode = 'radiation')

	# drop those with invalid nadir psa
	df_labs_processed_rad_final = df_labs_processed_rad.dropna(subset = ['nadir_psa_v2'])

	'''
	ii) patients who underwent radical prostatectomy
	'''
	single_rp_patients = pd.read_csv('../data/singlerp.csv')
	df_pathology_with_rp_all_feats = pd.read_csv('df_pathology_with_grade_pt_stage_margin.tsv', sep = '\t')
	empis_oi = set(df_pathology_with_rp_all_feats.EMPI.values) & set(df_labs.EMPI.values)
	empi_to_date_oi_dic = {empi:pd.to_datetime(date) for empi, date in zip(single_rp_patients.EMPI.values, single_rp_patients.prdate_parsed.values)}

	print('\n')
	print('Processing patients who received radical prostatectomy...')
	df_labs_processed_rp = get_psa_prior_to_rp(df_labs, empi_to_date_oi_dic, empis_oi, export_dic = False, mode = 'rp')
	# get rp with valid bcr
	df_rp_bcr = df_labs_processed_rp.dropna(subset = ['bcr_post_op'])
	df_count_stats = pd.value_counts(df_rp_bcr.bcr_post_op.values, sort = True)
	print(df_count_stats)
	print('BCR rates', np.round(df_count_stats[1]/df_count_stats.values.sum(), 3))
	empis_to_bcr_info_dic = get_bcr_info(df_labs_processed_rp, empis_to_bcr_info_dic, mode = 'rp')

	# export empi id to bcr info dictionary 
	f = open("../data/empis_to_bcr_info_dic.pkl", "wb")
	pickle.dump(empis_to_bcr_info_dic,f)
	f.close()
	breakpoint()

	print('\n')
	print('The number of patients who received both radiation and RP : ', len(set(df_labs_processed_rad_final.index) & set(df_rp_bcr.index)))

	# exclude patients who received both therapies for now:
	patient_empis_oi = set(df_labs_processed_rad_final.index) | set(df_rp_bcr.index) - (set(df_labs_processed_rad_final.index) & set(df_rp_bcr.index))
	df_lstm_ready = pd.DataFrame([], index = np.arange(120), columns = patient_empis_oi)

	for empi, pre_op_psa, pre_op_psa_time in tqdm(zip(df_labs_processed_rad_final.index, 
													  df_labs_processed_rad_final.pre_op_psa_trajectory_10_year_window.values, 
													  df_labs_processed_rad_final.pre_op_psa_trajectory_dates_relative.values), total = len(df_labs_processed_rad_final), desc = 'Transforming PSA trajectories for radiation patients...'):
		if empi in df_lstm_ready.columns:
			pre_op_psa_time_wrt_zero = np.asarray(pre_op_psa_time) + 120
			for psa, time in zip(pre_op_psa, pre_op_psa_time_wrt_zero):
				df_lstm_ready.at[time, empi] = psa

	for empi, pre_op_psa, pre_op_psa_time in tqdm(zip(df_rp_bcr.index, 
													  df_rp_bcr.pre_op_psa_trajectory_10_year_window.values, 
													  df_rp_bcr.pre_op_psa_trajectory_dates_relative.values), total = len(df_rp_bcr), desc = 'Transforming PSA trajectories for RP patients...'):
		if empi in df_lstm_ready.columns:
			pre_op_psa_time_wrt_zero = np.asarray(pre_op_psa_time) + 120
			for psa, time in zip(pre_op_psa, pre_op_psa_time_wrt_zero):
				df_lstm_ready.at[time, empi] = psa

	# filter out patients without any prior psa
	df_num_points_per_empi = pd.DataFrame(0, index = df_lstm_ready.columns, columns = ['num_points'])
	for col in df_lstm_ready.columns:
		df_num_points_per_empi.at[col, 'num_points'] = len(df_lstm_ready[col].dropna())
	zero_psa_empis = df_num_points_per_empi.loc[df_num_points_per_empi.num_points == 0].index
	df_num_points_per_empi_non_zero = df_num_points_per_empi.loc[set(df_num_points_per_empi.index) - set(zero_psa_empis)]
	# df_num_points_per_empi_non_zero
	print('\n')
	print('Stats : patients with non-zero PSA points')
	print(df_num_points_per_empi_non_zero.num_points.describe())

	df_prior_psa_stats = pd.DataFrame(0, index = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], columns = ['remaining number of patients'])
	for thresh in df_prior_psa_stats.index:
		df_prior_psa_stats.at[thresh, 'remaining number of patients'] = len(df_num_points_per_empi_non_zero.loc[df_num_points_per_empi_non_zero.num_points >= thresh])


	df_prior_psa_stats = pd.DataFrame(0, index = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], columns = ['remaining number of patients'])
	for thresh in df_prior_psa_stats.index:
		df_prior_psa_stats.at[thresh, 'remaining number of patients'] = len(df_num_points_per_empi.loc[df_num_points_per_empi.num_points >= thresh])
	print('Remaining number of patients at each threshold : ')
	print(df_prior_psa_stats)
	breakpoint()
	return

if __name__ == "__main__":
	main()