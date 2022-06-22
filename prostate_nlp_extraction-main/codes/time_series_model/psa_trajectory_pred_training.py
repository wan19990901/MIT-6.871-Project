import pandas as pd
import numpy as np
import pickle
from tqdm import tqdm

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import torch
import torch.nn as nn
import torch.autograd as autograd
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader

import pytorch_lightning as pl
from pytorch_lightning.callbacks import ModelCheckpoint, EarlyStopping
from pytorch_lightning.loggers import TensorBoardLogger
from pytorch_lightning.metrics.functional import accuracy

import matplotlib.pyplot as plt
from matplotlib import rc
import multiprocessing as mp

from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import classification_report, confusion_matrix

OBS_TIME_WINDOW = 121 # months
EVENT_OBS_WINDOW = 5 # Years
N_EPOCHS = 25 #250
BATCH_SIZE = 64

# create dataset module
class SurfaceDataset(Dataset):
	
	def __init__(self, sequences):
		self.sequences = sequences
	def __len__(self):
		return len(self.sequences)
	def __getitem__(self, idx):
		sequence, label = self.sequences[idx]
		return dict(
			sequence = torch.Tensor(sequence.to_numpy()),
			label = torch.tensor(label).long()
		)
class SurfaceDataModule(pl.LightningDataModule):
	
	def __init__(self, train_sequences, test_sequences, batch_size):
		super().__init__()
		self.train_sequences = train_sequences
		self.test_sequences = test_sequences
		self.batch_size = batch_size
		
	def setup(self, stage = None):
		self.train_dataset = SurfaceDataset(self.train_sequences)
		self.test_dataset = SurfaceDataset(self.test_sequences)

	def train_dataloader(self):
		return DataLoader(
			self.train_dataset,
			batch_size = self.batch_size,
			shuffle = True,
			num_workers = mp.cpu_count()
		)
	def val_dataloader(self):
		return DataLoader(
			self.test_dataset,
			batch_size = self.batch_size,
			shuffle = False,
			num_workers = mp.cpu_count()
		)
	def test_dataloader(self):
		return DataLoader(
			self.test_dataset,
			batch_size = self.batch_size,
			shuffle = False,
			num_workers = mp.cpu_count()
		)

# create training module
class SequenceModel(nn.Module):
	
	def __init__(self, n_features, n_classes, n_hidden = 256, n_layers = 3):
		super().__init__()
		
#         self.n_hidden = n_hidden
		self.lstm = nn.LSTM(
			input_size = n_features,
			hidden_size = n_hidden,
			num_layers = n_layers,
			batch_first = True,
			dropout = 0.75
		)
		
		self.classifier = nn.Linear(n_hidden, n_classes)
	
	def forward(self, x):
		self.lstm.flatten_parameters()
		_, (hidden, _) = self.lstm(x) #output of the hidden unit from last layer
		out = hidden[-1]
		return self.classifier(out)

class SurfacePredictor(pl.LightningModule):
	
	def __init__(self, n_features : int, n_classes : int):
		super().__init__() #call super constructor
		self.model = SequenceModel(n_features, n_classes)
		self.criterion = nn.CrossEntropyLoss()
	
	def forward(self, x, labels = None):
		output = self.model(x)
		loss = 0
		if labels is not None:
			loss = self.criterion(output, labels)
		return loss, output
	
	def training_step(self, batch, batch_idx):
		sequences = batch['sequence']
		labels = batch['label']
		loss, outputs = self(sequences, labels)
		predictions = torch.argmax(outputs, dim = 1)
		step_accuracy = accuracy(predictions, labels)
		
		self.log('train_loss', loss, prog_bar = True, logger = True)
		self.log('train_accuracy', step_accuracy, prog_bar = True, logger = True)
		return {'loss' : loss, 'accuracy' : step_accuracy}
	
	def validation_step(self, batch, batch_idx):
		sequences = batch['sequence']
		labels = batch['label']
		loss, outputs = self(sequences, labels)
		predictions = torch.argmax(outputs, dim = 1)
		step_accuracy = accuracy(predictions, labels)
		
		self.log('val_loss', loss, prog_bar = True, logger = True)
		self.log('val_accuracy', step_accuracy, prog_bar = True, logger = True)
		return {'loss' : loss, 'accuracy' : step_accuracy}

	def test_step(self, batch, batch_idx):
		sequences = batch['sequence']
		labels = batch['label']
		loss, outputs = self(sequences, labels)
		predictions = torch.argmax(outputs, dim = 1)
		step_accuracy = accuracy(predictions, labels)
		
		self.log('test_loss', loss, prog_bar = True, logger = True)
		self.log('test_accuracy', step_accuracy, prog_bar = True, logger = True)
		return {'loss' : loss, 'accuracy' : step_accuracy}
	
	def configure_optimizers(self):
		return optim.Adam(self.parameters(), lr = 0.0001)


def train_model(X_data, empis_to_bcr_info_dic):
	sequences = []
	FEATURE_COLUMNS = ['prior_psa', 'overall_grade', 'rp_ind', 'rad_ind']
	bcr_counter = 0
	for series_id, group in X_data.groupby('empi'):
		sequence_features = group[FEATURE_COLUMNS]
		if len(sequence_features['prior_psa'].dropna()) != 0:
			# breakpoint()
			sequence_features_filled = sequence_features.fillna(0)
			bcr_ind = empis_to_bcr_info_dic[series_id][0]
			event_date_minus_op = empis_to_bcr_info_dic[series_id][1]
			if event_date_minus_op < EVENT_OBS_WINDOW*365.25:
				label = empis_to_bcr_info_dic[series_id][0]
				sequences.append((sequence_features_filled, label))
				if label == 1:
					bcr_counter += 1
		# else:
			# breakpoint()
		# breakpoint()
	# split train and test
	# breakpoint()
	train_sequences, test_sequences = train_test_split(sequences, test_size = 0.2)
	print('num train sequences : ', len(train_sequences))
	print('num test sequences : ', len(test_sequences))
	# breakpoint()
	data_module = SurfaceDataModule(train_sequences, test_sequences, BATCH_SIZE)
	model = SurfacePredictor(n_features = len(FEATURE_COLUMNS), n_classes = 2)

	checkpoint_callback = ModelCheckpoint(
		dirpath = 'checkpoints/', 
		filename = 'best-checkpoint-{epoch:02d}-{val_loss:.2f}',
		save_top_k = 1,
		verbose = True,
		monitor = 'val_loss',
		mode = 'min'
	)

	logger = TensorBoardLogger('lightning_logs', name = 'surface')
	trainer = pl.Trainer(
		# logger = logger,
		callbacks = [checkpoint_callback],
		max_epochs = N_EPOCHS,
		progress_bar_refresh_rate = 30
	)

	trainer.fit(model, data_module)
	trainer.test(model)
	# breakpoint()
	trained_model = SurfacePredictor.load_from_checkpoint(
		trainer.checkpoint_callback.best_model_path,
		n_features = len(FEATURE_COLUMNS),
		n_classes = 2
	)

	trained_model.freeze()

	test_dataset = SurfaceDataset(test_sequences)

	predictions = []
	labels = []

	for item in tqdm(test_dataset):
		sequence = item["sequence"]
		label = item['label']
		
		_, output = trained_model(sequence.unsqueeze(dim = 0))
		prediction = torch.argmax(output, dim = 1)
		predictions.append(prediction.item())
		labels.append(label.item())

	print(classification_report(labels, predictions))

	cm = confusion_matrix(labels, predictions)
	df_cm = pd.DataFrame(
		cm, index = [0,1], columns = [0,1]
	)
	return

def main():
	load_pre_processed_data = input('Load pre-processed pathology data? (Y/N) : ') == 'Y'
	# get relevant data
	df_lstm_ready = pd.read_csv('../data/df_lstm_ready.csv')
	df_lstm_ready.set_index(df_lstm_ready.columns[0], inplace = True)
	# all empi patients with bcr info
	with open('../data/empis_to_bcr_info_dic.pkl', 'rb') as handle:
		empis_to_bcr_info_dic = pickle.load(handle)	
	# all empis in the data
	if not load_pre_processed_data:
		# all patients with overall grade info available
		df_pathology_with_biopsy_overall_grade = pd.read_csv('../data/df_pathology_with_biopsy_oi_overall_grade.csv', index_col = 'EMPI')
		
		all_empis_oi = [int(col) for col in df_lstm_ready.columns]

		empis_rp_rad_in_biopsy_reports = set(df_pathology_with_biopsy_overall_grade.index) & set(all_empis_oi)
		# df_pathology_with_biopsy_overall_grade.set_index('EMPI', inplace = True)
		df_pathology_with_biopsy_overall_grade['rp_rad_indicator'] = 0
		# df_pathology_with_biopsy_overall_grade.loc[empis_rp_in_biopsy_reports, 'rp_rad_indicator'] = 1
		df_pathology_with_biopsy_overall_grade['rp_rad_date'] = None#[empis_to_bcr_info_dic[empi][2] for empi in df_pathology_with_biopsy_overall_grade.index]
		df_pathology_with_biopsy_overall_grade['rp_rad_date_minus_report_date_in_days'] = None
		df_pathology_with_biopsy_overall_grade['relative_month'] = None
		# df_pathology_with_biopsy_overall_grade['Report_Date_Time'] = pd.to_datetime(df_pathology_with_biopsy_overall_grade['Report_Date_Time'].values)
		for idx, empi in tqdm(enumerate(df_pathology_with_biopsy_overall_grade.index), total = len(df_pathology_with_biopsy_overall_grade), desc = 'Getting RP/rad date...'):
			if empi in empis_rp_rad_in_biopsy_reports:
				rp_rad_date = empis_to_bcr_info_dic[empi][2]
				df_pathology_with_biopsy_overall_grade.iat[idx, list(df_pathology_with_biopsy_overall_grade.columns).index('rp_rad_indicator')] = 1
				df_pathology_with_biopsy_overall_grade.iat[idx, list(df_pathology_with_biopsy_overall_grade.columns).index('rp_rad_date')] = rp_rad_date
				report_date = df_pathology_with_biopsy_overall_grade.iat[idx, list(df_pathology_with_biopsy_overall_grade.columns).index('Report_Date_Time')]
				report_dt = pd.to_datetime(df_pathology_with_biopsy_overall_grade.iat[idx, list(df_pathology_with_biopsy_overall_grade.columns).index('Report_Date_Time')])
				df_pathology_with_biopsy_overall_grade.iat[idx, list(df_pathology_with_biopsy_overall_grade.columns).index('rp_rad_date_minus_report_date_in_days')] = (rp_rad_date - report_dt).days
				# relative times are rounded: (e.g.) -1 * np.round((RP_date - tuple_oi[1]).days/30.44)
				df_pathology_with_biopsy_overall_grade.iat[idx, list(df_pathology_with_biopsy_overall_grade.columns).index('relative_month')] = -1 * np.round((rp_rad_date - report_dt).days/30.44) + OBS_TIME_WINDOW

		# make sure not to include any post-rad-rp grade info
		df_pathology_with_biopsy_overall_grade_valid = df_pathology_with_biopsy_overall_grade.loc[df_pathology_with_biopsy_overall_grade.rp_rad_date_minus_report_date_in_days >= 0]
		df_pathology_with_biopsy_overall_grade_valid = df_pathology_with_biopsy_overall_grade_valid.loc[df_pathology_with_biopsy_overall_grade_valid.relative_month >= 0]
		empi_to_overall_grade_time_point_dic = {}
		for empi, overall_grade, time_point in zip(df_pathology_with_biopsy_overall_grade_valid.index, 
												   df_pathology_with_biopsy_overall_grade_valid.overall_grade_merged.values, 
												   df_pathology_with_biopsy_overall_grade_valid.relative_month.values):
			if empi not in empi_to_overall_grade_time_point_dic.keys():
				empi_to_overall_grade_time_point_dic[empi] = set([(overall_grade, time_point)])
			else:
				empi_to_overall_grade_time_point_dic[empi] = empi_to_overall_grade_time_point_dic[empi] | set([(overall_grade, time_point)])
		# f = open("../data/empi_to_overall_grade_time_point_dic.pkl", "wb") # cancer_type_to_eligible_feats_dic_Jan_27th_2021, cancer_type_to_eligible_feats_dic_Sep_25th
		# pickle.dump(empi_to_overall_grade_time_point_dic,f)
		# f.close()
	# else:
		# with open('../data/empi_to_overall_grade_time_point_dic.pkl', 'rb') as handle:
		# 	empi_to_overall_grade_time_point_dic = pickle.load(handle)	
		# breakpoint()
		# build trainable df
		# NUM_MONTHS = 120
		intital_flag = 1
		for empi in tqdm(df_lstm_ready.columns, desc = 'Getting feature df ready...'):
			df_features_local = pd.DataFrame([], index = np.arange(OBS_TIME_WINDOW),columns = ['empi', 'month', 'prior_psa', 'overall_grade', 'rp_ind', 'rad_ind'])
			df_oi = df_lstm_ready[empi].dropna()
			df_features_local['empi'] = empi
			df_features_local['month'] = np.arange(OBS_TIME_WINDOW)
			for idx, time_point, psa in zip(np.arange(len(df_oi)), df_oi.index, df_oi.values):
				if idx < len(df_oi) - 1:
					next_time_point = list(df_oi.index)[idx + 1]
				else: # at the last entry
					next_time_point = OBS_TIME_WINDOW
				for time in range(int(time_point), next_time_point):
					df_features_local.at[time, 'prior_psa'] = psa

			empi_int = int(empi)
			empi_rp_rad_info = empis_to_bcr_info_dic[empi_int][3]
			if empi_int in empi_to_overall_grade_time_point_dic.keys():
				for entry_info in list(empi_to_overall_grade_time_point_dic[empi_int]):
					overall_grade = entry_info[0]; time_point = entry_info[1];
					for time in range(int(time_point), OBS_TIME_WINDOW):
						df_features_local.at[time, 'overall_grade'] = overall_grade
						if empi_rp_rad_info == 'rp':
							df_features_local.at[time, 'rp_ind'] = 1
						elif empi_rp_rad_info == 'rad':
							df_features_local.at[time, 'rad_ind'] = 1
			   #  overall_grade = df_pathology_with_biopsy_overall_grade_valid[empi, 'overall_grade_merged']; time_point = ;
			   #  df_features_local['overall_grade'] =  
			if intital_flag:
				df_features_total = df_features_local
			else:
				df_features_total = pd.concat([df_features_total, df_features_local])

			if intital_flag == 1:
				intital_flag = 0
		df_features_total.to_csv('../data/df_features_total_lstm.csv')
	else:
		df_features_total = pd.read_csv('../data/df_features_total_lstm.csv')
	
	 # month
	# for val in 
	train_model(df_features_total, empis_to_bcr_info_dic)
	return

if __name__ == "__main__":
	main()