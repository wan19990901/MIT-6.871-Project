import numpy as np
import pandas as pd

import re
from tqdm import tqdm

import datefinder


def post_process_psa_result(df_psa_result):
    # creat psa dictionary based on df
    empi_to_psa_trajectory_date_dic = {}
    for empi, psa_val, date, comment in zip(df_psa_result.index, df_psa_result.psa_val.values,
                                            df_psa_result.date.values, df_psa_result.comment.values):
        empi_true = empi.split('_')[0]
        if empi_true in empi_to_psa_trajectory_date_dic.keys():
            if not pd.isnull(comment):
                empi_to_psa_trajectory_date_dic[empi_true] = empi_to_psa_trajectory_date_dic[empi_true] | {
                    (psa_val, date, '<')}
            else:
                empi_to_psa_trajectory_date_dic[empi_true] = empi_to_psa_trajectory_date_dic[empi_true] | {
                    (psa_val, date, '')}
        else:
            if not pd.isnull(comment):
                empi_to_psa_trajectory_date_dic[empi_true] = {(psa_val, date, '<')}
            else:
                empi_to_psa_trajectory_date_dic[empi_true] = {(psa_val, date, '')}

    # create df based on dictionary
    empi_initial = 0
    for empi, psa_vals in empi_to_psa_trajectory_date_dic.items():
        idx_counter = 0;
        df_psa_result_per_empi = pd.DataFrame([], index=np.arange(len(psa_vals)),
                                              columns=['empi', 'psa_val', 'date', 'comment'])
        for psa_tuple in psa_vals:
            df_psa_result_per_empi.iat[idx_counter, 0] = empi
            df_psa_result_per_empi.iat[idx_counter, 1] = psa_tuple[0]
            df_psa_result_per_empi.iat[idx_counter, 2] = psa_tuple[1]
            df_psa_result_per_empi.iat[idx_counter, 3] = psa_tuple[2]
            idx_counter += 1
        df_psa_result_per_empi.sort_values(by='date', inplace=True)
        if empi_initial == 0:
            df_psa_result_total = df_psa_result_per_empi
            empi_initial = 1
        else:
            df_psa_result_total = pd.concat([df_psa_result_total, df_psa_result_per_empi])

    df_psa_result_total.index = np.arange(len(df_psa_result_total))
    return df_psa_result_total


def fast_remove(str_oi):
    chars_to_remove = '():!@#$?psa'
    word_to_return = ''
    for char in str_oi:
        if not char in chars_to_remove:
            word_to_return = word_to_return + char
    return word_to_return


def remove_contents_in_braces(str_oi):
    word_to_return = ''
    for char in str_oi:
        if char == '[':
            break
        else:
            word_to_return = word_to_return + char
    return word_to_return


def get_dates_wo_noises(str_oi):
    word_to_return = []
    str_oi_split = str_oi.split('/')
    for str_ in str_oi_split:
        word_to_return.append(re.sub('[!@#$?psa]', '', str_))
    return '/'.join(word_to_return)


def get_psa_vals(df_prg_oi):
    df_result = pd.DataFrame([], columns=['psa_val', 'date', 'comment']);
    empi_psa_list = set();
    empi_idx_counter = 0
    empi_prev = ''
    for empi, report in tqdm(zip(df_prg_oi.index, df_prg_oi.Report_Text.values), total=len(df_prg_oi),
                             desc='Getting PSA vals...'):
        #     empi_idx_counter =
        if empi != empi_prev:
            empi_idx_counter = 0
        report_lower_case = report.lower().replace('*', '')
        if 'psa' in report_lower_case:
            report_split = report_lower_case.split(' ')
            for idx, word in enumerate(report_split):
                ineqaul_ind = 0
                if 'psa' in word:
                    potential_date_vals = report_split[idx - 10, idx + 10]
                    potential_dates_found = datefinder.find_dates(potential_date_vals)
                    for potential_date in potential_dates_found:
                        print(potential_date)
                    # checking for "XX (MM/DD/YYYY)", where XX is a PSA val and MM/DD/YYYY is a measurement date
                    if potential_date_val.count('/') == 2:
                        # record psa value
                        potential_psa_val = report_split[idx + 1]
                        if '<' in potential_psa_val:
                            # flag for < sign : note that some PSA values are recorded in inequality (e.g. < XX)
                            potential_psa_val_num = potential_psa_val.replace('<', '')
                            try:
                                df_result.at[str(empi) + '_' + str(empi_idx_counter), 'psa_val'] = float(
                                    potential_psa_val_num)
                                ineqaul_ind = 1
                            except:
                                # empi_idx_counter += 1
                                continue
                        else:
                            try:
                                df_result.at[str(empi) + '_' + str(empi_idx_counter), 'psa_val'] = float(
                                    potential_psa_val)
                            except:
                                # empi_idx_counter += 1
                                continue
                        # obtain psa measurement date after some post-processing
                        potential_date_val = fast_remove(potential_date_val)
                        if '[' in potential_date_val and ']' in potential_date_val:
                            potential_date_val = remove_contents_in_braces(potential_date_val)
                        try:
                            potential_date_val_final = pd.to_datetime(potential_date_val)
                            df_result.at[str(empi) + '_' + str(empi_idx_counter), 'date'] = potential_date_val_final
                            if ineqaul_ind:
                                df_result.at[str(empi) + '_' + str(empi_idx_counter), 'comment'] = '<'
                            empi_idx_counter += 1
                        except:
                            empi_idx_counter += 1
                            continue
        empi_prev = empi
    return df_result


def main():
    # load relevant files
    print('Loading relveant files...')
    df_prg_total = pd.read_csv('../data/progress_notes/df_prg_total.csv')
    df_prg_total.set_index('EMPI', inplace=True)

    df_pathology_biopsy_final = pd.read_csv('../data/df_pathology_biopsy_final.csv')
    df_pathology_biopsy_final.set_index('EMPI', inplace=True)
    df_labs_processed = pd.read_csv('../data/df_labs_processed_biopsy_based.csv')
    df_labs_processed.rename(columns={'Unnamed: 0': 'EMPI'}, inplace=True)
    df_labs_processed.set_index('EMPI', inplace=True)
    print('Loading complete!')
    print('\n')

    # get quick stats
    print('Number of unique patients in Progress reports : ', len(set(df_prg_total.index)))
    print('Number of unique patients in Biopsy reports : ', len(set(df_pathology_biopsy_final.index)))
    print('Number of unique patients in PSA test reports : ', len(set(df_labs_processed.index)))
    empi_in_pat_not_in_lab = set(df_pathology_biopsy_final.index) - set(df_labs_processed.index)
    empi_in_pat_not_in_lab_in_prg = empi_in_pat_not_in_lab & set(df_prg_total.index)
    print('Number of unique patients (not in PSA reports && in Biopsy reports && in Progress notes) : ',
          len(empi_in_pat_not_in_lab_in_prg))
    empi_in_psa_biopsy = set(df_prg_total.index) & set(df_pathology_biopsy_final.index)
    empi_in_psa_biopsy_progress = set(df_prg_total.index) & set(df_pathology_biopsy_final.index) & set(
        df_labs_processed.index)
    print('Number of unique patients (in Biopsy reports && in Progress notes) : ', len(empi_in_psa_biopsy))
    print('Number of unique patients (in PSA reports && in Biopsy reports && in Progress notes) : ',
          len(empi_in_psa_biopsy_progress))

    # choose cohort of interest
    df_prg_oi = df_prg_total.loc[empi_in_pat_not_in_lab_in_prg]
    df_result = get_psa_vals(df_prg_oi)
    breakpoint()
    return


if __name__ == "__main__":
    main()

