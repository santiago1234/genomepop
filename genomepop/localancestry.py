'''
The localancestry submodule of the genomepop module is designed to process
the output of local ancestry inference, providing methods for downstream
analysis of the data.
The primary focus of this submodule is to facilitate the
calculation of ancestry tract length distributions,
allowing for easy and accurate analysis of genetic variation.
The module provides a user-friendly interface for working with local ancestry
data and can be easily integrated into larger pipelines for
genetic analysis.
'''

import os

import pandas as pd
import numpy as np

test_msp = '/Users/santiagomedina/caapa-genomes/analysis/2301-LAI-TEST/test_script/query_results.msp'
# The non-sample columns of the msp file
_main_cols = ['#chm', 'spos', 'epos', 'sgpos', 'egpos', 'n snps']


def load_pop_codes(msp_file):
    """
    Reads the population codes from
    the msp_file.
    Args:
        msp_file: msp file generate by gnomix
    The population codes are included
    in the first line of the msp file.
    """
    mf = open(msp_file, 'r')
    header_codes = mf.readline()
    mf.close()
    pop_codes = (header_codes.replace('#Subpopulation order/codes: ',
                                      '').strip().split('\t'))
    return {int(b): a for (a, b) in [x.split('=') for x in pop_codes]}


def load_msp_file(msp_file):
    """
    Loads the msp_file as a pandas.DataFrame
    Note: The output file with contain the label codes.
    """
    msp = pd.read_table(msp_file, skiprows=1)
    p_codes = load_pop_codes(msp_file)
    samples_cols = [x for x in msp.columns if x not in _main_cols]
    msp[samples_cols] = msp[samples_cols].applymap(lambda x: p_codes[x])
    return msp


def get_individual(msp, individual):
    """
    Retrieves the data for a particular
    individual
    Args:
        msp: pd.DataFrame
        individual: str, sample name
    """
    ind_haps = [individual + '.' + str(x) for x in [0, 1]]
    if ind_haps[0] not in msp.columns:
        raise ValueError('Supplied individual not found')
    return msp.loc[:, _main_cols + ind_haps]


def tidy_individual(ind_data):
    """
    Args:
        ind_data: pd.DataFrame
    """
    # make long format data
    ind_data_l = ind_data.melt(id_vars=_main_cols, value_name='Ancestry')
    # happlotype make 0s to A and 1s to B.
    # That way is more intuitive
    haplo_mapper = {'0': 'A', '1': 'B'}
    # the last digit is the haplotype
    ind_data_l['Haplotype'] = ind_data_l.variable.str[-1].map(haplo_mapper)
    ind_data_l['Individual'] = ind_data_l.variable.str[:-2]
    # rename columns to match Tracts readme columns
    col_renamer = {'#chm': 'chrn'}
    return ind_data_l.rename(col_renamer, axis=1)


def is_ancestry_continuous(r_1, r_2, max_gap=1e6):
    """
    Params:
        r_1: pd.DataFrame, a single row pandas frame.
            this row contains the coulummn Ancestry
        r_2: same as r_1
        max_gap: maximum gap, if the distance between
            the windows is more than max_gap then they
            wont be merged
    Returns:
        lgl: True if r_1 and r_2 represent the same ancestry
    Note:
        I assume that r_2 follows from r_1. I need
        to write conditions to make sure that r_1 and r_2
        if merged together represent continuous chromosome
        fragments.
        Also if there is a gap of more than max_gap, the
        fragments are not considered to be continous
    """
    gap = r_2.spos.values[0] - r_1.epos.values[0]
    if gap > max_gap:
        return False

    res = r_1.Ancestry.values[0] == r_2.Ancestry.values[0]
    return res


def merge_ancestries(r_1, r_2):
    """
    Replace the end position (physical and recombination)
    in r_1 with the values in r_2.
    Args:
        r_1: row 1
        r_2: row 2
    Returns: merged (1 row) pandas.DataFrame
    """
    merged = r_1.reset_index(drop=True)
    r_2 = r_2.reset_index(drop=True)
    merged.loc[0, 'epos'] = r_2.loc[0, 'epos']
    merged.loc[0, 'egpos'] = r_2.loc[0, 'egpos']

    return merged


def collapse_windows_to_tracts(df_chm, max_gap):
    """
    If n continous windows are from the same ancestry,
    collpase to a single continous windows.
    Args:
        df_chm: pd.DataFrame
            Table with genomo coordinates and Ancestry assigments.
    Note:
        It is assumed that the input data represents one chromosome
        and only one Haplotype either maternal or paternal
    """
    # make sure df_chm is sorted by physical position
    df_chm = df_chm.sort_values(['spos'])

    current_ancstr_block = df_chm.iloc[[0]]
    tracts = list()

    for index in range(1, df_chm.shape[0]):
        current_row = df_chm.iloc[[index]]

        if is_ancestry_continuous(current_ancstr_block, current_row, max_gap):
            current_ancstr_block = merge_ancestries(current_ancstr_block,
                                                    current_row)

        else:
            tracts.append(current_ancstr_block)
            current_ancstr_block = current_row

    # add the last ancestry, the for loop does not include the last row

    tracts.append(current_ancstr_block)

    tracts = pd.concat(tracts).reset_index(drop=True)

    # compute tracts length
    tracts['len_bp'] = tracts.epos - tracts.spos
    tracts['len_cm'] = tracts.egpos - tracts.sgpos

    return tracts.reset_index(drop=True)
