# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 14:23:54 2015

@author: S.Brad
"""


###############################################################################
print("SECTION: INITIAL SETUP","\n")
###############################################################################

input_directory = 'F:\\Dropbox\\Research\\Hedge_Fund_Databases\\Data'
output_directory = 'F:\\Import_Data\\Data\\Eurekahedge'
function_directory = 'F:\\Dropbox\\Research_Methods\\R'


###############################################################################
print("SECTION: LIBRARIES","\n")
###############################################################################

from pandas import read_csv


###############################################################################
print("SECTION: SETUP","\n")
###############################################################################

identifier = 'Fund_ID'

final_folder_path = '\\'.join([output_directory,'Final_Expand3'])
melt_folder_path = '\\'.join([output_directory,'NAV_AUM_melt'])
revision_folder_path = '\\'.join([output_directory,'Revision'])

Dead_cols_keep = ['pull',identifier,'Dead_Date','date','min_date','max_date','bad_min','bad_max']

Dead_dates = read_csv(''.join([final_folder_path,'\\','EurekahedgeHF_NAV_AUM_Ret','.csv']))
#print(Dead_dates)
#print(Dead_dates['age'])
