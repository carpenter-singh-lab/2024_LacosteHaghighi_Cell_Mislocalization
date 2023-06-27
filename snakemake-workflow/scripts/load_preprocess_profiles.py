import numpy as np
import scipy.spatial
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from random import sample,choices
from scipy.stats import pearsonr
import os
from sklearn import preprocessing

import sys
sys.path.insert(0, '/home/ubuntu/workspace_SingleCell/SingleCell_Morphological_Analysis/') 

# from singlecell.read import read_single_cell_sql
from singlecell.preprocess import handle_nans, extract_cpfeature_names
from singlecell.process.normalize_funcs import *
# from singlecell.preprocess.filter_out_untransfected import extract_singlecell_transfection_labels

#####################################################################
 loading configuration parameters from yaml file
with open('config.yaml') as file:
    config_params = yaml.safe_load(file)

def read_merge_preprocess_meanProfiles(annot_df_batches,read_pop_params):
    
#     annot_df_batches=annot_df.copy()
#     annot_df_batches=annot_df_sets.copy()
    feature_params=read_pop_params['feature_scaling_params_dict']

    rootDir=read_pop_params['dirs_params_dict']['rootDir']
    profiles_folder_in_workspace=read_pop_params['dirs_params_dict']['profiles_folder_in_workspace']

    batchNames=annot_df_batches['batch'].unique().tolist()

    df_batch_plates_annot_ls=[]
#     dfTransSummary_batch_ls=[]
    cp_features_analysis_batch_ls=[]

    for bi in range(len(batchNames)):
        batchName=batchNames[bi]
        annot_df=annot_df_batches[annot_df_batches['batch']==batchName].reset_index(drop=True)
        metaDataPlates=annot_df['Metadata_Plate'].unique()

    #         profiles_address=rootDir+'/backend/'+profType+'PerWells'+ndd+'/'    # for set1 and set2

        listOfPlates=annot_df['Metadata_Plate'].unique().tolist()

        df_batch_plates_ls=[]
        for p in listOfPlates:
    #         profiles_address_old=rootDir+'/'+profiles_folder_in_workspace+'/'+batchName+'/'+p+'/population_profiles_'+p+'.csv.gz'
            profiles_address=rootDir+'/'+profiles_folder_in_workspace+'/'+batchName+'/'+p+'/'+p+'.csv.gz'

    #         os.rename(profiles_address_old, profiles_address)

            print(profiles_address)
            if os.path.exists(profiles_address):
                perplate_df=pd.read_csv(profiles_address)

                if feature_params['zscored_profiles'][0]:
                    if 'zscored' in perplate_df.columns:
                        perplate_df=perplate_df[perplate_df['zscored']==feature_params['zscored_profiles'][1]].\
                        reset_index(drop=True)
                    else:
                        raise ValueError('This batch doent have zscored column! '
                                         'Probably no untransfected cells were assigned! check and carefully change'
                                         'zscored_profiles enabled flag to False!')
                else:
                    if 'zscored' in perplate_df.columns:
                        perplate_df=perplate_df[perplate_df['zscored'].isnull()].reset_index(drop=True);

                perplate_df=perplate_df[perplate_df['normalization']==feature_params['sc_per_plate_scaling']].\
                reset_index(drop=True)            


            df_batch_plates_ls.append(perplate_df[perplate_df['n_transf']>0].reset_index(drop=True))

        df_batch=pd.concat(df_batch_plates_ls,ignore_index=True)    

        cp_features, cp_features_analysis = extract_cpfeature_names.extract_cpfeature_names(df_batch);
    #   df_batch_plates, cp_features_analysis_batch = handle_nans.handle_nans(df,cp_features_analysis,\
    #                                                    thrsh_null_ratio=0.05,fill_na_method='drop-rows');  

        df_batch_plates, cp_features_analysis_batch = handle_nans.handle_nans(df_batch,cp_features_analysis,\
                                                       thrsh_null_ratio=0);          

    #         print('df_batch_plates1',df_batch_plates['Category'].unique())
    #         print("cp_features_analysis_batch",cp_features_analysis_batch)
        cp_features_analysis_batch_ls.append(cp_features_analysis_batch)    

        df_batch_plates_annot_ls.append(df_batch_plates)


    df_annot_batches=pd.concat(df_batch_plates_annot_ls,ignore_index=True)
    # dfTransSummary_batches=pd.concat(dfTransSummary_batch_ls,ignore_index=True)

    feats_intersection_batches=cp_features_analysis_batch_ls[0]
    for bi in range(1,len(cp_features_analysis_batch_ls)):
        feats_intersection_batches=list(set(feats_intersection_batches) &\
                                               set(cp_features_analysis_batch_ls[bi]))
        print(len(feats_intersection_batches),len(cp_features_analysis_batch_ls[bi]))
        


        
    if feature_params['post_scale_all_profiles'][0]:
#             dataScaled=scaler.fit_transform(df2.loc[:,cpFeatures4scale])
#             df_scaled[cpFeatures4scale]=dataScaled
        df_scaled_annot_batches = standardize_df_columns(df_annot_batches,feats_intersection_batches,\
                                                         feature_params['post_scale_all_profiles'][1])
    else:
        df_scaled_annot_batches=df_annot_batches.copy()    
    
    
    df_scaled_annot_batches, feats_intersection_batches = \
    handle_nans.handle_nans(df_scaled_annot_batches,feats_intersection_batches, thrsh_null_ratio=0)    

    ######################################################
#     cpFeats_P=[p for p in feats_intersection_batches if 'Protein' in p]
#     cpFeats_NP=[p for p in feats_intersection_batches if ('Protein' not in p) and ('DsRed' not in p)]
    
    prot_ch_suffix = read_pop_params['protein_channel_suffix']

    cpFeats_P=[p for p in feats_intersection_batches if prot_ch_suffix in p]
    cpFeats_NP=[p for p in feats_intersection_batches if ('prot_ch_suffix' not in p)]
 

    return df_scaled_annot_batches, feats_intersection_batches,cpFeats_P,cpFeats_NP



# main script
def main():

    rootDir = config_params['rootDir']
    batch_12 = config_params['batch_12']
    channels_used = config_params['channels_used']
    intFeatures = config_params['intFeatures']
    feature_scaling_params_dict = config_params['feature_scaling_params_dict']
    dirs_params_dict = config_params['dirs_params_dict']

    read_pop_params = {
        'dirs_params_dict': dirs_params_dict,
        'feature_scaling_params_dict': feature_scaling_params_dict,
        'protein_channel_suffix': channels_used[0]
    }

    annot_df = pd.read_csv(os.path.join(rootDir, "image_level_annot_df.csv"))

    df_scaled_annot, cpFeats_A, cpFeats_P, cpFeats_NP = read_merge_preprocess_meanProfiles(annot_df, read_pop_params)

    dfTransSummary = df_scaled_annot[annot_df.columns.tolist() + ['n_transf', 'n_untransf', 'transf_Ratio']]
    df_scaled_annot = df_scaled_annot[df_scaled_annot['Metadata_transfection'] == False]
    df_scaled_annot = df_scaled_annot[df_scaled_annot['n_transf'] > 5].reset_index(drop=True)

    # Now you can continue processing df_scaled_annot, cpFeats_A, cpFeats_P, cpFeats_NP according to your needs

if __name__ == "__main__":
    main()