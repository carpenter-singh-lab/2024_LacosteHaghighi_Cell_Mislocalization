import os
import pandas as pd
import sys
import time
import yaml
from functools import reduce
import gc
import numpy as np
import pickle
import sys
sys.path.insert(0, '/home/ubuntu/workspace_SingleCell/SingleCell_Morphological_Analysis/') 

from singlecell.read import read_single_cell_sql
from singlecell.preprocess import handle_nans, extract_cpfeature_names
from singlecell.process.normalize_funcs import *
from singlecell.preprocess.filter_out_untransfected import extract_singlecell_transfection_labels
from singlecell.preprocess.filter_out_edge_single_cells import edgeCellFilter
from singlecell.save.save_pandas_dfs import saveDF_to_CSV_GZ_no_timestamp


# Import or define all the necessary functions you've used such as read_single_cell_sql, 
# handle_nans, edgeCellFilter, standardize_df_columns, extract_singlecell_transfection_labels,
# fill_enrichment_profile, saveDF_to_CSV_GZ_no_timestamp, generate_zscored_to_untransfected_profiles, etc.
# You can import them if they are defined in other files.
# For example:
# from utils import read_single_cell_sql, handle_nans, edgeCellFilter

def main(annot_df, root_dir, config_file):
    with open(config_file, 'r') as f:
        all_params = yaml.safe_load(f)

    listOfBatchPlates = annot_df.Metadata_batch_Plate.unique().tolist()

    for bp in listOfBatchPlates:
        b = bp.split('-')[0]
        all_params['transfection_params_dict']['precomputed_params'] = precomputed_params_batches[b]
        generate_population_profiles(bp, annot_df, root_dir, all_params)



def generate_population_profiles(batchPlateName,annot_df,rootPath, params):
    """ 
    This function reads single cell sql files and save single cells for each well in "singlecell_profiles" folder
  
    Inputs: 
    ++ batchPlateName (str) batch and plate info in batch_plate form
    ++ annot_df   (pandas df) metadata file 
    ++ rootPath (str) project root path for ex, '/home/ubuntu/bucket/projects/2017_10_19_Profiling_rare_ORFs/workspace' 
    ++ perPlateScaleFlag (dtype: int): a 0/1 flag indicating if we want to scale single cell\
    features per plate before saving or not
    ++ applyTransFilt (int) a 0/1 flag indicating if we want to save all the cells or just the transfected cells

    """    
    transfection_params_dict=params['transfection_params_dict']
    enrichement_profiles_params=params['enrichement_profiles_params']
    feature_scaler=params['feature_scaling_params_dict']['feature_scaler']
    
    
    batchName=annot_df[annot_df['Metadata_batch_Plate']==batchPlateName]['batch'].tolist()[0]
    plateName=annot_df[annot_df['Metadata_batch_Plate']==batchPlateName]['Metadata_Plate'].tolist()[0]   

    annot_df_plate=annot_df[annot_df['Metadata_batch_Plate']==batchPlateName].reset_index(drop=True)

    print("Processing ",batchName,plateName)

    ## Set path to save pupulation profiles
    save_population_profiles_path=rootPath+'/population_profiles/'+batchName+'/'+plateName
    save_sc_profiles_path=rootPath+'/singlecell_profiles/'+batchName+'/'+plateName
    save_enrichment_profiles_path=rootPath+'/enrichment_profiles/'+batchName+'/'+plateName
    

    ## path for sqlite files outputs of cell profiler and cytominer packages
    fileName=rootPath+"/backend/"+batchName+"/"+plateName+"/"+plateName+".sqlite"

    if not os.path.exists(fileName):
        print(fileName+' do not exist!')
        return
    
    os.system("mkdir -p "+save_population_profiles_path)    
    os.system("mkdir -p "+save_sc_profiles_path)        
    os.system("mkdir -p "+save_enrichment_profiles_path)   
    
    
    ## read single cell sql files
    start_time = time.time()
    sc_df = read_single_cell_sql.readSingleCellData_sqlalch(fileName,["cells", "cytoplasm", "nuclei"])

    ## filter features with low std and nans for than a threshold
    cp_features, cp_features_analysis_0 =  extract_cpfeature_names.extract_cpfeature_names(sc_df);
    sc_df, cp_features_analysis = handle_nans.handle_nans(sc_df,cp_features_analysis_0,fill_na_method='drop-rows');   
    print(sc_df[cp_features_analysis].shape)
    
    ## filter cells on the edge
    sc_df,_=edgeCellFilter(sc_df)
    print(sc_df[cp_features_analysis].shape)

    ## scale features per plate for analysis
    standardized_per_plate_sc = standardize_df_columns(sc_df,cp_features_analysis,feature_scaler)


    ## apply transfection filtering of single cell level of data based on the parameters
    # transfection_params_dict={'Method':'single_intensity_feature_thrsh',\
    #                               'intensity_feature_to_use':'Cells_Intensity_MeanIntensity_DsRed',\
    #                               'thresholding_method': 'precomputed_batch_specific_thrsh',\
    #                               'precomputed_params': [precomputed_bottom_thresh , precomputed_top_thresh , data_norm_used]} 

    if transfection_params_dict:
#         if "precomputed_params" in transfection_params_dict.keys() and transfection_params_dict["precomputed_params"]=='raw':
        
        transfection_labels = extract_singlecell_transfection_labels(sc_df,transfection_params_dict)

        sc_df['transfection_status']=transfection_labels
        standardized_per_plate_sc['transfection_status']=transfection_labels

    else:
        sc_df['transfection_status']=1
        standardized_per_plate_sc['transfection_status']=1

        
    if params['save_single_cells']:
        transfected_sc_df=sc_df[sc_df['transfection_status']==1].reset_index(drop=True)
        standardized_transfected_sc_df = standardize_df_columns(transfected_sc_df,cp_features_analysis,feature_scaler)

        wells=standardized_transfected_sc_df['Metadata_Well'].unique().tolist();
        
        for w in wells:
            fileNameToSave=save_sc_profiles_path+'/'+plateName+'_'+w;
#             standardized_transfected_sc_df[standardized_transfected_sc_df['Metadata_Well']==w].\
#             reset_index(drop=True).to_pickle(fileNameToSave);
            saveDF_to_CSV_GZ_no_timestamp(standardized_transfected_sc_df[standardized_transfected_sc_df\
                                                   ['Metadata_Well']==w].reset_index(drop=True),fileNameToSave)

        
        
        
    # create enrichement profiles
    if enrichement_profiles_params:
    #     clustering_model_file = './utils/kmeans_20clusters_fitted_on_population_profiles.sav'
#         clustering_model_file = './utils/kmeans_20clusters_fitted_on_sc_profiles_3.sav'
        clustering_model_file = enrichement_profiles_params['clustering_model_file']
        clustering_model_dict = pickle.load(open(clustering_model_file, 'rb'))

        transfected_sc_df=sc_df\
        [sc_df['transfection_status']==1].reset_index(drop=True)    

    #     transfected_standardized_per_plate_sc=standardized_per_plate_sc\
    #     [standardized_per_plate_sc['transfection_status']==1].reset_index(drop=True)
        enrichment_profile_df=transfected_sc_df.groupby('Metadata_Well').size().\
        reset_index().rename(columns={0:'count'})

        enrichment_profile_df['Metadata_Plate']=plateName

        for fi in ['p','np']:
            n_clusts=clustering_model_dict['n_clusts']
            cpFeats=clustering_model_dict['model_'+fi+'_feats']
            model_cls=clustering_model_dict['model_'+fi]
            scaler_pre_cls=clustering_model_dict['scaler_'+fi]
            enrichment_profile_df_tmp=transfected_sc_df.groupby('Metadata_Well').\
            apply(lambda x: fill_enrichment_profile(np.unique(model_cls.predict(scaler_pre_cls.transform(x[cpFeats])),\
                                                              return_counts=True),n_clusts)).reset_index()
            enrichment_profile_df=pd.concat([enrichment_profile_df,\
                                             enrichment_profile_df_tmp.apply(lambda x:x[0], axis=1, result_type='expand').\
                      add_prefix(fi+'_')],axis=1)

        fileNameToSave=save_enrichment_profiles_path+'/enrichment_profiles_'+plateName;
        saveDF_to_CSV_GZ_no_timestamp(enrichment_profile_df,fileNameToSave)
#         enrichment_profile_df.to_pickle(fileNameToSave);
    

    ## count per-well transfection summaries
    transfection_summary =\
    sc_df[['Metadata_Plate','Metadata_Well','transfection_status']].\
    groupby(['Metadata_Plate','Metadata_Well']).\
    apply(lambda r: pd.Series({ 
    "n_transf": sum(r["transfection_status"]==1), 
    "n_untransf": sum(r["transfection_status"]!=1),
    "transf_Ratio":sum(r["transfection_status"]==1)/r.shape[0]})).reset_index()


    save_population_profiles_zscored_by_same_well_untransfected=False
    ## save population profiles z scored by same well untransfected cells
    if save_population_profiles_zscored_by_same_well_untransfected:
        print('sfds')
        standardized_per_plate_sc.groupby('Metadata_Well').\
        apply(lambda g: zscore_by_untrasfected(g))


    ## generate population level profiles    
    population_profiles=sc_df.groupby(['Metadata_Plate','Metadata_Well','transfection_status'])\
    .mean().reset_index();
    print(sc_df[cp_features_analysis].shape,population_profiles.shape)
    population_profiles_standardized=standardized_per_plate_sc.\
    groupby(['Metadata_Plate','Metadata_Well','transfection_status']).\
    mean().reset_index();
    
    print(sc_df[cp_features_analysis].shape)
    print(population_profiles[cp_features_analysis].shape)
    
    ## augment profiles with metadata and transfection summary
    population_prof_aug = reduce(lambda  left,right: \
                                pd.merge(left,right,on=['Metadata_Plate','Metadata_Well'],
                                how='outer'), [population_profiles,transfection_summary,annot_df_plate])
    population_prof_standardized_aug = reduce(lambda  left,right: \
                                pd.merge(left,right,on=['Metadata_Plate','Metadata_Well'],
                                how='outer'), [population_profiles_standardized,transfection_summary,annot_df_plate])


    ## save raw averaged profiles together with perplate scaled and then per well averaged profiles
    population_prof_aug['normalization']='raw'
    population_prof_standardized_aug['normalization']='sc_scaled_per_plate'

    

    profile_list_2save=[population_prof_aug,population_prof_standardized_aug]
    
    ## generate zscored population profiles to untransfected population profiles

    if transfection_params_dict:
        population_prof_aug_zscored = generate_zscored_to_untransfected_profiles(population_prof_aug,\
                                                                                 cp_features_analysis,feature_scaler)
        population_prof_standardized_aug_zscored = generate_zscored_to_untransfected_profiles(\
                                  population_prof_standardized_aug,cp_features_analysis,feature_scaler)    

        profile_list_2save+=[population_prof_aug_zscored, population_prof_standardized_aug_zscored]
        
        
    population_profiles_2save=pd.concat(profile_list_2save,ignore_index=True)


#     fileNameToSave=save_population_profiles_path+'/population_profiles_'+plateName;
    fileNameToSave=save_population_profiles_path+'/'+plateName;
    saveDF_to_CSV_GZ_no_timestamp(population_profiles_2save,fileNameToSave)
#     population_profiles_2save.to_pickle(fileNameToSave);

    gc.collect()

    print("--- %s seconds ---" % (time.time() - start_time))
    return 



if __name__ == "__main__":
    annot_df_file = sys.argv[1]
    root_dir = sys.argv[2]
    config_file = sys.argv[3]
    
    annot_df = pd.read_csv(annot_df_file)

    main(annot_df, root_dir, config_file)