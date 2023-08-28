"""
@author: mhaghigh
"""

import pandas as pd
import numpy as np
import seaborn as sns

# from sklearn import preprocessing
import sklearn.preprocessing as sp
import pickle

# from imblearn.over_sampling import SMOTE
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, confusion_matrix
from sklearn.cluster import AgglomerativeClustering, SpectralClustering, KMeans
import matplotlib.pyplot as plt
import os
import time
import gc
from functools import reduce

import sys

sys.path.insert(
    0, "/home/ubuntu/workspace_SingleCell/SingleCell_Morphological_Analysis/"
)

from singlecell.read import read_single_cell_sql
from singlecell.preprocess import handle_nans, extract_cpfeature_names
from singlecell.process.normalize_funcs import *
from singlecell.preprocess.filter_out_untransfected import (
    extract_singlecell_transfection_labels,
)
from singlecell.preprocess.filter_out_edge_single_cells import edgeCellFilter
from singlecell.save.save_pandas_dfs import saveDF_to_CSV_GZ_no_timestamp

################################################################################


################################################################################
def enrichment_clustering_model_random_single_cell(annot_df, rootPath, params):
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
    transfection_params_dict = params["transfection_params_dict"]
    feature_scaler = params["feature_scaling_params_dict"]["feature_scaler"]

    listOfBatchPlates = annot_df.Metadata_batch_Plate.unique().tolist()

    trans_sc_allPlates = []
    for batchPlateName in listOfBatchPlates:
        batchName = annot_df[annot_df["Metadata_batch_Plate"] == batchPlateName][
            "batch"
        ].tolist()[0]
        plateName = annot_df[annot_df["Metadata_batch_Plate"] == batchPlateName][
            "Metadata_Plate"
        ].tolist()[0]

        annot_df_plate = annot_df[
            annot_df["Metadata_batch_Plate"] == batchPlateName
        ].reset_index(drop=True)

        ## path for sqlite files outputs of cell profiler and cytominer packages
        fileName = (
            rootPath
            + "/backend/"
            + batchName
            + "/"
            + plateName
            + "/"
            + plateName
            + ".sqlite"
        )

        ## read single cell sql files
        start_time = time.time()
        sc_df = read_single_cell_sql.readSingleCellData_sqlalch(
            fileName, ["cells", "cytoplasm", "nuclei"]
        )

        ## filter features with low std and nans for than a threshold
        (
            cp_features,
            cp_features_analysis_0,
        ) = extract_cpfeature_names.extract_cpfeature_names(sc_df)
        sc_df, cp_features_analysis = handle_nans.handle_nans(
            sc_df, cp_features_analysis_0, fill_na_method="drop-rows"
        )
        print(sc_df[cp_features_analysis].shape)

        ## filter cells on the edge
        sc_df, _ = edgeCellFilter(sc_df)
        print("edge removed: ", sc_df[cp_features_analysis].shape)

        ## scale features per plate for analysis
        standardized_per_plate_sc = standardize_df_columns(
            sc_df, cp_features_analysis, feature_scaler
        )

        if transfection_params_dict:
            transfection_labels = extract_singlecell_transfection_labels(
                sc_df, transfection_params_dict
            )

            sc_df["transfection_status"] = transfection_labels
            standardized_per_plate_sc["transfection_status"] = transfection_labels

        else:
            sc_df["transfection_status"] = 1
            standardized_per_plate_sc["transfection_status"] = 1

        # create enrichement profiles

        transfected_sc_df = (
            sc_df[sc_df["transfection_status"] == 1].reset_index(drop=True).sample(1000)
        )

        trans_sc_allPlates.append(transfected_sc_df[cp_features_analysis])

    return trans_sc_allPlates


################################################################################
def generate_population_profiles(batchPlateName, annot_df, rootPath, params):
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
    transfection_params_dict = params["transfection_params_dict"]
    enrichement_profiles_params = params["enrichement_profiles_params"]
    feature_scaler = params["feature_scaling_params_dict"]["feature_scaler"]

    batchName = annot_df[annot_df["Metadata_batch_Plate"] == batchPlateName][
        "batch"
    ].tolist()[0]
    plateName = annot_df[annot_df["Metadata_batch_Plate"] == batchPlateName][
        "Metadata_Plate"
    ].tolist()[0]

    annot_df_plate = annot_df[
        annot_df["Metadata_batch_Plate"] == batchPlateName
    ].reset_index(drop=True)

    print("Processing ", batchName, plateName)

    ## Set path to save pupulation profiles
    save_population_profiles_path = (
        rootPath + "/population_profiles/" + batchName + "/" + plateName
    )
    save_sc_profiles_path = (
        rootPath + "/singlecell_profiles/" + batchName + "/" + plateName
    )
    save_enrichment_profiles_path = (
        rootPath + "/enrichment_profiles/" + batchName + "/" + plateName
    )

    ## path for sqlite files outputs of cell profiler and cytominer packages
    fileName = (
        rootPath
        + "/backend/"
        + batchName
        + "/"
        + plateName
        + "/"
        + plateName
        + ".sqlite"
    )

    if not os.path.exists(fileName):
        print(fileName + " do not exist!")
        return

    os.system("mkdir -p " + save_population_profiles_path)
    os.system("mkdir -p " + save_sc_profiles_path)
    os.system("mkdir -p " + save_enrichment_profiles_path)

    ## read single cell sql files
    start_time = time.time()
    sc_df = read_single_cell_sql.readSingleCellData_sqlalch(
        fileName, ["cells", "cytoplasm", "nuclei"]
    )

    ## filter features with low std and nans for than a threshold
    (
        cp_features,
        cp_features_analysis_0,
    ) = extract_cpfeature_names.extract_cpfeature_names(sc_df)
    sc_df, cp_features_analysis = handle_nans.handle_nans(
        sc_df, cp_features_analysis_0, fill_na_method="drop-rows"
    )
    print(sc_df[cp_features_analysis].shape)

    ## filter cells on the edge
    sc_df, _ = edgeCellFilter(sc_df)
    print(sc_df[cp_features_analysis].shape)

    ## scale features per plate for analysis
    standardized_per_plate_sc = standardize_df_columns(
        sc_df, cp_features_analysis, feature_scaler
    )

    ## apply transfection filtering of single cell level of data based on the parameters
    # transfection_params_dict={'Method':'single_intensity_feature_thrsh',\
    #                               'intensity_feature_to_use':'Cells_Intensity_MeanIntensity_DsRed',\
    #                               'thresholding_method': 'precomputed_batch_specific_thrsh',\
    #                               'precomputed_params': [precomputed_bottom_thresh , precomputed_top_thresh , data_norm_used]}

    if transfection_params_dict:
        #         if "precomputed_params" in transfection_params_dict.keys() and transfection_params_dict["precomputed_params"]=='raw':

        transfection_labels = extract_singlecell_transfection_labels(
            sc_df, transfection_params_dict
        )

        sc_df["transfection_status"] = transfection_labels
        standardized_per_plate_sc["transfection_status"] = transfection_labels

    else:
        sc_df["transfection_status"] = 1
        standardized_per_plate_sc["transfection_status"] = 1

    if params["save_single_cells"]:
        transfected_sc_df = sc_df[sc_df["transfection_status"] == 1].reset_index(
            drop=True
        )
        standardized_transfected_sc_df = standardize_df_columns(
            transfected_sc_df, cp_features_analysis, feature_scaler
        )

        wells = standardized_transfected_sc_df["Metadata_Well"].unique().tolist()

        for w in wells:
            fileNameToSave = save_sc_profiles_path + "/" + plateName + "_" + w
            #             standardized_transfected_sc_df[standardized_transfected_sc_df['Metadata_Well']==w].\
            #             reset_index(drop=True).to_pickle(fileNameToSave);
            saveDF_to_CSV_GZ_no_timestamp(
                standardized_transfected_sc_df[
                    standardized_transfected_sc_df["Metadata_Well"] == w
                ].reset_index(drop=True),
                fileNameToSave,
            )

    # create enrichement profiles
    if enrichement_profiles_params:
        #     clustering_model_file = './utils/kmeans_20clusters_fitted_on_population_profiles.sav'
        #         clustering_model_file = './utils/kmeans_20clusters_fitted_on_sc_profiles_3.sav'
        clustering_model_file = enrichement_profiles_params["clustering_model_file"]
        clustering_model_dict = pickle.load(open(clustering_model_file, "rb"))

        transfected_sc_df = sc_df[sc_df["transfection_status"] == 1].reset_index(
            drop=True
        )

        #     transfected_standardized_per_plate_sc=standardized_per_plate_sc\
        #     [standardized_per_plate_sc['transfection_status']==1].reset_index(drop=True)
        enrichment_profile_df = (
            transfected_sc_df.groupby("Metadata_Well")
            .size()
            .reset_index()
            .rename(columns={0: "count"})
        )

        enrichment_profile_df["Metadata_Plate"] = plateName

        for fi in ["p", "np"]:
            n_clusts = clustering_model_dict["n_clusts"]
            cpFeats = clustering_model_dict["model_" + fi + "_feats"]
            model_cls = clustering_model_dict["model_" + fi]
            scaler_pre_cls = clustering_model_dict["scaler_" + fi]
            enrichment_profile_df_tmp = (
                transfected_sc_df.groupby("Metadata_Well")
                .apply(
                    lambda x: fill_enrichment_profile(
                        np.unique(
                            model_cls.predict(scaler_pre_cls.transform(x[cpFeats])),
                            return_counts=True,
                        ),
                        n_clusts,
                    )
                )
                .reset_index()
            )
            enrichment_profile_df = pd.concat(
                [
                    enrichment_profile_df,
                    enrichment_profile_df_tmp.apply(
                        lambda x: x[0], axis=1, result_type="expand"
                    ).add_prefix(fi + "_"),
                ],
                axis=1,
            )

        fileNameToSave = (
            save_enrichment_profiles_path + "/enrichment_profiles_" + plateName
        )
        saveDF_to_CSV_GZ_no_timestamp(enrichment_profile_df, fileNameToSave)
    #         enrichment_profile_df.to_pickle(fileNameToSave);

    ## count per-well transfection summaries
    transfection_summary = (
        sc_df[["Metadata_Plate", "Metadata_Well", "transfection_status"]]
        .groupby(["Metadata_Plate", "Metadata_Well"])
        .apply(
            lambda r: pd.Series(
                {
                    "n_transf": sum(r["transfection_status"] == 1),
                    "n_untransf": sum(r["transfection_status"] != 1),
                    "transf_Ratio": sum(r["transfection_status"] == 1) / r.shape[0],
                }
            )
        )
        .reset_index()
    )

    save_population_profiles_zscored_by_same_well_untransfected = False
    ## save population profiles z scored by same well untransfected cells
    if save_population_profiles_zscored_by_same_well_untransfected:
        print("sfds")
        standardized_per_plate_sc.groupby("Metadata_Well").apply(
            lambda g: zscore_by_untrasfected(g)
        )

    ## generate population level profiles
    population_profiles = (
        sc_df.groupby(["Metadata_Plate", "Metadata_Well", "transfection_status"])
        .mean()
        .reset_index()
    )
    print(sc_df[cp_features_analysis].shape, population_profiles.shape)
    population_profiles_standardized = (
        standardized_per_plate_sc.groupby(
            ["Metadata_Plate", "Metadata_Well", "transfection_status"]
        )
        .mean()
        .reset_index()
    )

    print(sc_df[cp_features_analysis].shape)
    print(population_profiles[cp_features_analysis].shape)

    ## augment profiles with metadata and transfection summary
    population_prof_aug = reduce(
        lambda left, right: pd.merge(
            left, right, on=["Metadata_Plate", "Metadata_Well"], how="outer"
        ),
        [population_profiles, transfection_summary, annot_df_plate],
    )
    population_prof_standardized_aug = reduce(
        lambda left, right: pd.merge(
            left, right, on=["Metadata_Plate", "Metadata_Well"], how="outer"
        ),
        [population_profiles_standardized, transfection_summary, annot_df_plate],
    )

    ## save raw averaged profiles together with perplate scaled and then per well averaged profiles
    population_prof_aug["normalization"] = "raw"
    population_prof_standardized_aug["normalization"] = "sc_scaled_per_plate"

    profile_list_2save = [population_prof_aug, population_prof_standardized_aug]

    ## generate zscored population profiles to untransfected population profiles

    if transfection_params_dict:
        population_prof_aug_zscored = generate_zscored_to_untransfected_profiles(
            population_prof_aug, cp_features_analysis, feature_scaler
        )
        population_prof_standardized_aug_zscored = (
            generate_zscored_to_untransfected_profiles(
                population_prof_standardized_aug, cp_features_analysis, feature_scaler
            )
        )

        profile_list_2save += [
            population_prof_aug_zscored,
            population_prof_standardized_aug_zscored,
        ]

    population_profiles_2save = pd.concat(profile_list_2save, ignore_index=True)

    #     fileNameToSave=save_population_profiles_path+'/population_profiles_'+plateName;
    fileNameToSave = save_population_profiles_path + "/" + plateName
    saveDF_to_CSV_GZ_no_timestamp(population_profiles_2save, fileNameToSave)
    #     population_profiles_2save.to_pickle(fileNameToSave);

    gc.collect()

    print("--- %s seconds ---" % (time.time() - start_time))
    return


def fill_enrichment_profile(unique_counts_tuple, n_clusts):
    enrichment_profile = np.zeros((n_clusts), dtype=int)
    enrichment_profile[unique_counts_tuple[0]] = unique_counts_tuple[1]
    return enrichment_profile


def generate_zscored_to_untransfected_profiles(
    population_prof_aug, cp_features_analysis, scaling_method="Standard"
):
    sc_control_df = population_prof_aug.loc[
        population_prof_aug["transfection_status"] != 1
    ].reset_index(drop=True)
    sc_control_df_stringent = population_prof_aug.loc[
        population_prof_aug["transfection_status"] == 0
    ].reset_index(drop=True)

    sc_df = population_prof_aug.loc[
        population_prof_aug["transfection_status"] == 1
    ].reset_index(drop=True)

    #     zscore_df_columns_by_control(sc_control_df, sc_df,cp_features_analysis,scaling_method):
    population_prof_aug_zscored = zscore_df_columns_by_control(
        sc_control_df, sc_df, cp_features_analysis, scaling_method
    )

    population_prof_aug_zscored["zscored"] = "untransfected"

    population_prof_aug_zscored_stringent = zscore_df_columns_by_control(
        sc_control_df_stringent, sc_df, cp_features_analysis, scaling_method
    )

    population_prof_aug_zscored_stringent["zscored"] = "untransfected_stringent"

    return pd.concat(
        [population_prof_aug_zscored, population_prof_aug_zscored_stringent],
        ignore_index=True,
    )


# -----------------------------------------------------------------------------------------------------
def saveRawIntensityFeatures(batchPlateName, annot_df, rootPath, channels_used):
    from sqlalchemy import create_engine

    #     plateName=listOfPlates[0]
    #     batchName=annot_df[annot_df['Metadata_Plate']==plateName]['batch'].tolist()[0]

    batchName = annot_df[annot_df["Metadata_batch_Plate"] == batchPlateName][
        "batch"
    ].tolist()[0]
    plateName = annot_df[annot_df["Metadata_batch_Plate"] == batchPlateName][
        "Metadata_Plate"
    ].tolist()[0]

    print(batchName, plateName)

    start_time = time.time()
    fileName = (
        rootPath
        + "/backend/"
        + batchName
        + "/"
        + plateName
        + "/"
        + plateName
        + ".sqlite"
    )

    #     print(os.path.exists(fileName))
    if not os.path.exists(fileName):
        print(fileName + " do not exist!")
        return

    selected_fs = [
        "_Intensity_IntegratedIntensity",
        "_Intensity_MeanIntensity",
        "_Intensity_StdIntensity",
        "_Intensity_MaxIntensity",
        "_Intensity_UpperQuartileIntensity",
    ]

    compartments = ["Cells", "Nuclei"]

    #     channels_used=['Protein','DsRed']
    selected_features = [f + "_" + ch for ch in channels_used for f in selected_fs]

    sql_file = "sqlite:////" + fileName
    engine = create_engine(sql_file)
    conn = engine.connect()

    # compartments=selected_feature.split("_")[0]
    #     compartments='Cells'
    plateDf_list = []
    for comp in compartments:
        selected_features_comp = [comp + f for f in selected_features]
        selected_feature = " ,".join(selected_features_comp)
        query_cols = (
            "TableNumber, ImageNumber, ObjectNumber, " + selected_feature
        )  # +", "+f2
        #         print(query_cols)
        compartment_query = "select {} from {}".format(query_cols, comp)
        plateDf_list.append(pd.read_sql(sql=compartment_query, con=conn))

    plateDf = reduce(
        lambda left, right: pd.merge(
            left, right, on=["TableNumber", "ImageNumber", "ObjectNumber"]
        ),
        plateDf_list,
    )
    #     plateDf=plateDf.dropna()

    img_query = "select * from {}".format("Image")
    plateImageDf = pd.read_sql(sql=img_query, con=conn)

    plateDfwMeta = pd.merge(plateDf, plateImageDf, on=["TableNumber", "ImageNumber"])
    #     print(plateDfwMeta.columns)

    selected_features_all = [
        c + f + "_" + ch
        for c in compartments
        for ch in channels_used
        for f in selected_fs
    ]

    transfectionDetectionFeatures = plateDfwMeta.loc[
        :,
        selected_features_all
        + ["Metadata_batch_Plate", "Metadata_Plate", "Metadata_Well"],
    ]
    #     os.system("mkdir "+rootPath+'/plate_raw_intensity_features/'+batchName)
    os.makedirs(rootPath + "/plate_raw_intensity_features/" + batchName, exist_ok=True)
    transfectionDetectionFeatures.to_pickle(
        rootPath
        + "/plate_raw_intensity_features/"
        + batchName
        + "/df_intensityFeatures_"
        + plateName
    )
    print(
        "Reading sql, saving int df --- %s minutes ---"
        % ((time.time() - start_time) / 60)
    )
    return


def readNormalizeIntensityFeatures(annot_df, rootDir, channels_used):
    listOfBatchPlates = annot_df.Metadata_batch_Plate.unique().tolist()
    scaler_minmax = sp.MinMaxScaler(feature_range=(0, 1))

    df_inten = pd.DataFrame()
    df_inten_scaled_perPlate = pd.DataFrame()

    for bp in listOfBatchPlates:  # [0:1]:
        batch = annot_df[annot_df["Metadata_batch_Plate"] == bp]["batch"].tolist()[0]
        p = annot_df[annot_df["Metadata_batch_Plate"] == bp]["Metadata_Plate"].tolist()[
            0
        ]

        fileNameToSave = (
            rootDir
            + "/plate_raw_intensity_features/"
            + batch
            + "/df_intensityFeatures_"
            + p
        )
        if os.path.exists(fileNameToSave):
            intFeaturesDf = pd.read_pickle(fileNameToSave, compression="infer")
            df_inten = df_inten.append(intFeaturesDf, ignore_index=True)
            df_inten_scaled0 = intFeaturesDf.copy()
            intFeatures = intFeaturesDf.columns[
                intFeaturesDf.columns.str.contains("|".join(channels_used))
            ].tolist()
            #         sfsdfdsfds
            for ifi in intFeatures:
                qpi_up = intFeaturesDf[ifi].quantile(0.999)
                qpi_low = intFeaturesDf[ifi].quantile(0.01)
                intFeaturesDf[ifi] = intFeaturesDf[ifi].clip(qpi_low, qpi_up)

            #     dataScaled=scaler0.fit_transform(intFeaturesDf.loc[:,intFeatures])
            df_inten_scaled0[intFeatures] = scaler_minmax.fit_transform(
                intFeaturesDf.loc[:, intFeatures]
            )
            df_inten_scaled_perPlate = df_inten_scaled_perPlate.append(
                df_inten_scaled0, ignore_index=True
            )

    #     print(df_inten.shape)

    df_inten = pd.merge(
        df_inten, annot_df, how="inner", on=["Metadata_Plate", "Metadata_Well"]
    )
    df_inten_scaled_perPlate = pd.merge(
        df_inten_scaled_perPlate,
        annot_df,
        how="inner",
        on=["Metadata_Plate", "Metadata_Well"],
    )

    return df_inten, df_inten_scaled_perPlate


def plot_intensity_features(listOfBatches, df_ls, params):
    from matplotlib import rcParams

    rcParams["patch.force_edgecolor"] = False

    intFeatures = params["intFeature"]

    fig, axes = plt.subplots(
        2,
        len(listOfBatches),
        figsize=(len(listOfBatches) * 5, 2 * 5),
        sharex=params["sharex_enabled"],
    )
    # fig.suptitle('_'.join(intFeatures[i].split('_')[2:]));

    for bi in range(len(listOfBatches)):
        for si, dfi, sn in zip([0, 1], df_ls, ["raw", "normalized"]):
            #     for si,dfi in enumerate([df_inten,df_inten_scaled_perPlate]):
            df_inten_b = dfi[
                (dfi["batch"] == listOfBatches[bi])
                & (~dfi[params["hue_column"]].isnull())
            ].reset_index(drop=True)
            #### Method 1 -a

            if params["hue_column"]:
                qpi_up = (
                    df_inten_b.loc[
                        df_inten_b[params["hue_column"]] == True, intFeatures
                    ]
                    .quantile(0.99)
                    .round(4)
                )
                qpi_low = (
                    df_inten_b.loc[
                        df_inten_b[params["hue_column"]] == True, intFeatures
                    ]
                    .quantile(0.4)
                    .round(4)
                )

                print(sn, listOfBatches[bi], qpi_low, qpi_up)

            if df_inten_b[intFeatures].min() == 0:
                df_inten_b[intFeatures] = (
                    df_inten_b[intFeatures] + params["zero_offset_4logscale_plot"]
                )
                qpi_up = df_inten_b.loc[
                    df_inten_b[params["hue_column"]] == True, intFeatures
                ].quantile(0.99)

            if params["hue_column"]:
                sns.histplot(
                    data=df_inten_b,
                    x=intFeatures,
                    bins=params["hist_n_bins"],
                    stat="density",
                    hue=params["hue_column"],
                    element="step",
                    common_norm=False,
                    legend=True,
                    log_scale=params["log_scale_enabled"],
                    ax=axes[si, bi],
                )
                #             perc95=qpi_up+zero_offset_4logscale_plot
                axes[si, bi].axvline(x=qpi_up, linestyle=":", color="r")
            else:
                sns.histplot(
                    data=df_inten_b,
                    x=intFeatures,
                    bins=params["hist_n_bins"],
                    stat="density",
                    element="step",
                    common_norm=False,
                    legend=True,
                    log_scale=params["log_scale_enabled"],
                    ax=axes[si, bi],
                )

        axes[0, bi].set_title(listOfBatches[bi])
    axes[0, 0].set_ylabel("raw")
    axes[1, 0].set_ylabel("scaled")
    return fig
