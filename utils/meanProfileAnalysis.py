import numpy as np
import scipy.spatial
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from random import sample, choices
from scipy.stats import pearsonr
import os
from sklearn import preprocessing

import sys

sys.path.insert(
    0, "/home/ubuntu/workspace_SingleCell/SingleCell_Morphological_Analysis/"
)

# from singlecell.read import read_single_cell_sql
from singlecell.preprocess import handle_nans, extract_cpfeature_names
from singlecell.process.normalize_funcs import *

# from singlecell.preprocess.filter_out_untransfected import extract_singlecell_transfection_labels

#####################################################################


def read_merge_preprocess_meanProfiles(annot_df_batches, read_pop_params):
    #     annot_df_batches=annot_df.copy()
    #     annot_df_batches=annot_df_sets.copy()
    feature_params = read_pop_params["feature_scaling_params_dict"]

    rootDir = read_pop_params["dirs_params_dict"]["rootDir"]
    profiles_folder_in_workspace = read_pop_params["dirs_params_dict"][
        "profiles_folder_in_workspace"
    ]

    batchNames = annot_df_batches["batch"].unique().tolist()

    df_batch_plates_annot_ls = []
    #     dfTransSummary_batch_ls=[]
    cp_features_analysis_batch_ls = []

    for bi in range(len(batchNames)):
        batchName = batchNames[bi]
        annot_df = annot_df_batches[annot_df_batches["batch"] == batchName].reset_index(
            drop=True
        )
        metaDataPlates = annot_df["Metadata_Plate"].unique()

        #         profiles_address=rootDir+'/backend/'+profType+'PerWells'+ndd+'/'    # for set1 and set2

        listOfPlates = annot_df["Metadata_Plate"].unique().tolist()

        df_batch_plates_ls = []
        for p in listOfPlates:
            #         profiles_address_old=rootDir+'/'+profiles_folder_in_workspace+'/'+batchName+'/'+p+'/population_profiles_'+p+'.csv.gz'
            profiles_address = (
                rootDir
                + "/"
                + profiles_folder_in_workspace
                + "/"
                + batchName
                + "/"
                + p
                + "/"
                + p
                + ".csv.gz"
            )

            #         os.rename(profiles_address_old, profiles_address)

            print(profiles_address)
            if os.path.exists(profiles_address):
                perplate_df = pd.read_csv(profiles_address)

                if feature_params["zscored_profiles"][0]:
                    if "zscored" in perplate_df.columns:
                        perplate_df = perplate_df[
                            perplate_df["zscored"]
                            == feature_params["zscored_profiles"][1]
                        ].reset_index(drop=True)
                    else:
                        raise ValueError(
                            "This batch doent have zscored column! "
                            "Probably no untransfected cells were assigned! check and carefully change"
                            "zscored_profiles enabled flag to False!"
                        )
                else:
                    if "zscored" in perplate_df.columns:
                        perplate_df = perplate_df[
                            perplate_df["zscored"].isnull()
                        ].reset_index(drop=True)

                perplate_df = perplate_df[
                    perplate_df["normalization"]
                    == feature_params["sc_per_plate_scaling"]
                ].reset_index(drop=True)

            df_batch_plates_ls.append(
                perplate_df[perplate_df["n_transf"] > 0].reset_index(drop=True)
            )

        df_batch = pd.concat(df_batch_plates_ls, ignore_index=True)

        (
            cp_features,
            cp_features_analysis,
        ) = extract_cpfeature_names.extract_cpfeature_names(df_batch)
        #   df_batch_plates, cp_features_analysis_batch = handle_nans.handle_nans(df,cp_features_analysis,\
        #                                                    thrsh_null_ratio=0.05,fill_na_method='drop-rows');

        df_batch_plates, cp_features_analysis_batch = handle_nans.handle_nans(
            df_batch, cp_features_analysis, thrsh_null_ratio=0
        )

        #         print('df_batch_plates1',df_batch_plates['Category'].unique())
        #         print("cp_features_analysis_batch",cp_features_analysis_batch)
        cp_features_analysis_batch_ls.append(cp_features_analysis_batch)

        df_batch_plates_annot_ls.append(df_batch_plates)

    df_annot_batches = pd.concat(df_batch_plates_annot_ls, ignore_index=True)
    # dfTransSummary_batches=pd.concat(dfTransSummary_batch_ls,ignore_index=True)

    feats_intersection_batches = cp_features_analysis_batch_ls[0]
    for bi in range(1, len(cp_features_analysis_batch_ls)):
        feats_intersection_batches = list(
            set(feats_intersection_batches) & set(cp_features_analysis_batch_ls[bi])
        )
        print(len(feats_intersection_batches), len(cp_features_analysis_batch_ls[bi]))

    if feature_params["post_scale_all_profiles"][0]:
        #             dataScaled=scaler.fit_transform(df2.loc[:,cpFeatures4scale])
        #             df_scaled[cpFeatures4scale]=dataScaled
        df_scaled_annot_batches = standardize_df_columns(
            df_annot_batches,
            feats_intersection_batches,
            feature_params["post_scale_all_profiles"][1],
        )
    else:
        df_scaled_annot_batches = df_annot_batches.copy()

    df_scaled_annot_batches, feats_intersection_batches = handle_nans.handle_nans(
        df_scaled_annot_batches, feats_intersection_batches, thrsh_null_ratio=0
    )

    ######################################################
    #     cpFeats_P=[p for p in feats_intersection_batches if 'Protein' in p]
    #     cpFeats_NP=[p for p in feats_intersection_batches if ('Protein' not in p) and ('DsRed' not in p)]

    prot_ch_suffix = read_pop_params["protein_channel_suffix"]

    cpFeats_P = [p for p in feats_intersection_batches if prot_ch_suffix in p]
    cpFeats_NP = [p for p in feats_intersection_batches if (prot_ch_suffix not in p)]

    return df_scaled_annot_batches, feats_intersection_batches, cpFeats_P, cpFeats_NP


def read_merge_preprocess_meanProfiles2(annot_df_batches, read_pop_params):
    #     annot_df=annot_df_3.copy()

    feature_params = read_pop_params["feature_scaling_params_dict"]

    rootDir = read_pop_params["dirs_params_dict"]["rootDir"]
    profiles_folder_in_backend = read_pop_params["dirs_params_dict"][
        "profiles_folder_in_backend"
    ]

    batchNames = annot_df_batches["batch"].unique().tolist()

    df_batch_plates_annot_ls = []
    dfTransSummary_batch_ls = []
    cp_features_analysis_batch_ls = []
    for bi in range(len(batchNames)):
        batchName = batchNames[bi]
        annot_df = annot_df_batches[annot_df_batches["batch"] == batchName].reset_index(
            drop=True
        )
        metaDataPlates = annot_df["Metadata_Plate"].unique()

        #         profiles_address=rootDir+'/backend/'+profType+'PerWells'+ndd+'/'    # for set1 and set2
        profiles_address = (
            rootDir + "/" + profiles_folder_in_backend + "/" + batchName + "/"
        )  # for last batch

        # listOfPlates0=os.listdir(rootDir+'/backend/wellsSingleCells/')[1:]
        listOfPlates0 = os.listdir(profiles_address)  # for set1 and set2
        #     listOfPlates0=os.listdir(rootDir+'/backend/population_profiles/'+batchName+'/') # for last batch
        ######################### load mean profiles

        #         strProf='df_'+normalization
        strProf = read_pop_params["dirs_params_dict"]["profiles_prefix_str"]

        listOfPlates = [
            p for p in listOfPlates0 if p.split(strProf)[1:] in metaDataPlates
        ]

        #         scaler_Plate= preprocessing.RobustScaler()
        feature_scaler = "Standard"

        #         df = pd.DataFrame();
        df_batch_plates_ls = []
        for p in listOfPlates:  # [0:1]: ['df_RC4_IF_05']:
            # for p in ['df_Replicate_3_s']:
            fileNameToSave = profiles_address + p
            #     if os.path.exists(fileNameToSave):
            print(fileNameToSave)
            transfectedMeanPerWell = pd.read_pickle(fileNameToSave, compression="infer")
            print("transfectedMeanPerWell", transfectedMeanPerWell.shape)
            transfectedMeanPerWell = transfectedMeanPerWell[
                transfectedMeanPerWell["tranfection_status"].isin([0, 1])
            ].reset_index(drop=True)
            print("transfectedMeanPerWell_transd", transfectedMeanPerWell.shape)
            #             print(transfectedMeanPerWell.shape,transfectedMeanPerWell[transfectedMeanPerWell['Variant'].isnull()].shape)
            #             print('transfectedMeanPerWell',transfectedMeanPerWell['Category'].unique())
            if transfectedMeanPerWell.shape[0] > 0:
                if feature_params["per_plate_scaling"][0]:
                    #                 if scaleMeanProfilesForEachPlate:
                    #                     print('hereeeee')
                    (
                        cp_features,
                        cp_features_analysis0,
                    ) = extract_cpfeature_names.extract_cpfeature_names(
                        transfectedMeanPerWell
                    )
                    (
                        transfectedMeanPerWell,
                        cp_features_analysis,
                    ) = handle_nans.handle_nans(
                        transfectedMeanPerWell, cp_features_analysis0
                    )
                    transfectedMeanPerWell = standardize_df_columns(
                        transfectedMeanPerWell,
                        cp_features_analysis,
                        feature_params["per_plate_scaling"][1],
                    )

                df_batch_plates_ls.append(transfectedMeanPerWell)

        df = pd.concat(df_batch_plates_ls, ignore_index=True)

        #                 df=df.append(transfectedMeanPerWell, ignore_index=True,sort=True)

        #         print('df_batch_plates0',df['Category'].unique())

        ############################ merge mean profiles and anotations
        (
            cp_features,
            cp_features_analysis,
        ) = extract_cpfeature_names.extract_cpfeature_names(df)
        #         df_batch_plates, cp_features_analysis_batch = handle_nans.handle_nans(df,cp_features_analysis,\
        #                                                    thrsh_null_ratio=0.05,fill_na_method='drop-rows');

        df_batch_plates, cp_features_analysis_batch = handle_nans.handle_nans(
            df, cp_features_analysis, thrsh_null_ratio=0
        )

        #         print('df_batch_plates1',df_batch_plates['Category'].unique())
        #         print("cp_features_analysis_batch",cp_features_analysis_batch)
        cp_features_analysis_batch_ls.append(cp_features_analysis_batch)
        #         cols2remove=[i for i in df.columns.tolist() if df[i].isnull().sum(axis=0)>.05*df.shape[0]]
        #     print(cols2remove)
        #         df2=df.drop(cols2remove, axis=1);

        if ("rep" in df_batch_plates.columns.tolist()) and (
            "rep" in annot_df.columns.tolist()
        ):
            df_batch_plates = df_batch_plates.drop("rep", axis=1)

        # cols2remove=[i for i in df2.columns.tolist() if df2[i].isnull().sum(axis=0)>0]
        # df2 = df2.fillna(df2.median())

        #         df2 = df2.interpolate()

        # df2_2=pd.merge(df2, annot_df, how='inner',on=['Metadata_Plate','Metadata_Well']);

        #         scaler = preprocessing.StandardScaler()
        #         df_scaled = df2.copy()

        #         cpFeatures=df2.columns[df2.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")].tolist()
        #         locFeature2beremoved=list(filter(lambda x: "_X" in x or "_Y" in x , cpFeatures))
        #         corFeature2beremoved=list(filter(lambda x: "Correlation" in x , cpFeatures))
        #         cpFeatures4scale=list(set(cpFeatures)-set(locFeature2beremoved)-set(corFeature2beremoved))
        #     cpFeatures4scale=list(set(cpFeatures)-set(locFeature2beremoved))

        #         if postScaleAllMeanProfiles:
        # #             dataScaled=scaler.fit_transform(df2.loc[:,cpFeatures4scale])
        # #             df_scaled[cpFeatures4scale]=dataScaled
        #             df_scaled = standardize_df_columns(df2,cp_features_analysis,feature_scaler)
        #         else:
        #             df_scaled=df2.copy()

        cols_to_use = annot_df.columns.difference(df_batch_plates.columns).tolist()
        print(cols_to_use)

        #         print('df_batch_plates',df_batch_plates['Category'].unique(),annot_df['Category'].unique())
        df_batch_plates_annot = pd.merge(
            df_batch_plates,
            annot_df[cols_to_use + ["batch", "Metadata_Plate", "Metadata_Well"]],
            how="inner",
            on=["batch", "Metadata_Plate", "Metadata_Well"],
        )

        #         print('left',pd.merge(df_batch_plates, annot_df[cols_to_use+['batch','Metadata_Plate','Metadata_Well']],\
        #                                  how='left',on=['batch','Metadata_Plate','Metadata_Well']).shape)
        #         print('inner',pd.merge(df_batch_plates, annot_df[cols_to_use+['batch','Metadata_Plate','Metadata_Well']],\
        #                                  how='inner',on=['batch','Metadata_Plate','Metadata_Well']).shape)

        #         print(df_batch_plates_annot['Category'].unique())
        # df4=df4.rename(columns={"Metadata_Location":"manual_Annot","rep_x":"rep"})

        dfTransSummary = df_batch_plates_annot[
            annot_df.columns.tolist() + ["n_transf", "n_untransf", "transf_Ratio"]
        ]

        dfTransSummary.loc[:, ["n_transf", "n_untransf", "transf_Ratio"]] = (
            dfTransSummary[["n_transf", "n_untransf", "transf_Ratio"]].fillna(0).values
        )

        dfTransSummary_batch_ls.append(dfTransSummary)
        df_batch_plates_annot_ls.append(df_batch_plates_annot)

    df_annot_batches = pd.concat(df_batch_plates_annot_ls, ignore_index=True)
    dfTransSummary_batches = pd.concat(dfTransSummary_batch_ls, ignore_index=True)

    feats_intersection_batches = cp_features_analysis_batch_ls[0]
    for bi in range(1, len(cp_features_analysis_batch_ls)):
        feats_intersection_batches = list(
            set(feats_intersection_batches) & set(cp_features_analysis_batch_ls[bi])
        )
        print(len(feats_intersection_batches), len(cp_features_analysis_batch_ls[bi]))

    if feature_params["post_scale_all_profiles"][0]:
        #             dataScaled=scaler.fit_transform(df2.loc[:,cpFeatures4scale])
        #             df_scaled[cpFeatures4scale]=dataScaled
        df_scaled_annot_batches = standardize_df_columns(
            df_annot_batches,
            feats_intersection_batches,
            feature_params["post_scale_all_profiles"][1],
        )
    else:
        df_scaled_annot_batches = df_annot_batches.copy()

    #     dfTransSummary[['n_transf','n_untransf','transf_Ratio']].fillna(value=0, inplace=True)

    ######################################################
    #     cpFeats_P=[p for p in feats_intersection_batches if 'Protein' in p]
    #     cpFeats_NP=[p for p in feats_intersection_batches if ('Protein' not in p) and ('DsRed' not in p)]

    cpFeats_P = [p for p in feats_intersection_batches if "GFP" in p]
    cpFeats_NP = [p for p in feats_intersection_batches if ("GFP" not in p)]

    return (
        df_scaled_annot_batches,
        feats_intersection_batches,
        cpFeats_P,
        cpFeats_NP,
        dfTransSummary_batches,
    )


def read_merge_EnrichmentProfiles(
    rootDir,
    annot_df_batches,
    scaleMeanProfilesForEachPlate,
    postScaleAllMeanProfiles,
    profile_version,
):
    #     annot_df=annot_df_3.copy()

    batchNames = annot_df_batches["batch"].unique().tolist()

    df_batch_plates_annot_ls = []
    #     dfTransSummary_batch_ls=[]
    #     cp_features_analysis_batch_ls=[]
    for bi in range(len(batchNames)):
        batchName = batchNames[bi]
        annot_df = annot_df_batches[annot_df_batches["batch"] == batchName].reset_index(
            drop=True
        )
        metaDataPlates = annot_df["Metadata_Plate"].unique()

        #         profiles_address=rootDir+'/backend/enrichment_profiles_sc_fit/'+batchName+'/' # for last batch
        profiles_address = (
            rootDir + "/backend/" + profile_version + "/" + batchName + "/"
        )  # for last batch

        listOfPlates0 = os.listdir(profiles_address)
        ######################### load mean profiles

        #         strProf='df_'+normalization
        strProf = "enrichment_profiles_"

        listOfPlates = [
            p for p in listOfPlates0 if p.split(strProf)[1:] in metaDataPlates
        ]

        #         scaler_Plate= preprocessing.RobustScaler()
        feature_scaler = "Standard"

        #         df = pd.DataFrame();
        df_batch_plates_ls = []
        for p in listOfPlates:  # [0:1]: ['df_RC4_IF_05']:
            # for p in ['df_Replicate_3_s']:
            fileNameToSave = profiles_address + p
            #     if os.path.exists(fileNameToSave):
            print(fileNameToSave)
            transfectedMeanPerWell = pd.read_pickle(fileNameToSave, compression="infer")
            if transfectedMeanPerWell.shape[0] > 0:
                if scaleMeanProfilesForEachPlate:
                    (
                        cp_features,
                        cp_features_analysis0,
                    ) = extract_cpfeature_names.extract_cpfeature_names(
                        transfectedMeanPerWell
                    )
                    (
                        transfectedMeanPerWell,
                        cp_features_analysis,
                    ) = handle_nans.handle_nans(
                        transfectedMeanPerWell, cp_features_analysis0
                    )
                    transfectedMeanPerWell = standardize_df_columns(
                        transfectedMeanPerWell, cp_features_analysis, feature_scaler
                    )

                df_batch_plates_ls.append(transfectedMeanPerWell)

        df_batch_plates = pd.concat(df_batch_plates_ls, ignore_index=True)

        #                 df=df.append(transfectedMeanPerWell, ignore_index=True,sort=True)

        ############################ merge mean profiles and anotations

        cp_features_analysis = df_batch_plates.columns[
            df_batch_plates.columns.str.contains("p_")
        ].tolist()

        #         cp_features_analysis_batch_ls.append(cp_features_analysis_batch)
        #         cols2remove=[i for i in df.columns.tolist() if df[i].isnull().sum(axis=0)>.05*df.shape[0]]
        #     print(cols2remove)
        #         df2=df.drop(cols2remove, axis=1);

        if ("rep" in df_batch_plates.columns.tolist()) and (
            "rep" in annot_df.columns.tolist()
        ):
            df_batch_plates = df_batch_plates.drop("rep", axis=1)

        cols_to_use = annot_df.columns.difference(df_batch_plates.columns).tolist()
        print(cols_to_use)

        df_batch_plates["batch"] = batchName
        #         print('df_batch_plates',df_batch_plates['Category'].unique(),annot_df['Category'].unique())
        df_batch_plates_annot = pd.merge(
            df_batch_plates,
            annot_df[cols_to_use + ["Metadata_Plate", "Metadata_Well"]],
            how="inner",
            on=["batch", "Metadata_Plate", "Metadata_Well"],
        )

        #         df_batch_plates_annot=pd.merge(df_batch_plates, annot_df[cols_to_use+['Metadata_Plate','Metadata_Well']],\
        #                                  how='left',on=['Metadata_Plate','Metadata_Well']);
        # df4=df4.rename(columns={"Metadata_Location":"manual_Annot","rep_x":"rep"})

        #         dfTransSummary=df_batch_plates_annot[annot_df.columns.tolist()+['n_transf','n_untransf','transf_Ratio']];

        #         dfTransSummary.loc[:,['n_transf','n_untransf','transf_Ratio']]=dfTransSummary[['n_transf','n_untransf','transf_Ratio']].fillna(0).values

        #         dfTransSummary_batch_ls.append(dfTransSummary)
        df_batch_plates_annot_ls.append(df_batch_plates_annot)

    df_annot_batches = pd.concat(df_batch_plates_annot_ls, ignore_index=True)
    #     dfTransSummary_batches=pd.concat(dfTransSummary_batch_ls,ignore_index=True)

    if postScaleAllMeanProfiles:
        #             dataScaled=scaler.fit_transform(df2.loc[:,cpFeatures4scale])
        #             df_scaled[cpFeatures4scale]=dataScaled
        df_scaled_annot_batches = standardize_df_columns(
            df_annot_batches, cp_features_analysis, feature_scaler
        )
    else:
        df_scaled_annot_batches = df_annot_batches.copy()

    #     dfTransSummary[['n_transf','n_untransf','transf_Ratio']].fillna(value=0, inplace=True)

    ######################################################
    #     cpFeats_P=[p for p in feats_intersection_batches if 'Protein' in p]
    #     cpFeats_NP=[p for p in feats_intersection_batches if ('Protein' not in p) and ('DsRed' not in p)]

    cpFeats_NP = ["np_" + str(i) for i in range(20)]
    cpFeats_P = ["p_" + str(i) for i in range(20)]

    return df_scaled_annot_batches, cp_features_analysis, cpFeats_P, cpFeats_NP


def categoricalRepCor(inDf, pertColName, featColNames, plotEnabled):
    """
    Calculates replicate correlation versus across purtburtion correlations

    This function takes the input dataframe and output/plot replicate correlations.

    Parameters:
    inDf   (pandas df): input dataframe contains metadata and features
    pertColName  (str): The column based on which we define replicates of a purturbation
    featColNames(list): The column based on which we define replicates of a purturbation
    plotEnabled (bool): If True or 1, plots the curves

    Returns:
    list: [repC,randC]

    """
    df = inDf.copy()
    uniqPert = df[pertColName].unique().tolist()
    repC = []
    #     df_repCor= pd.DataFrame(columns=[pertColName,'R1','R2','Plate1','Plate2','cc'])
    for u in uniqPert:
        df1 = df[df[pertColName] == u].drop_duplicates().reset_index(drop=True)
        #         df2=df[df[pertColName]!=u].drop_duplicates().reset_index(drop=True)

        repCorrPurtbs = df1.loc[:, featColNames].T.corr()
        triuIndices = np.triu_indices(repCorrPurtbs.shape[0], k=1)
        repCorr = list(repCorrPurtbs.values[triuIndices])

        #         repCorr=np.sort(np.unique(df1.loc[:,featColNames].T.corr().values))[:-1].tolist()
        #         repC=repC+repCorr
        #         print(np.median(repCorr))
        repC = repC + [np.median(repCorr)]

    ##################### correlation with random
    #     print(randC)
    randC_v2 = []
    for i in range(10):
        uniqeSamplesFromEachPurt = inDf.groupby(pertColName)[featColNames].apply(
            lambda s: s.sample(1)
        )
        corrMatAcrossPurtbs = uniqeSamplesFromEachPurt.loc[:, featColNames].T.corr()
        randCorrVals = list(
            corrMatAcrossPurtbs.values[
                np.triu_indices(corrMatAcrossPurtbs.shape[0], k=1)
            ]
        )
    randC_v2 = randC_v2 + randCorrVals

    return repC, randC_v2


#####################################################################


def generateFeatureGroupsMap_medianZscore(MutF, WtF, saveResultDir, Channelss):
    #     Channelss=['Protein','Mito','ER','DNA']
    featureGroups = ["Texture", "Intensity", "RadialDistribution"]

    compartmentsStr = "Cells_|Nuclei_|Cytoplasm_"
    Mat = np.zeros((len(featureGroups), len(Channelss)))
    scaler = preprocessing.StandardScaler()

    for ch in range(len(Channelss)):
        for f in range(len(featureGroups)):
            xM = MutF.loc[
                :,
                MutF.columns.str.contains(Channelss[ch])
                & MutF.columns.str.contains(featureGroups[f])
                & MutF.columns.str.contains(compartmentsStr),
            ]
            xW = WtF.loc[
                :,
                WtF.columns.str.contains(Channelss[ch])
                & WtF.columns.str.contains(featureGroups[f])
                & WtF.columns.str.contains(compartmentsStr),
            ]
            #             print(Channelss[ch],featureGroups[f],xW.shape,xM.shape)
            if xW.shape[1] > 0 and xM.shape[1] > 0:
                s0 = scaler.fit(xW.values)
                zScored = s0.transform(xM.values)
                zScored_averageAcrossSamples = np.mean(zScored, 0)
                #                 print(zScored_averageAcrossSamples.shape)
                Mat[f, ch] = np.median(abs(zScored_averageAcrossSamples))
    #                 Mat[f,ch]=np.median(zScored_averageAcrossSamples);
    #             if xM.shape[0]>1:
    #                 xM=xM.median();
    #             diff=xM.values-xW.values;
    #             print(xM.shape,xW.shape)
    #             Mat[f,ch]=np.linalg.norm(diff, 2);

    #     print(Mat)
    fig, axes = plt.subplots(2, 1, gridspec_kw={"height_ratios": [1, 3]})
    sns.heatmap(Mat, cmap=sns.cubehelix_palette(8), ax=axes[1])
    axes[1].xaxis.set_ticklabels(Channelss)
    axes[1].yaxis.set_ticklabels(featureGroups, rotation=360 + 45)
    plt.tight_layout()

    cellAreas = ["Cells", "Nuclei", "Cytoplasm"]
    MatArea = np.zeros((1, len(cellAreas)))
    for ca in range(len(cellAreas)):
        xM = MutF.loc[
            :,
            MutF.columns.str.contains("AreaShape")
            & MutF.columns.str.contains(cellAreas[ca]),
        ]
        xW = WtF.loc[
            :,
            WtF.columns.str.contains("AreaShape")
            & WtF.columns.str.contains(cellAreas[ca]),
        ]
        #         if xM.shape[0]>1:
        #             xM=xM.median();
        #         MatArea[0,ca]=np.linalg.norm(xM.values-xW.values, 2);

        if xW.shape[1] > 0 and xM.shape[1] > 0:
            s0 = scaler.fit(xW.values)
            zScored = s0.transform(xM.values)
            zScored_averageAcrossSamples = np.mean(zScored, 0)
            #             print(zScored_averageAcrossSamples.shape)
            MatArea[0, ca] = np.median(abs(zScored_averageAcrossSamples))

    #     print(MatArea)
    sns.heatmap(MatArea, cmap=sns.cubehelix_palette(8), ax=axes[0])
    axes[0].xaxis.tick_top()
    axes[0].xaxis.set_ticklabels(cellAreas, rotation=360 + 45)
    axes[0].yaxis.set_ticklabels(["AreaShape"], rotation=360 + 45)
    plt.tight_layout()
    if 1:
        fig.savefig(saveResultDir + "/featureCategoryImportance-medianAbsZscore.png")
        plt.close("all")
    return Mat, MatArea


#######################################################


def generateFeatureGroupsMap_l2ofDiff(MutF, WtF, saveResultDir, Channelss):
    #     Channelss=['Protein','Mito','ER','DNA']
    featureGroups = ["Texture", "Intensity", "RadialDistribution"]

    compartmentsStr = "Cells_|Nuclei_|Cytoplasm_"
    Mat = np.zeros((len(featureGroups), len(Channelss)))
    scaler = preprocessing.StandardScaler()

    for ch in range(len(Channelss)):
        for f in range(len(featureGroups)):
            xM = MutF.loc[
                :,
                MutF.columns.str.contains(Channelss[ch])
                & MutF.columns.str.contains(featureGroups[f])
                & MutF.columns.str.contains(compartmentsStr),
            ]
            xW = WtF.loc[
                :,
                WtF.columns.str.contains(Channelss[ch])
                & WtF.columns.str.contains(featureGroups[f])
                & WtF.columns.str.contains(compartmentsStr),
            ]
            #             print(Channelss[ch],featureGroups[f],xW.shape,xM.shape)
            if xW.shape[1] > 0 and xM.shape[1] > 0:
                Mat[f, ch] = np.linalg.norm(xM.mean().values - xW.mean().values, 2)

    #     print(Mat)
    fig, axes = plt.subplots(2, 1, gridspec_kw={"height_ratios": [1, 3]})
    sns.heatmap(Mat, cmap=sns.cubehelix_palette(8), ax=axes[1])
    axes[1].xaxis.set_ticklabels(Channelss)
    axes[1].yaxis.set_ticklabels(featureGroups, rotation=360 + 45)
    plt.tight_layout()

    cellAreas = ["Cells", "Nuclei", "Cytoplasm"]
    MatArea = np.zeros((1, len(cellAreas)))
    for ca in range(len(cellAreas)):
        xM = MutF.loc[
            :,
            MutF.columns.str.contains("AreaShape")
            & MutF.columns.str.contains(cellAreas[ca]),
        ]
        xW = WtF.loc[
            :,
            WtF.columns.str.contains("AreaShape")
            & WtF.columns.str.contains(cellAreas[ca]),
        ]
        #         if xM.shape[0]>1:
        #             xM=xM.median();
        #         MatArea[0,ca]=np.linalg.norm(xM.values-xW.values, 2);

        if xW.shape[1] > 0 and xM.shape[1] > 0:
            MatArea[0, ca] = np.linalg.norm(xM.mean().values - xW.mean().values, 2)

    #     print(MatArea)
    sns.heatmap(MatArea, cmap=sns.cubehelix_palette(8), ax=axes[0])
    axes[0].xaxis.tick_top()
    axes[0].xaxis.set_ticklabels(cellAreas, rotation=360 + 45)
    axes[0].yaxis.set_ticklabels(["AreaShape"], rotation=360 + 45)
    plt.tight_layout()
    if 1:
        fig.savefig(saveResultDir + "/featureCategoryImportance-l2ofDiff.png")
        plt.close("all")
    return Mat, MatArea


#######################################################


def featureImportanceEmbeded(embeddingDf, MutF, WtF, saveResultDir):
    cat2a, cat1a = MutF[embeddingDf.columns], WtF[embeddingDf.columns]
    if cat1a.shape[0] > 1:
        cat1a = cat1a.median()

    if cat2a.shape[0] > 1:
        cat2a = cat2a.median()

    embeddingDf.loc[2, :] = cat2a.values - cat1a.values
    embeddingDfT = embeddingDf.T
    embeddingDfT.loc[:, 3] = embeddingDfT.loc[:, 2].abs()
    absFeatureImportance = embeddingDfT.sort_values(3, ascending=False)[:20]

    xs, ys = (
        absFeatureImportance.loc[:, 0].values,
        absFeatureImportance.loc[:, 1].values,
    )
    #     print(ys)
    #     ys=10*(ys+abs(np.min(ys)))
    #     ys=np.power(ys,.1)
    #     print(ys)
    #     ind=np.argsort(ys)

    ys_2 = np.sort(ys)
    ind2 = [np.where(ys[i] == ys_2)[0][0] for i in range(len(ys))]
    for i in range(1, len(ys_2)):
        if ys_2[i] - ys_2[i - 1] < 0.06:
            ys_2[i:] = ys_2[i:] + 0.06
    ys = ys_2[ind2]

    maxAbsFeatImp = absFeatureImportance.loc[:, 3].max()
    ys2 = np.copy(ys)
    names = absFeatureImportance.index.tolist()
    fig, axes = plt.subplots(figsize=(9, 7))
    count = 0
    for x, y, name in zip(xs, ys, names):
        colorr = "blue"
        if embeddingDfT.loc[name, 2] < 0:
            colorr = "red"

        axes.scatter(x, y, c=colorr)
        axes.text(
            x,
            y,
            name,
            color=colorr,
            fontsize=6 + ((embeddingDfT.loc[name, 3] / maxAbsFeatImp) * 4),
        )
    #              linespacing=xx.loc[name,3]*10)
    #     plt.tight_layout(w_pad=10)
    axes.set_xlim([np.min(xs) - 0.1, np.max(xs) + 1.5])
    axes.xaxis.set_major_locator(plt.NullLocator())
    axes.yaxis.set_major_locator(plt.NullLocator())
    saveFormat = ".png"
    #'.png'
    if 1:
        plt.tight_layout()
        fig.savefig(saveResultDir + "/featureImportance" + saveFormat)
        #         os.system('rm ../typeA-B/WT-MUT-'+namess+'/featureImportance.eps')
        plt.close("all")

    return
