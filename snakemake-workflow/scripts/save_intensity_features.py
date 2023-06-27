import pandas as pd
import os
import sys
from sqlalchemy import create_engine
from functools import reduce
import time

def saveRawIntensityFeatures(batchPlateName, annot_df, rootPath, channels_used):
    batchName = annot_df[annot_df['Metadata_batch_Plate'] == batchPlateName]['batch'].tolist()[0]
    plateName = annot_df[annot_df['Metadata_batch_Plate'] == batchPlateName]['Metadata_Plate'].tolist()[0]
    
    start_time = time.time()
    fileName = os.path.join(rootPath, "backend", batchName, plateName, f"{plateName}.sqlite")

    if not os.path.exists(fileName):
        print(fileName + ' does not exist!')
        return

    selected_fs = ['_Intensity_IntegratedIntensity', '_Intensity_MeanIntensity',
                   '_Intensity_StdIntensity', '_Intensity_MaxIntensity',
                   '_Intensity_UpperQuartileIntensity']

    compartments = ['Cells', 'Nuclei']
    selected_features = [f + '_' + ch for ch in channels_used for f in selected_fs]

    sql_file = "sqlite:///" + fileName
    engine = create_engine(sql_file)
    conn = engine.connect()

    plateDf_list = []
    for comp in compartments:
        selected_features_comp = [comp + f for f in selected_features]
        selected_feature = ' ,'.join(selected_features_comp)
        query_cols = "TableNumber, ImageNumber, ObjectNumber, " + selected_feature
        compartment_query = "select {} from {}".format(query_cols, comp)
        plateDf_list.append(pd.read_sql(sql=compartment_query, con=conn))

    plateDf = reduce(lambda left, right: pd.merge(left, right, on=["TableNumber", "ImageNumber", "ObjectNumber"]), plateDf_list)

    img_query = "select * from {}".format("Image")
    plateImageDf = pd.read_sql(sql=img_query, con=conn)

    plateDfwMeta = pd.merge(plateDf, plateImageDf, on=["TableNumber", "ImageNumber"])

    selected_features_all = [c + f + '_' + ch for c in compartments for ch in channels_used for f in selected_fs]

    transfectionDetectionFeatures = plateDfwMeta.loc[:, selected_features_all + ["Metadata_batch_Plate",
                                                                                 "Metadata_Plate", "Metadata_Well"]]
    os.makedirs(rootPath + '/plate_raw_intensity_features/' + batchName, exist_ok=True)
    transfectionDetectionFeatures.to_pickle(rootPath + '/plate_raw_intensity_features/' + batchName +
                                           '/df_intensityFeatures_' + plateName)
    print("Reading sql, saving int df --- %s minutes ---" % ((time.time() - start_time) / 60))


rootDir = sys.argv[1]
batch_12 = sys.argv[2]
channels_used = ['Protein']

annot_df = pd.read_csv(os.path.join(rootDir, "metadata", "reprocessed", f"{batch_12}.csv"))

listOfBatchPlates = annot_df.Metadata_batch_Plate.unique().tolist()
for bp in listOfBatchPlates:
    saveRawIntensityFeatures(bp, annot_df, rootDir, channels_used)
