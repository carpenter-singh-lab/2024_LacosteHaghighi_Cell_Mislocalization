from sklearn.metrics.pairwise import cosine_similarity
from sklearn.manifold import MDS
import matplotlib.pyplot as plt
import os
import pandas as pd
from .meanProfileAnalysis import (
    generateFeatureGroupsMap_medianZscore,
    featureImportanceEmbeded,
)


def generate_feature_importance_figures(df_scaled_annot, cpFeats_A, save_dir, channels):
    list_of_plates = df_scaled_annot.Metadata_Plate.unique().tolist()

    incompPairs = []

    for p in list_of_plates:
        df_scaled_annot_plateX = df_scaled_annot[
            df_scaled_annot["Metadata_Plate"] == p
        ].reset_index(drop=True)

        df_scaled_trans = df_scaled_annot_plateX[
            df_scaled_annot_plateX[cpFeats_A].isnull().sum(axis=1) == 0
        ][cpFeats_A].T
        dist = 1 - cosine_similarity(df_scaled_trans)
        mds = MDS(n_components=2, dissimilarity="precomputed")
        pos = mds.fit_transform(dist)

        embeddingDf = pd.DataFrame(index=range(3), columns=cpFeats_A)
        embeddingDf.loc[0:1, :] = pos.T

        plt.ioff()
        dfGM = df_scaled_annot_plateX[
            df_scaled_annot_plateX[cpFeats_A].isnull().sum(axis=1) == 0
        ].copy()
        genes = dfGM.Gene.unique().tolist()

        for g in range(len(genes)):
            df_g = dfGM[dfGM["Gene"] == genes[g]]
            wt_mts = df_g.Metadata_Sample.unique().tolist()
            wt_mts.sort(key=len)

            if len(wt_mts) > 1 and len(wt_mts[0].split(" ")) < 2:
                for i in range(1, len(wt_mts)):
                    cat1a = dfGM[dfGM["Metadata_Sample"] == wt_mts[0]].reset_index(
                        drop=True
                    )
                    cat2a = dfGM[dfGM["Metadata_Sample"] == wt_mts[i]].reset_index(
                        drop=True
                    )

                    pairName = wt_mts[i].replace(" ", "_")
                    save_path = os.path.join(save_dir, "WT-MUT-" + pairName, p)
                    os.makedirs(save_path, exist_ok=True)
                    #                     save_path = save_dir + '/WT-MUT-' + pairName + "/" + p
                    #                     os.system("mkdir -p " + save_dir)
                    generateFeatureGroupsMap_medianZscore(
                        cat2a, cat1a, save_path, channels
                    )
                    featureImportanceEmbeded(embeddingDf, cat2a, cat1a, save_path)
            else:
                #                 metaa = df_g.Metadata_Plate.unique().tolist()[0] + '-' + df_g.Metadata_Well.tolist()[0]
                incompPairs.append(wt_mts)

    return incompPairs
