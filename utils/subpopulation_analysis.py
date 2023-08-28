import numpy as np
import scipy.spatial
import pandas as pd
import os
from singlecell.process import normalize_funcs
from singlecell.visualize import visualize_n_SingleCell, cluster


def generate_subpopulation_analysis_figures(
    df_rep_level_scaled, feats, root_dir, save_dir, params
):
    list_of_batch_plates = df_rep_level_scaled["Metadata_batch_Plate"].unique().tolist()
    incompPairs = []

    for bp in list_of_batch_plates:
        b = bp.split("-")[0]
        p = bp.split("-")[1]
        df_scaled_annot_plateX = df_rep_level_scaled[
            df_rep_level_scaled["Metadata_batch_Plate"] == bp
        ].reset_index(drop=True)

        #     dfGM = df_scaled_annot_plateX[df_scaled_annot_plateX[cpFeats_A].isnull().sum(axis=1) == 0].copy()
        genes = df_scaled_annot_plateX["Gene"].unique().tolist()

        for g in range(len(genes)):
            print(list_of_batch_plates.index(bp), bp, g, genes[g])
            df_g = df_scaled_annot_plateX[df_scaled_annot_plateX["Gene"] == genes[g]]
            wt_mts = df_g.Metadata_Sample.unique().tolist()
            wt_mts.sort(key=len)

            if len(wt_mts) > 1 and len(wt_mts[0].split(" ")) < 2:
                wells_wt = df_g.loc[
                    df_g["Metadata_Sample"] == wt_mts[0], "Metadata_Well"
                ].values.tolist()
                wells_wt_ls = []
                for wi in range(len(wells_wt)):
                    wt_sc_path = (
                        root_dir
                        + "/singlecell_profiles/"
                        + b
                        + "/"
                        + p
                        + "/"
                        + p
                        + "_"
                        + wells_wt[wi]
                        + ".csv.gz"
                    )
                    wells_wt_ls.append(pd.read_csv(wt_sc_path))
                wells_wt_df = pd.concat(wells_wt_ls, ignore_index=True)
                wells_wt_df["label"] = "Wild-type"

                for i in range(1, len(wt_mts)):
                    wells_mt = df_g.loc[
                        df_g["Metadata_Sample"] == wt_mts[i], "Metadata_Well"
                    ].values.tolist()

                    wells_mt_ls = []
                    for wi in range(len(wells_mt)):
                        mt_sc_path = (
                            root_dir
                            + "/singlecell_profiles/"
                            + b
                            + "/"
                            + p
                            + "/"
                            + p
                            + "_"
                            + wells_mt[wi]
                            + ".csv.gz"
                        )
                        wells_mt_ls.append(pd.read_csv(mt_sc_path))
                    wells_mt_df = pd.concat(wells_mt_ls, ignore_index=True)
                    wells_mt_df["label"] = "Mutant"

                    wtANDmtDf = pd.concat([wells_wt_df, wells_mt_df], ignore_index=True)

                    pairName = wt_mts[i].replace(" ", "_")
                    save_path = os.path.join(save_dir, "WT-MUT-" + pairName, p)

                    os.makedirs(save_path, exist_ok=True)

                    feature_scaler = "Robust"
                    standardized_sc_df = normalize_funcs.standardize_df_columns(
                        wtANDmtDf, feats, feature_scaler
                    )

                    standardized_sc_df = cluster.add_clustering_index_column(
                        standardized_sc_df, feats, "kmeans", params["n_clusters"]
                    )

                    params["results_dir"] = save_path

                    cluster.subpopulation_clustering_visualization(
                        standardized_sc_df, feats, params
                    )

            else:
                incompPairs.append(wt_mts)

    return incompPairs
