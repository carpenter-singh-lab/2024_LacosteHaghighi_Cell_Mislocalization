"""
@author: mhaghigh
"""
import pandas as pd


def impact_score_wt_mt(
    df_rep_level_scaled, repCorr_df_avg, featColNames_ls, wt_mt_cols
):
    """
    This function calculates wt mt impact scores for treatment level profiles

    Inputs:
    df_rep_level_scaled   (pd df): input dataframe contains metadata and features
    repCorr_df_avg  (pd df): replicate correlation of samples
    featColNames_ls (list): [cpFeats_P,cpFeats_NP]

    Output:
    impact_scores_trt_profs (pd df):  with columns for score between wt and mt (cc_p,cc_np) and
                                      replicate correlations for mt (RepCor_p,RepCor_np)
                                      and its wt (wt_RepCor_p,wt_RepCor_np)
    """

    impact_scores_df_chs_ls = []
    for f, ch in zip(featColNames_ls, ["p", "np"]):
        impact_corr_mat = (
            df_rep_level_scaled.groupby([wt_mt_cols[0], wt_mt_cols[1]])
            .mean()[f]
            .T.corr()
        )

        genes_variant_size = impact_corr_mat.groupby(wt_mt_cols[0]).size().reset_index()
        genes_with_variant = list(
            set(
                genes_variant_size.loc[
                    genes_variant_size[0] > 1, wt_mt_cols[0]
                ].tolist()
            )
            & set(df_rep_level_scaled[wt_mt_cols[1]].unique().tolist())
        )
        impact_scores_df_ls = []
        for g in genes_with_variant:
            per_gene_df = impact_corr_mat.loc[g][g][g].reset_index()
            #     per_gene_df['Gene']=per_gene_df.columns[1]
            per_gene_df[wt_mt_cols[0]] = g
            per_gene_df["wt_RepCor_" + ch] = repCorr_df_avg.loc[
                repCorr_df_avg[wt_mt_cols[1]] == g, "RepCor_" + ch
            ].values[0]

            impact_scores_df_ls.append(per_gene_df.rename(columns={g: "cc_" + ch}))

        impact_scores_df = pd.concat(impact_scores_df_ls, ignore_index=True)
        impact_scores_df_chs_ls.append(impact_scores_df)

    impact_scores_df_chs = (
        pd.concat(impact_scores_df_chs_ls, axis=1).T.drop_duplicates().T
    )
    impact_scores_trt_profs = impact_scores_df_chs[
        impact_scores_df_chs[wt_mt_cols[1]] != impact_scores_df_chs[wt_mt_cols[0]]
    ].reset_index(drop=True)

    impact_scores_trt_profs = pd.merge(
        impact_scores_trt_profs, repCorr_df_avg, how="inner", on=[wt_mt_cols[1]]
    )
    impact_scores_trt_profs["IS_p"] = impact_scores_trt_profs["cc_p"].apply(
        lambda x: (1 - x) / 2
    )
    impact_scores_trt_profs["IS_np"] = impact_scores_trt_profs["cc_np"].apply(
        lambda x: (1 - x) / 2
    )
    return impact_scores_trt_profs


def impact_score_wt_mt_perplate(
    df_rep_level_scaled, repCorr_df_avg, featColNames_ls, wt_mt_cols
):
    """
    This function calculates wt mt impact scores for treatment level profiles

    Inputs:
    df_rep_level_scaled   (pd df): input dataframe contains metadata and features
    repCorr_df_avg  (pd df): replicate correlation of samples
    featColNames_ls (list): [cpFeats_P,cpFeats_NP]

    Output:
    impact_scores_perplate (pd df):  with columns for score between wt and mt (cc_p,cc_np) and
                                      replicate correlations for mt (RepCor_p,RepCor_np)
                                      and its wt (wt_RepCor_p,wt_RepCor_np)
    """
    impact_scores_df_chs_ls = []

    plates = df_rep_level_scaled.Metadata_batch_Plate.unique()

    for f, ch in zip(featColNames_ls, ["p", "np"]):
        impact_scores_df_ls = []
        for p in plates:
            df_rep_plate = df_rep_level_scaled[
                df_rep_level_scaled["Metadata_batch_Plate"] == p
            ].reset_index(drop=True)
            impact_corr_mat = (
                df_rep_plate.groupby([wt_mt_cols[0], wt_mt_cols[1]]).mean()[f].T.corr()
            )

            genes_variant_size = (
                impact_corr_mat.groupby(wt_mt_cols[0]).size().reset_index()
            )
            genes_with_variant = list(
                set(genes_variant_size[genes_variant_size[0] > 1].Gene.tolist())
                & set(df_rep_plate.Metadata_Sample_Unique.unique().tolist())
            )

            for g in genes_with_variant:
                per_gene_df = impact_corr_mat.loc[g][g][g].reset_index()
                #             sfgsg
                #     per_gene_df['Gene']=per_gene_df.columns[1]
                per_gene_df[wt_mt_cols[0]] = g
                per_gene_df["Metadata_batch_Plate"] = p

                impact_scores_df_ls.append(per_gene_df.rename(columns={g: "cc_" + ch}))

        #         fsdfgsfs
        if impact_scores_df_ls:
            impact_scores_df = pd.concat(impact_scores_df_ls, ignore_index=True)
            #             impact_scores_df['Metadata_batch_Plate']=p
            impact_scores_df_chs_ls.append(impact_scores_df)

    impact_scores_df_chs = (
        pd.concat(impact_scores_df_chs_ls, axis=1).T.drop_duplicates().T
    )
    impact_scores_perplate = impact_scores_df_chs[
        impact_scores_df_chs[wt_mt_cols[1]] != impact_scores_df_chs[wt_mt_cols[0]]
    ].reset_index(drop=True)

    impact_scores_perplate = pd.merge(
        impact_scores_perplate, repCorr_df_avg, how="inner", on=[wt_mt_cols[1]]
    )
    impact_scores_perplate["IS_p"] = impact_scores_perplate["cc_p"].apply(
        lambda x: (1 - x) / 2
    )
    impact_scores_perplate["IS_np"] = impact_scores_perplate["cc_np"].apply(
        lambda x: (1 - x) / 2
    )

    return impact_scores_perplate
