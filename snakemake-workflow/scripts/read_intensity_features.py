import os
import sys
import pandas as pd
import seaborn as sns
from sklearn import preprocessing as sp
import matplotlib.pyplot as plt

def readNormalizeIntensityFeatures(annot_df,rootDir,channels_used):
    
    listOfBatchPlates=annot_df.Metadata_batch_Plate.unique().tolist();
    scaler_minmax = sp.MinMaxScaler(feature_range=(0,1))

    df_inten = pd.DataFrame();
    df_inten_scaled_perPlate = pd.DataFrame();

    for bp in listOfBatchPlates: #[0:1]:

        batch=annot_df[annot_df['Metadata_batch_Plate']==bp]['batch'].tolist()[0]
        p=annot_df[annot_df['Metadata_batch_Plate']==bp]['Metadata_Plate'].tolist()[0]   

        fileNameToSave=rootDir+'/plate_raw_intensity_features/'+batch+'/df_intensityFeatures_'+p;
        if os.path.exists(fileNameToSave):
            intFeaturesDf=pd.read_pickle(fileNameToSave, compression='infer');    
            df_inten=df_inten.append(intFeaturesDf, ignore_index=True)  
            df_inten_scaled0 = intFeaturesDf.copy()
            intFeatures=intFeaturesDf.columns[intFeaturesDf.columns.str.contains('|'.join(channels_used))].tolist()
    #         sfsdfdsfds
            for ifi in intFeatures:
                qpi_up=intFeaturesDf[ifi].quantile(0.999)
                qpi_low=intFeaturesDf[ifi].quantile(0.01)
                intFeaturesDf[ifi]=intFeaturesDf[ifi].clip(qpi_low, qpi_up)

        #     dataScaled=scaler0.fit_transform(intFeaturesDf.loc[:,intFeatures])
            df_inten_scaled0[intFeatures]=scaler_minmax.fit_transform(intFeaturesDf.loc[:,intFeatures])
            df_inten_scaled_perPlate =df_inten_scaled_perPlate.append(df_inten_scaled0, ignore_index=True)  

#     print(df_inten.shape)   

    df_inten=pd.merge(df_inten, annot_df, how='inner', on=['Metadata_Plate','Metadata_Well']);
    df_inten_scaled_perPlate=pd.merge(df_inten_scaled_perPlate, annot_df, how='inner', on=['Metadata_Plate','Metadata_Well']);
    
    return df_inten,df_inten_scaled_perPlate



def plot_intensity_features(listOfBatches,df_ls,params):
    from matplotlib import rcParams
    rcParams['patch.force_edgecolor'] = False

    intFeatures=params['intFeature']
    
    fig, axes = plt.subplots(2,len(listOfBatches), figsize=(len(listOfBatches)*5,2*5),sharex=params['sharex_enabled'])
    # fig.suptitle('_'.join(intFeatures[i].split('_')[2:]));

    for bi in range(len(listOfBatches)):
        for si,dfi,sn in zip([0,1],df_ls,['raw','normalized']):
    #     for si,dfi in enumerate([df_inten,df_inten_scaled_perPlate]):
            df_inten_b=dfi[(dfi['batch']==listOfBatches[bi]) & (~dfi[params['hue_column']].isnull())].reset_index(drop=True)
            #### Method 1 -a

            if params['hue_column']:
                qpi_up=df_inten_b.loc[df_inten_b[params['hue_column']]==True,intFeatures].quantile(0.99).round(4)
                qpi_low=df_inten_b.loc[df_inten_b[params['hue_column']]==True,intFeatures].quantile(0.4).round(4)

                print(sn,listOfBatches[bi],qpi_low,qpi_up)        

            if df_inten_b[intFeatures].min()==0:
                df_inten_b[intFeatures]=df_inten_b[intFeatures]+params['zero_offset_4logscale_plot']
                qpi_up=df_inten_b.loc[df_inten_b[params['hue_column']]==True,intFeatures].quantile(0.99)


            if params['hue_column']:
                sns.histplot(data=df_inten_b,x=intFeatures, bins=params['hist_n_bins'],stat="density",hue=params['hue_column'],\
                     element="step",common_norm=False,legend=True,log_scale=params['log_scale_enabled'],ax=axes[si,bi])  
    #             perc95=qpi_up+zero_offset_4logscale_plot
                axes[si,bi].axvline(x=qpi_up,linestyle=':',color="r")
            else:

                sns.histplot(data=df_inten_b,x=intFeatures, bins=params['hist_n_bins'],stat="density",\
                     element="step",common_norm=False,legend=True,log_scale=params['log_scale_enabled'],ax=axes[si,bi])

        axes[0,bi].set_title(listOfBatches[bi]);
    axes[0,0].set_ylabel('raw');
    axes[1,0].set_ylabel('scaled');
    return fig 


# Main script here
rootDir = sys.argv[1]
batch_12 = sys.argv[2]

channels_used = ['Protein']

annot_df = pd.read_csv(os.path.join(rootDir, "metadata", "reprocessed", f"{batch_12}.csv"))

df_inten, df_inten_scaled_perPlate = readNormalizeIntensityFeatures(annot_df, rootDir, channels_used)

intFeatures = ['Cells_Intensity_UpperQuartileIntensity_Protein',
               'Cells_Intensity_MeanIntensity_Protein',
               'Cells_Intensity_MaxIntensity_Protein']

plot_int_params = {'log_scale_enabled': True,
                   'sharex_enabled': False,
                   'hist_n_bins': 1000,
                   'zero_offset_4logscale_plot': 0.00001,
                   'intFeature': intFeatures[1],
                   'hue_column': 'Metadata_transfection'}

listOfBatches = annot_df.batch.unique().tolist()

for inti in intFeatures:
    plot_int_params['intFeature'] = inti
    fig = plot_intensity_features(listOfBatches, [df_inten, df_inten_scaled_perPlate], plot_int_params)

    saveDir = os.path.join(rootDir, "results", "intensity_dists", batch_12)
    os.makedirs(saveDir, exist_ok=True)
    fig.savefig(os.path.join(saveDir, f"{inti}.png"))
