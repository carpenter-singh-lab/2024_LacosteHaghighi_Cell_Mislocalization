# Pervasive mislocalization of pathogenic coding variants underlying human diseases


 [Link to Paper](will be added after submision to bioarchive)
 
 
 
 # Table of Contents

- [Abstract](#toc-abstract)
- [Dataset](#toc-dataset)
  - [Summary](#toc-summary)
  - [Download](#toc-download)
    - [Downloading the data residing on s3 cellpainting-gallery bucket](#toc-extracting)
  - [Folder Structure](#toc-folder-structure)
  - [Data Level descriptions](#toc-data-descriptions)
  
- [Analysis](#toc-analysis)
  - [Transfection Detection](#toc-trans-dec)
    - Filter out untransfected (unperturbed) single cells as the first analysis step
  - [Perturbation Level Profile](#toc-pert-prof)
    - processing CellProfiler single cell outputs to generate per-well level profiles for each perturbation
  - [Techinal replicate reproducibility](#toc-tech-rep)
    - as a measure of profile quality
  - [Protein Localization](#toc-prot-loc)
    - Supervised population level classification of protein localization 
  - [WT/MT impact-score](#toc-prot-loc)
    - Supervised population level classification of protein localization 
  - [Variant/Variant+Treatment reversal-score](#toc-prot-loc)
    - Supervised population level classification of protein localization 
  - [Protein Localization](#toc-prot-loc)
    - Supervised population level classification of protein localization 

  
# <a id="toc-abstract"></a>Abstract

TBA

# <a id="toc-dataset"></a>Dataset

## <a id="toc-summary"></a>Summary
This dataset contains eight batches of data each having various properties.
                

| Batch   | Description  | cell painting channels | protein marker channels |
| ------- | ------------ | ------------------- | ------------------------------------- |
| PILOT_1     | initial WT/MT screen | `Mito`,`ER`,`DNA`                | `Protein`                                   |
| Cancer_Mutations_Screen     | follow up WT/MT screen | `Mito`,`ER`,`DNA`                 | `Protein`                                   |
| Common_Variants | follow up WT/MT screen | `Mito`,`ER`,`DNA`                 | `Protein`                                   |
| Kinase_Plates   | follow up WT/MT screen | `Mito`,`ER`,`DNA`                 | `Protein`                                   |
| Replicates_Original_Screen   | replicate of intial WT/MT screen | `Mito`,`ER`,`DNA`                 | `Protein`                                   |
| 2021_05_21_QualityControlPathwayArrayedScreen   | compound screen | `ER`,`DNA`                 | `Protein`,`DsRed`                                 |
| 2022_01_12_Batch1   | compound screen | `DNA`,`Lysosomes`                 | `Protein`,`DsRed`                                      |
| 2022_01_12_Batch2   | compound screen | `DNA`                 | `Protein`,`DsRed`                                      |



## <a id="toc-download"></a>Download

They can be downloaded at no cost and no need for registration of any sort, using the command:

```bash
aws s3 sync \
  --no-sign-request \
  s3://cellpainting-gallery/cpg0026-lacoste_haghighi-rare-diseases/broad/ .
```

- AWS CLI installation instructions can be found [here](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html).

## <a id="toc-folder-structure"></a>Folder Structure

The parent structure of the dataset is as follows.

```
cellpainting-gallery
└── cpg0026-lacoste_haghighi-rare-diseases
    └── broad
        ├── images
        │   ├── PILOT_1
        │   │   ├── illum
        │   │   ├── unprojected_images
        │   │   └── images [structure](https://github.com/broadinstitute/cellpainting-gallery/blob/main/folder_structure.md#images-folder-structure)
        │   ├── Cancer_Mutations_Screen 
        │   ├── Common_Variants
        │   ├── Kinase_Plates
        │   ├── Replicates_Original_Screen
        │   ├── 2021_05_21_QualityControlPathwayArrayedScreen 
        │   ├── 2022_01_12_Batch1     
        │   └── 2022_01_12_Batch2
        └── workspace
            ├── analysis
            ├── backend
            ├── load_data_csv
            ├── metadata
            └── profiles
```

- [images folder structure](https://github.com/broadinstitute/cellpainting-gallery/blob/main/folder_structure.md#images-folder-structure)
- [analysis folder structure](https://github.com/broadinstitute/cellpainting-gallery/blob/main/folder_structure.md#analysis-folder-structure)
- [backend folder structure](https://github.com/broadinstitute/cellpainting-gallery/blob/main/folder_structure.md#backend-folder-structure)
- [load_data_csv folder structure](https://github.com/broadinstitute/cellpainting-gallery/blob/main/folder_structure.md#load_data_csv-folder-structure)


## <a id="toc-data-descriptions"></a>Data Level descriptions
### Images
- Images (unprojected)
   - File names pattern:  r(n)c(n)f(n)p(n)-c(n)sk1fk1fl1.tiff
     -  where (n) is a number describing each variable (the letters). 

     - r = row
     - c = column
     - f = field
     - p = position in the z-stack
     - c = channel
     
### CellProfiler generated single-cell profiles and cell outlines
  - This follows standard (link) CP outputs 

### Preprocessed perturbation level profiles
  - Population profiles
  - Enrichment profiles

## <a id="toc-metadata"></a>Metadata


| Column  | Description  | Batches which have this column | CP<br/>`normalized_variable_selected` |
| ------- | ------------ | ------------------- | ------------------------------------- |
| Metadata_Batch     | 977 | 1565                | 727                                   |
| Metadata_Plate     | 977 | 1565                | 727                                   |
| Metadata_Well | 977 | 1570                | 601                                   |
| Metadata_Sample   | 978 | 1569                | 291                                   |
| Variant   | 978 | 1569                | 291                                   |
| Gene   | 978 | 1569                | 291                                   |
| MT   | 978 | 1569                | 291                                   |
| Treatment   | 978 | 1569                | 291                                   |
| Metadata_Sample_Unique   | 978 | 1677                | 63                                    |
| Metadata_Location    | 978 | 1670                | 119                                   |

Metadata_Sample
Index(['Metadata_Plate', 'Metadata_Well', 'Metadata_Sample',
       'Metadata_Location', 'Metadata_Notes', 'Metadata_Efficiency', 'rep',
       'Seq ok?', 'Loc 1', 'Loc 2', 'Loc 3', 'Loc 4', 'Loc 5', 'batch', 'DT',
       'Gene', 'Mutation', 'Class', 'Loc 7', 'Loc 6', 'Metadata_transfection',
       'Metadata_batch_Plate', 'Metadata_Sample_Unique', 'MT', 'Variant']
       
       
'Metadata_Plate', 'Metadata_Well', 'Sample name', 'Treatment',
       '3xF tagged construct expressed?', 'In redo plates?', 'batch',
       'control', 'Variant', 'Gene', 'MT', 'Metadata_Sample_Unique',
       'Metadata_batch_Plate', 'rep'
       
       
'Metadata_Plate', 'Metadata_Well', 'Treatment', 'Variant', 'Category',
       'Different antibody used in 647?', 'Gene', 'MT',
       'Metadata_Sample_Unique', 'batch', 'Metadata_Plate2',
       'Metadata_batch_Plate', 'rep', 'control'   


# <a id="toc-analysis"></a>Analysis

## <a id="toc-trans-dec"></a>Transfection Detection
- In general, transfection detection is performed by single cell single feature thresholding method for most of the batches. The raw per-plate single-cell feature values are subjected to truncation at the 0.999th percentile of their respective distributions within each plate.
Subsequently, these clipped per-plate values are normalized to a range of 0 to 1 by employing the per-plate min-max scaling technique. For batch `2022_01_12_Batch1` and `2022_01_12_Batch2`, segmentation is done on `DsRed` channel in image analysis step and any cell that was detected based on that channels was assumed to be transfected. Therefore, CellProfiler outputs for the single cell profiles are just for transfected cells for this batch of data. Single Cell feature used for transfection detection for the remaining batches of data was `Cells_Intensity_IntegratedIntensity_Protein` except for `2021_05_21_QualityControlPathwayArrayedScreen` batch in which we used `Cells_Intensity_UpperQuartileIntensity_DsRed`.

   * **Transfected/Untransfected Profiles:**

## <a id="toc-pert-prof"></a>Perturbation Level Profile

- **Population (Average) Profiles:**
  - Population profiles per-well are form by averaging the single-cell transfected cells (labeled in a prior step) profiles
- **Subpopulation (Enrichment) Profiles:**
  - For extracting subpopulation or enrichement profiles, first we perform a kmeans (k=20) clustering model to a subsampled (1000 cell per plate) population of transfected single cells from all plates. Subpopulation or enrichement profiles for each well are then formed by per-cluster abundance of single cells in that well.


## <a id="toc-tech-rep"></a>Techinal replicate reproducibility



## <a id="toc-prot-loc"></a>Protein Localization
- `Manual annotations`: for the follwing batches of data, we have per-well annotation derived by biologist's visual inspection of data and 




