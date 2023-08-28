# Pervasive mislocalization of pathogenic coding variants underlying human diseases


 [Link to Paper](will be added after submision to bioarchive)
 
 
 
 # Table of Contents

- [Abstract](#toc-abstract)
- [Dataset](#toc-dataset)
  - [Summary](#toc-summary)
  - [Download](#toc-download)
    - [Downloading the data residing on s3 cellpainting-gallery bucket](#toc-extracting)
  - [Folder Structure](#toc-folder-structure)
  - [Descriptions of each data level](#toc-data-descriptions)
  
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
                

| Batch   | Description  | cell painting channels | protein marker channels | Size (max projected images/ single cell sqlite files) |
| ------- | ------------ | ------------------- | ------------------------------------- | ------------------------------------- |
| PILOT_1     | initial WT/MT screen | `Mito`,`ER`,`DNA`                | `Protein`                                          | 1.04 TB /  74.12 GB | 
| Cancer_Mutations_Screen     | follow up WT/MT screen | `Mito`,`ER`,`DNA`                 | `Protein`                       | 144.5 GB /  14.61 GB          |
| Common_Variants | follow up WT/MT screen | `Mito`,`ER`,`DNA`                 | `Protein`                                   | 56.08 GB /  4.6 GB |
| Kinase_Plates   | follow up WT/MT screen | `Mito`,`ER`,`DNA`                 | `Protein`                                   | 84.06 GB /  8.26 GB
| Replicates_Original_Screen   | replicate of intial WT/MT screen | `Mito`,`ER`,`DNA`                 | `Protein`            | 376.69 GB / 32.88 GB          |
| 2021_05_21_QualityControlPathwayArrayedScreen   | compound screen | `ER`,`DNA`                 | `Protein`,`DsRed`         |  699.91 GB /
| 2022_01_12_Batch1   | compound screen | `DNA`,`Lysosomes`                 | `Protein`,`DsRed`                              |  223.68 GB /
| 2022_01_12_Batch2   | compound screen | `DNA`                 | `Protein`,`DsRed`                                          |  83.97 GB /



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
        │   │   ├── images_unprojected
        │   │   └── images
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
            │   ├── reprocessed
            │   └── raw
            └── profiles
                 ├── singlecell_profiles (Transfected per well single cells)
                 ├── enrichment_profiles
                 └── population_profiles (Average of transfected and untransfected per well)


```

General structure of the images, analysis, backend and load_data_csv for all projects in cellpainting-gallery can be found in the below links:
    - [images folder structure](https://github.com/broadinstitute/cellpainting-gallery/blob/main/folder_structure.md#images-folder-structure)
    - [analysis folder structure](https://github.com/broadinstitute/cellpainting-gallery/blob/main/folder_structure.md#analysis-folder-structure)
    - [backend folder structure](https://github.com/broadinstitute/cellpainting-gallery/blob/main/folder_structure.md#backend-folder-structure)
    - [load_data_csv folder structure](https://github.com/broadinstitute/cellpainting-gallery/blob/main/folder_structure.md#load_data_csv-folder-structure)


## <a id="toc-data-descriptions"></a>Data Level descriptions
### Images
- [General images folder structure for all projects in cellpainting-gallery](https://github.com/broadinstitute/cellpainting-gallery/blob/main/folder_structure.md#images-folder-structure)
- Image file name pattern (`images_unprojected` subfolder)
   - File names pattern:  r(n)c(n)f(n)p(n)-c(n)sk1fk1fl1.tiff
         -  (where (n) is a number describing each variable (the letters).)

     - r = row
     - c = column
     - f = field
     - p = position in the z-stack
     - c = channel
    
- Image file name pattern (`images` subfolder)
  - same as above except that p(n) is always p01
  - These images are max projected images and have been used for inputs to CellProfiler software for feature extraction

### CellProfiler software outputs
   - pipelines can be found at load_data_csv folder which follows [the standard load_data_csv folder structure](https://github.com/broadinstitute/cellpainting-gallery/blob/main/folder_structure.md#load_data_csv-folder-structure)
     
  - CellProfiler generated single-cell profiles and cell outlines
    - This follows standard (link) CP outputs
      - General structure of analysis and backend for all projects in cellpainting-gallery can be found in the below links:
        - [analysis folder structure](https://github.com/broadinstitute/cellpainting-gallery/blob/main/folder_structure.md#analysis-folder-structure)
        - [backend folder structure](https://github.com/broadinstitute/cellpainting-gallery/blob/main/folder_structure.md#backend-folder-structure)

### Preprocessed perturbation level/single cell level profiles
  - Population profiles `enrichment_profiles`
    - population profiles are formed by first detecting the transferected single cells per well and then aggregation of the transfected single cells per well 
  - Enrichment profiles
    - Enrichment profiles are basically the histogram of cell counts per clusters defined in an experiment level
   
```
cellpainting-gallery
└── cpg0026-lacoste_haghighi-rare-diseases
    └── broad
        ├── images
        └── workspace
            └── profiles
                 ├── singlecell_profiles (Transfected per well single cells)
                 ├── enrichment_profiles
                 └── population_profiles (Average of transfected and untransfected per well)

```


      
## <a id="toc-metadata"></a>Metadata
- each batch has a raw annotation file provided by the wet lab and a reprocessed and standardized version that is used as the input for the analysis

```
cellpainting-gallery
└── cpg0026-lacoste_haghighi-rare-diseases
    └── broad
        └── workspace
            └── metadata
               ├── reprocessed
               └── raw
```



| Column  | Description  | Batches which have this column |
| ------- | ------------ | ------------------- |
| Metadata_Batch| one of eight batches | all                |
| Metadata_Plate| plate key | all                |
| Metadata_Well | well key | all                |
| Metadata_Sample   | WT+MT string | all                |
| Variant   | WT+MT string | all                |
| Gene   | WT string | all                |
| MT   | MT string | all                |
| Treatment   | chemical perturbation added to genetic perturbation  | 2021_05_21_QualityControlPathwayArrayedScreen, 2022_01_12_Batch1,  2022_01_12_Batch2              |
| Metadata_Sample_Unique   | 978 | 1677                |
| Metadata_Location    | primary location of protein by visual annotation| PILOT_1                |

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
- The uniformity of single treatment profiles across distinct experimental batches serves as an indicator of data quality. We assess this uniformity through the following procedure:

   - Replicate-level profiles per plate undergo standardization.
   - The Pearson correlation coefficient is computed for each pair of replicate-level profiles corresponding to the same perturbation.
   - The distribution of these coefficients for each dataset and modality is depicted in Supplementary Figure 1 as red curves.
   - The blue curves, juxtaposed to the red curves, represent the null distribution, which displays the correlation coefficient between profile pairs    originating from distinct perturbations.
   - The non-zero dotted vertical line on the right signifies the 90th percentile of the null distribution.
   - Perturbations exhibiting an average replicate correlation exceeding the 90th percentile of the null distribution are deemed high-quality data points, suitable for subsequent analyses.


## <a id="toc-prot-loc"></a>Protein Localization
- `Manual annotations`: for the follwing batches of data, we have per-well annotation captured by biologist's visual inspection of data and ...




