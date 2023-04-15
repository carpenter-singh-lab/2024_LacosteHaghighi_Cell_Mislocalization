# Pervasive mislocalization of pathogenic coding variants underlying human diseases


 [Link to Paper](will be added after submision to bioarchive)
 
 
 
 # Table of Contents

- [Abstract](#toc-abstract)
- [Dataset](#toc-dataset)
  - [`Download`](#toc-download)
    - [Downloading the data residing on s3 cellpainting-gallery bucket](#toc-extracting)
  - [`Folder Structure`](#toc-folder-structure)
  - [`Data Level descriptions`](#toc-data-descriptions)
  
- [Analysis](#toc-analysis)
  - [`Techinal replicate reproducibility`](#toc-tech-rep)
    - as a measure of profile quality
  - [`Transfection Detection`](#toc-trans-dec)
    - TBA 

  - [`Protein Localization`](#toc-prot-loc)
    - TBA 
 
  
# <a id="toc-abstract"></a>Abstract

TBA

# <a id="toc-dataset"></a>Dataset
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

# <a id="toc-analysis"></a>Analysis
## <a id="toc-tech-rep"></a>`Techinal replicate reproducibility`

## <a id="toc-trans-dec"></a>`Transfection Detection`

## <a id="toc-prot-loc"></a>`Protein Localization`
