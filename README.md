# Pervasive mislocalization of pathogenic coding variants underlying human diseases


 [Link to Paper](will be added after submision to bioarchive)
 
 
 
 # Table of Contents

- [Abstract](#toc-abstract)
- [Dataset](#toc-dataset)
  - [Download](#toc-download)
    - [Downloading the data residing on s3 cellpainting-gallery bucket](#toc-extracting)
  - [Folder Structure](#toc-folder-structure)
  - [Data Level descriptions](#toc-data-descriptions)
  
- [Analysis](#toc-analysis)
  - [Techinal replicate reproducibility](#toc-tech-rep)
    - as a measure of profile quality
  - [Transfection Detection](#toc-trans-dec)
    - TBA 

  - [Protein Localization](#toc-prot-loc)
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


## <a id="toc-tech-rep"></a>Techinal replicate reproducibility



## <a id="toc-prot-loc"></a>Protein Localization





