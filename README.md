# rhiziis

High throughput in situ imaging reveals complex ecological behaviour of giant marine mixotrophic protists

Thelma Panaïotis PhD Thesis

Laboratoire d'Océanographie de Villefranche (UMR 7093)

## Organisation of repository

### Folders

-   `data` where data lives

-   `figures` figures for the paper

-   `lib` useful stuff

-   `model_weights` CNN weights for classification of Rhizaria (generated from training and used for prediction)

-   `plot` as the name suggests, plots

-   `figures` figures for the paper


### Scripts

Scripts order is self-explanatory.

#### Data generation

-   `00.download_data.py` Download data from the Rhizaria EcoTaxa project (6133).

-   `01.download_images.py` Download images from the Rhizaria EcoTaxa project (6133).

-   `02.prepare_learning_set.R` Prepare the learning set from all validated images to train a CNN classifier.

-   `03.train_cnn.py` Train the CNN.

-   `04.predict_vignettes.py` Predict all vignettes. These predictions are then imported into the Ecotaxa project 6133.

-   `05.generate_calibration_set.R` Prepare a calibration set to assess prediction quality. This dataset will be imported to EcoTaxa for manual inspection.

-   `06.download_calibration_data.py` Download data from the calibration set (Ecotaxa project 6521).

-   `07.calibrate_cnn.Rmd` Perform calibration, i.e. set a confidence threshold for each taxa.

-   `08.process_kp_annotations.R` Read, reformat and save keypoint annotations.

#### Data preparation

Start here to reproduce the analyses with the data shared on XXX.

-   `09.threshold_predictions.Rmd` Apply confidence thresholds to predictions, but keep validated objects.

-   `10.concentration_interpolation.Rmd` Compute plankton concentrations, interpolate plankton and environmental data.

#### Data analysis

-   `11.analyses_transects.R` Perform all analyses and plots on transect data.

-   `12.analyses_orientation.R` Analyse orientation of vacuoles in solitary Collodaria, phaeodium in Phaeodaria, and Arthracanthida.

### Data

#### External data

-   Plankton data `ecotaxa_export_all_rhizaria_cc.tsv` is shared on XXX.

-   Environmental data `isiis_images_env.parquet` is shared on XXX.

#### Local data

-   `data_git/taxa_threshold_rhizaria.csv` Score thresholds for each taxa.

-   `data_git/taxonomy_match.csv` Table of taxonomy matches to rename objects.

-   `data_git/08.keypoint_annotations.csv` Keypoint tool annotations.
