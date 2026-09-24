# hicTADClassifier

Predicts TAD boundaries with a supervised classifier.

```text
usage: hicTADClassifier --matrices MATRICES [MATRICES ...] --out_file OUT_FILE [OUT_FILE ...]
                       [--normalization_method {obs_exp,range}]
                       [--saved_classifier SAVED_CLASSIFIER] [--unselect_border_cases]
                       [--threads THREADS] [--chromosomes CHROMOSOMES [CHROMOSOMES ...]]
                       [--help] [--version]
```

Uses Supervised Learning to call TAD boundaries. One or multiple HiC-Matrices can be passed, from which a BED file will be produced containing the predicted boundary positions. By default, a EasyEnsembleClassifier as described in Liu et al.: "Exploratory Undersampling for Class-Imbalance Learning" will be used to call TADs. Internally this classifier relies on Resampling, Boosting and Bagging. Passed matrices will be range normalized by default. Alternatively, obs/exp normalization can be used. Currently, only classifiers for 10kb resolution are implemented. For building own classifiers or tune existing ones, hicTrainClassifier can be used and passed with the saved_classifer argument. A simple usage example can be seen here:

This tool runs in Python and is part of the HiCExplorer 4 package.

## Required arguments

| Flag | Meaning |
|---|---|
| `--matrices MATRICES`, `-m MATRICES` | HiC-Matrix file or list of files for input. Only COOLER files are supported! |
| `--out_file OUT_FILE`, `-o OUT_FILE` | output file path for predictions |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--normalization_method {obs_exp,range}`, `-n {obs_exp,range}` | set the normalization mode, with which the passed matrices will be normalized. If not set, matrices will be range normalized |
| `--saved_classifier SAVED_CLASSIFIER` | Default classifier are available for 10kb, 25kb, 50kb and 100kb resolution. Do not set this parameter to use the default models. Pass a self-trained classifier (from hicTrainTADClassifier) to load a non-default model. |
| `--unselect_border_cases` | set whether genes at the border of the matrices will not be predicted |
| `--threads THREADS`, `-t THREADS` | number of threads used |
| `--chromosomes CHROMOSOMES` | Chromosomes to include in the analysis. If not set, all chromosomes are included. |
| `--version` | show program's version number and exit |
