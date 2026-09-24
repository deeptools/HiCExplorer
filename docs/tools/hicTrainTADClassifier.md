# hicTrainTADClassifier

Trains and tests the classifiers used by hicTADClassifier.

```text
usage: hicTrainTADClassifier --mode {train_new,train_existing,train_test,predict_test} --matrices
                            MATRICES [MATRICES ...] --domain_file DOMAIN_FILE [DOMAIN_FILE ...]
                            --out_file OUT_FILE [--normalization_method {obs_exp,range}]
                            [--resolution RESOLUTION] [--threshold THRESHOLD]
                            [--leniency LENIENCY] [--saved_classifier SAVED_CLASSIFIER]
                            [--unselect_border_cases]
                            [--protein_file PROTEIN_FILE [PROTEIN_FILE ...]] [--threads THREADS]
                            [--chromosomes CHROMOSOMES [CHROMOSOMES ...]]
                            [--concatenate_before_resample] [--estimators_per_step [5-1000]]
                            [--resampling_method {undersample_cluster_centroids,undersample_random,passed_method}]
                            [--alternative_resampling_method ALTERNATIVE_RESAMPLING_METHOD]
                            [--distance [5-30]] [--impute_value IMPUTE_VALUE]
                            [--alternative_classifier ALTERNATIVE_CLASSIFIER] [--use_cleanlab]
                            [--help] [--version]
```

Check out hicTADClassifier for calling TADs using our default classifiers or your trained classifier.

This program can be used to train and test new and existing classifiers for hicTADClassifier. These classifiers can later be run to call boundaries for TADs. There are four modes available: train_new, train_existing, train_test and predict_test. By default, an EasyEnsembleClassifier as described in Liu et al.: "Exploratory Undersampling for Class-Imbalance Learning" will be trained, but you can pass any sklearn classifier that allows for a warm start. You may also vary the resampling method and a range of hyperparameters to fine tune the model. Do mind to set the correct normalization method and resolution for the classifier. The program will check and raise warnings, when resolutions and normalization methods are mixed up. Also, a protein track file in the narrowPeak format with a threshold value may be passed to filter out low quality boundaries.


train_test mode: this is a convenience function, where a single matrix/domains set can be passed to quickly assert the performance of a new classifier. Nothing will be saved from this mode, instead, the classifier will be trained on 80% of the data and tested on the remaining 20%. The output will be a performance report. A quick usage example can be seen here:

train_new mode: this mode allows the training of a new classifier. Note that range of optional arguments, that can be used to fine tune. The resulting classifier will be pickled at the specified out_file. A quick example can be seen here, where we varied the feature distance:

train_existing mode: train the classifier specified in saved_classifier on new data. When not setting the saved_classifier, the preset classifiers will be used as preset. The output will be the classifier trained with additional data and more internal estimators.

predict_test mode: predict using an existing classifier and produce a classification report. The difference in using this over hicTADClassifier is, that this version will predict on a balanced test set. Normally, HiC-Matrices contain a lot more non-boundaries than boundaries, which skews the classification report to the point, where it does not contain usefull information anymore. By passing a domain file produced by another TAD Caller, hicTrainClassifier will build a test set using the boundaries of this domain file and will pick at random as many non-boundaries from the passed matrix. Use this over hicTADClassifier to produce a meaningful output, but not for TAD calling.

This tool runs in Python. Its machine-learning dependencies are optional; install them with `pip install "hicexplorer[tadclassifier]"`.

## Required arguments

| Flag | Meaning |
|---|---|
| `--mode {train_new,train_existing,train_test,predict_test}`, `-mo {train_new,train_existing,train_test,predict_test}` | choice of program |
| `--matrices MATRICES`, `-m MATRICES` | HiC-Matrix file or list of files for input. Only COOLER files are supported! |
| `--domain_file DOMAIN_FILE`, `-d DOMAIN_FILE` | domain file or list of files containing tad boundaries |
| `--out_file OUT_FILE`, `-o OUT_FILE` | output file for either the classification report or the saved classifier |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--normalization_method {obs_exp,range}`, `-n {obs_exp,range}` | set the normalization mode, with which the passed matrices will be normalized |
| `--resolution RESOLUTION`, `-r RESOLUTION` | resolution in bases of the classifier |
| `--threshold THRESHOLD` | threshold for protein quality check |
| `--leniency LENIENCY` | leniency for protein quality check. Widens peaks of protein file by leniency*resolution |
| `--saved_classifier SAVED_CLASSIFIER` | pickled classifier to be trained or used for prediction |
| `--unselect_border_cases` | set whether genes at the border of the matrix up to set distance will not be used for training and testing |
| `--protein_file PROTEIN_FILE` | provide a bed file for TAD quality control |
| `--threads THREADS`, `-t THREADS` | number of threads used |
| `--chromosomes CHROMOSOMES` | Chromosomes to include in the analysis. If not set, all chromosomes are included. |
| `--concatenate_before_resample` | whether features build from matrix list are concatenated and resampled together or resampled separatly per matrix. Not important for random undersampling, but alter for other resampling methods and check if performance increases. |
| `--estimators_per_step [5-1000]` | how many estimators are added in each training step for the classifier (new classifier) |
| `--resampling_method {undersample_cluster_centroids,undersample_random,passed_method}` | the method used to resample the training set(new classifier) |
| `--alternative_resampling_method ALTERNATIVE_RESAMPLING_METHOD` | pass alternative resampling method from imblearn library (new classifier) |
| `--distance [5-30]` | max distance between TADs to be used in calculation (new classifier) |
| `--impute_value IMPUTE_VALUE` | non-numerical float values in matrix will be replaced by this value (new classifier) |
| `--alternative_classifier ALTERNATIVE_CLASSIFIER` | pass custom classifier, needs to implement warm_start (new classifier) |
| `--use_cleanlab` | use Confident Learning with the cleanlab library (new classifier) |
| `--version` | show program's version number and exit |
