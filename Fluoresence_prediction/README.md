# Fluorescence Prediction

This folder contains an example notebook for predicting fluorescence.

For running the notebook from top to bottom, raw data (zipped) is needed.
Intermediate files (`promoters_TRAIN_FACS.csv, promoters_VAL_FACS.csv, train_dat.csv, val_dat.csv`) are provided as well and loaded in the notebook.

The notebook (`Fluorescence_Prediction.Rmd`) can be ran using RStudio or the command line.
e.g. filling in the FILENAME
```
R -e rmarkdown::render"('Fluorescence_Prediction.Rmd',output_file='FILENAME.html')"
```

The output of this notebook is provided in `Fluorescence_Prediction.md` and `Fluorescence_Prediction.html`.


`Fluorescence_Prediction.md` is also rendered directly in Github as a notebook with all output inside.
