Included is three scripts for the quantification of select marker genes on endothelium in whole slide immunoflourescence images from the human pancreas


processing_script.py includes script using scikit-image to segment endothelial regions using VECAD intensity, and quanitfy the raw and normalized expression of the select marker along with other attributes at the object and patch level. The output is a folder containing csv files for each patch, and a csv file combining each of these into one whole dataframe an individual marker

filtering_script.R includes script for QC filtering and scaling intensity by area. The output is a csv file with the QC filtered quantification for each file and folder (one of which is Diabetic one of which is Non-Diabetic)

plotting_script.py includes script for plotting and t-test of QC fliltered dataframe for comparison
