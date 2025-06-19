# IR-Fingerprint

IR-Fingerprint is a newly proposed method of constructing molecular fingerprints, that integrates both molecular structure of a compound and its spectral data (MIR spectrum). At present, the method has been employed to develop logP prediction models using a variety of algorithms, with Support Vector Regression demonstrating the greatest efficacy. Fingerprints generated in this manner demonstrated acceptable results, although they were inferior to conventional fingerprints such as Morgan or MACCS fingerprints.

The following items are included in this repository:
- Python script analyzing molecular structure and generation of theoretical fingerprint ("TheoreticalFingerprint.py")
- R script analyzing MIR spectrum of a given compound and generating molecular fingerprint ("BaseFingerprint.R") 
- R script constructing QSAR models and performing feature selection ("FeatureSelection.R")
- R script performing clustering based on proposed fingerprint ("Clustering.R")
- Database of used molecules with logP values
- CSV files with generated fingerprints
