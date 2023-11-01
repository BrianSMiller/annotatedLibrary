# annotatedLibrary  

This repository contains code used in the creation of the IWC-SORP Blue and Fin Whale Library of annotated sounds. 

The library itself is publicly available and accessible from the Australian Antarctic Data Centre ([Miller et al 2020](http://data.aad.gov.au/metadata/AcousticTrends_BlueFinLibrary)). There is a corresponding publication in Scientific Reports as well, [Miller et al 2021](http://www.nature.com/articles/s41598-020-78995-8).

Rather than being an archive of exemplary recordings, like the [Watkins Marine Mammal Sound Database](https://cis.whoi.edu/science/B/whalesounds/index.cfm), this Annotated Library is intended to contain a sample of sounds that is representative of real-world long-term underwater recordings. Such representative recordings are more appropriate for training and testing detection algorithms than curated 'best of' recordings. 

The code included in this repository is a grab-bag of different functions that were developed during creation of the IWC-SORP Anntotated Library, as well as functions used to analyse these data for the resulting publications. Some of the functionality provided by this code includes:
  * Functions for sub-sampling long-term datasets
  * Functions for reading Raven Selection Tables
  * Functions for reading Koogu Detections (which are essentially a sub-type of Raven Selection Table)
  * Functions for reading PAMGuard detections from a Sqlite database (for the Ishmael Spectrogram Correlation Detector and Ishmael Energy Detector)
  * Functions for estimating the signal-to-noise ratio of annotations and detections

## Installation
1) Download and install dependencies. The code in this repository depends on: 
  * (soundFolder)[https://github.com/BrianSMiller/soundFolder] to streamline accessing of wav files.
  * (mkqlite)[https://mksqlite.sourceforge.net/] for accessing PAMGuard sqlite databases
2) Download the repository to a location on your computer
3) Add the location of the downloaded files to the Matlab path

## References
Miller, B.S., Stafford, K.M., Van Opzeeland, I., Harris, D., Samaran, F., Širović, A., Buchan, S., Findlay, K., Balcazar, N., Nieukirk, S., Leroy, E.C., Aulich, M., Shabangu, F.W., Dziak, R.P., Lee, W., Hong, J., 2020. An annotated library of underwater acoustic recordings for testing and training automated algorithms for detecting Antarctic blue and fin whale sounds. Version 1. Australian Antarctic Data Centre. DOI: 10.26179/5e6056035c01b. http://data.aad.gov.au/metadata/records/AcousticTrends_BlueFinLibrary

Miller, B. S., IWC SORP Acoustic Trends Working Group, Kathleen M. Stafford, Ilse Van Opzeeland, Danielle Harris, Flore Samaran, Ana Širović, Susannah Buchan, Ken Findlay, Balcazar, N., Nieukirk, S., Leroy, E. C., et al. (2021). An open access dataset for developing automated detectors of Antarctic baleen whale sounds and performance evaluation of two commonly used detectors. Scientific Reports 11, 806. DOI: 10.1038/s41598-020-78995-8.
