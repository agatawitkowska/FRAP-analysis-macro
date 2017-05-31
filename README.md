FRAP analysis macro v1.0.0
http://dx.doi.org/10.5281/zenodo.376619

This package contains:
- FRAP_analysis_Fiji.ijm
- FRAP_analysis_Octave.m

This macro was developed by Agata Witkowska for the study
Witkowska A. & Jahn R., Biophysical Journal (2017)
"Rapid SNARE-mediated fusion of liposomes and chromaffin granules with giant unilamellar vesicles"
http://dx.doi.org/10.1016/j.bpj.2017.03.010

These macros were developed for the analysis of FRAP data with circular bleach spot generated with ZEN 2010 software (Zeiss) and is optimized to work with a native file format for this software ".czi".
Macros are intended to work together consecutively - first the user starts FRAP_analysis_Fiji.ijm within Fiji program and the file with code for GNU Octave that will be executed after completing this macro. The ".m" file should be located in user's ImageJ directory. For full execution on Windows GNU Octave has to be added to the command line.
Data analysis groups multiple FRAP traces into groups with the same bleach parameters and analyses groups of data instead of individual files.
Required GNU Octave packages: optim, io.
