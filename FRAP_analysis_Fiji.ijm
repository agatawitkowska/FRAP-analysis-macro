// FRAP_analysis_Fiji.ijm v1.0.0
// This macro was developed by Agata Witkowska for the study
// Witkowska A. & Jahn R., Biophysical Journal (2017)
// "Rapid SNARE-mediated fusion of liposomes and chromaffin granules with giant unilamellar vesicles"
// http://dx.doi.org/10.1016/j.bpj.2017.03.010

// This Fiji macro was developed for the analysis of FRAP data with circular bleach spot generated with ZEN 2010 software (Zeiss) and is optimized to work with a native file format for this software ".czi"
// The file with code for GNU Octave that will be executed after completing this macro should be located in your ImageJ directory.
// For full execution of this macro you need GNU Octave to be added to your command line


files = newArray(0);
filesOctave = newArray(0)
dir = getDirectory("Select the directory"); // here you select a directory where your FRAP data files are located
list = getFileList(dir);
for (i = 0; i < list.length; i++) {
  if (endsWith(list[i], ".czi") && indexOf(list[i], "FRAP") >= 0) { // the files that will be analysed should contain ".czi" extension and word "FRAP" in their name
    files = Array.concat(files, list[i]);
  }
}

setBatchMode(true);

for (f = 0; f < files.length; f++) {

  run("Bio-Formats", "open='" + dir + files[f] + "' color_mode=Default view=Hyperstack stack_order=XYCZT"); // opening files

  image = replace(getTitle(), ".czi", "");
  height = getHeight();
  getPixelSize(unit, pixelWidth, pixelHeight);
  pxwidth = 1 * pixelWidth;

  ROIx = 1 * getInfo("Layer|Circle|Geometry|CenterX #1"); // gathering info about the bleach spot
  ROIy = 1 * getInfo("Layer|Circle|Geometry|CenterY #1");
  ROIr = 1 * getInfo("Layer|Circle|Geometry|Radius #1");

  x = ROIx - ROIr;
  y = ROIy - ROIr;

  makeOval(x, y, 2 * ROIr, 2 * ROIr); // reproducing bleach ROI from ZEN
  run("ROI Manager...");
  roiManager("Add");
  roiManager('select', 0);
  roiManager("Rename", "FRAP");
  run("Set Measurements...", "mean display redirect=None decimal=5");

  makeRectangle(0, height - 6, 6, 6); // generation of a background ROI - here you can change its size and localization
  roiManager("Add");
  roiManager('select', 1);
  roiManager("Rename", "background");

  roiManager('select', newArray(0, 1));
  roiManager("Multi Measure");
  selectWindow("Results");

  lastRow = nResults;
  setResult("Label", lastRow, "bleach frame");
  setResult("Mean(FRAP)", lastRow, "frame time");
  setResult("Mean(background)", lastRow, "pixel size");
  setResult("new", lastRow, "FRAP radius");

  lastRow = nResults;
  setResult("Label", lastRow, getInfo("Experiment|AcquisitionBlock|MultiTrackSetup|TrackSetup|BleachSetup|BleachParameterSet|StartNumber #1"));
  setResult("Mean(FRAP)", lastRow, getInfo("Information|Image|Channel|LaserScanInfo|FrameTime #1"));
  setResult("Mean(background)", lastRow, pixelWidth);
  setResult("new", lastRow, ROIr * pxwidth);

  outfile = dir + image + "_Results.txt";
  saveAs("Results", outfile); // saving results file for each image separately with table including FRAP and background values as well as information about the bleach frame, frame time etc. in the last row of the file
  filesOctave = Array.concat(filesOctave, outfile);

  roiManager("Deselect");
  roiManager("Delete");
  selectWindow("Results");
  run("Close");
  run("Close All");
}

setBatchMode(false);

imagejdir = getDirectory("imagej"); // here your Octave script file should be located

Array.print(filesOctave);
selectWindow("Log");
saveAs("Text", imagejdir + "Log.txt"); // printing log table with results file list for Octave processing
run("Close");

exec("cmd", "/c", "octave --persist FRAP_analysis_Octave.m"); // execution of the Octave script - make sure file name and localization are correct