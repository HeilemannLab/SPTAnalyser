// @author: Claudia Catapano, Marina Dietz
// Research group Heilemann
// Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

// This macro loads in a directory and the areas of all *.roi files within the directory and all subdirectories will be measured and saved as *.csv file in the same directory as the respective *.roi file

// to use the script, open one acquired *.tif movie first, then run the script in Fiji


run("Set Measurements...", "area redirect=None decimal=9");
run("Set Scale...", "distance=0 known=0 unit=pixel global");

// Funktion, die rekursiv alle .roi-Dateien in einem Verzeichnis und seinen Unterordnern sucht
function getROIFiles(directory) {
    list = getFileList(directory);
    roiFiles = newArray(0);
    
    for (i = 0; i < list.length; i++) {
        if (endsWith(list[i], ".roi")) {
            roiFiles = Array.concat(roiFiles, directory + list[i]);
        } else if (File.isDirectory(directory + list[i])) {
            subDir = directory + list[i] + "/";
            subFiles = getROIFiles(subDir);
            roiFiles = Array.concat(roiFiles, subFiles);
        }
    }
    return roiFiles;
}

// Prompt user to select a directory
input = getDirectory("Select a directory containing ROI files");

// Get all .roi files in the selected directory and its subdirectories
roiList = getROIFiles(input);

// Create an empty ROI Manager
roiManager("reset");
newResults = newArray(0);

// Loop through each .roi file found
for (i = 0; i < roiList.length; i++) {
    path = roiList[i];
    
    // Open the current ROI file
    roiManager("Open", path);
    
    roiManager("Select", roiManager("Count") - 1); // select the last added ROI
    
    // Measure the current ROI
    run("Measure");
    
    // Get the measurements as a string
    measurements = getResult("Label", 0) + "\t" + getResult("Area", 0) + "\t" + getResult("Mean", 0) + "\t" + getResult("StdDev", 0) + "\t" + getResult("Min", 0);
    
    // Get the original folder of the ROI file
    roiFolderPath = File.getParent(path);
    
    // Save measurements to a CSV file in the same folder as the original ROI file
    resultFile = replace(File.getName(path), ".roi", ".csv");
    saveAs("Measurements", roiFolderPath + File.separator + resultFile);
    
    // Close the ROI Manager to clear the selections
    roiManager("Deselect");
    run("Clear Results");
    close("*Results*");
}

waitForUser("Macro Finished!", "All files in the chosen directory were processed.");