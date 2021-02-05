/*
 * Original code by James Shelford, modified by Stephen Royle
 * Requires ROI_csv.py to be in Fiji.app/plugins
 * Make measurements of spindle parameters in Fiji for analysis in R
 */


macro "Spindle analytics" {
	inputFolder = getDirectory("Choose input folder");
	//outputFolder = getDirectory("Choose output folder for the results");
	outputFolder = inputFolder;
	list = getFileList(inputFolder);
	
	for(i = 0; i < list.length; i ++) {
		path = inputFolder + list[i];
		if(checkPath(inputFolder, list[i], outputFolder) == 1)	{
			open(path);
		
			outputPath = outputFolder+list[i];
	 		// The following two lines remove the file extension
			fileExtension = lastIndexOf(outputPath,"."); 
			if(fileExtension != -1)
				outputPath = substring(outputPath,0,fileExtension);
			
			// Z project
			getDimensions(ww, hh, cc, ss, ff);
			if( ss > 1)
				run("Z Project...", "projection=[Max Intensity]");
			for (j = 0; j < cc; j++) {
				Stack.setChannel(j+1);
				run("Enhance Contrast", "saturated=0.35");
			}
			// Select the centrosomes
			roiManager("reset");
			setTool("multipoint");
			run("ROI Manager...");
			Stack.setChannel(cc); // set to last channel
			waitForUser("Select the two centrosomes");
			roiManager("Add");
			roiManager("Select", 0);
			roiManager("rename", outputPath + "_centrosomes");
			run("Select None");
			
			// Draw around the cell
			setTool("freehand");
			waitForUser("Draw around the cell");
			roiManager("Add");
			roiManager("select", 1);
			roiManager("rename", outputPath + "_outline");
			run("Select None");
			
			// Draw a line through the DNA 
			setTool("line");
			Stack.setChannel(1); // set to last channel
			waitForUser("Draw a line through the DNA");
			roiManager("Add");
			roiManager("select", 2);
			roiManager("rename", outputPath + "_DNA");
			run("Select None");
			
			// Running Erick's code to convert ROI to .csv
			run("ROI csv", "output = " + outputFolder );
			run("Close All");
		}
	}
}

function checkPath(dir, name, outdir) {
	path = dir + name;
	if(endsWith(path,".tif")) {
		outputPath = outdir + name;
		// The following two lines remove the file extension
		fileExtension = lastIndexOf(outputPath,"."); 
		if(fileExtension != -1)
			outputPath = substring(outputPath,0,fileExtension);
		// form potential file path to check for
		possibleFilePath = outputPath + "_centrosomes.csv";
		list = getFileList(outdir);
		for(i = 0; i < list.length; i ++) {
			if(outdir + list[i] == possibleFilePath) {
				return 0; // it has a companion csv file - skip
			}
		}
		return 1; // it doesn't have a companion csv file
	}
	else {
		return -1; // it wasn't a tif file
	}
}
