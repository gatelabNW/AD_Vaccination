// Ensures there are no overlapping regions between plaque and vessel distance maps
// Vessel distance maps take priority and are erased from plaque distance maps
// If there is a plaque present within the distance map of a vessel, that binary 
// plaque signal is preserved

function erasingOverlapDistanceMaps(vessel5thLayer, plaqueBinaryDir, plaqueDistDir, newPlaqueDistDir) {
	print("Function running...");
	vessel_FL = getFileList(vessel5thLayer);
	plaqueBin_FL = getFileList(plaqueBinaryDir);
	plaqueDist_FL = getFileList(plaqueDistDir)
	
	for (i = 0; i < vessel_FL.length; i++) {
		open(vesselLayer100 + vessel_FL[i]);
		vessel = getTitle();
		sampleName = vessel.substring(7,9); // change if necessary
		print("Processing " + sampleName);
		
		open(plaqueDisDir + plaqueDist_FL[i]);
		plaqueDist = getTitle();
		
		open(plaqueBinaryDir + plaqueBin_FL[i]);
		plaqueBin = getTitle();
		
		selectImage(vessel);
		run("Select None");
		run("Create Selection");
		roiManager("Add"); // 0
		run("Select None");
		print("Created selection of 5th layer of vessel distance map");
		
		selectImage(plaqueBin);
		run("Select None");
		run("Create Selection");
		roiManager("Add"); // 1
		run("Select None");
		print("Created selection of binary plaques");
		
		// removing binary plaque signals so these are lost during the vessel distance map removal
		selectImage(plaqueDist);
		roiManager("Select", 1);
		run("Clear");
		run("Select None");
		
		// removing vessel distance map signal from plaque distance map
		roiManager("Select", 0);
		run("Clear");
		run("Select None");
		
		// restoring binary plaque signal on map
		run("Invert LUT");
		roiManager("Select", 1);
		run("Clear");
		run("Select None");
		run("Invert LUT");
		
		print("Cleared vessel distance maps from plaque distance maps");
		print("Retained all plaques within vessel distance maps");
		selectImage(plaqueDist);
		rename(sampleName + "_plaqueDM_minus_vesselDM_plaques_retained");
		newTitle = getTitle();
		print(getTitle());

		saveAs("tiff", newPlaqueDistDir + newTitle);
		print("Saved isolated plaque distance map\n");
		close("*");
		close("ROI Manager");
	}
}

// run it!
vesselDistance5thLayer = getDirectory("Navigate to the folder holding the 5th layer of the vessel distance maps");
plaqueBinary = getDirectory("Navigate to the binary plaque signal (that was used for distance map creation)");
plaqueDistanceMaps = getDirectory("Navigate to the binary plaque distance maps");
output = getDirectory("Navigate to the desired output folder");
erasingOverlapDistanceMaps(vesselDistance5thLayer, plaqueBinary, plaqueDistanceMaps, output);

print("Script finished");