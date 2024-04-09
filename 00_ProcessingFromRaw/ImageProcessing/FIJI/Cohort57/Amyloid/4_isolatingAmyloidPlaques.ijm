// ImageJ macro for using the cleaned amyloid vessel channels to isolate
// the amyloid plaque signals in the original binary IF amyloid images

function isolatingAmyloidPlaques(cleanedVesselProbMap, amyloidBinary, outputDir) {
	print("Function running...");
	vessel_FL = getFileList(cleanedVesselProbMap);
	allAmyloid_FL = getFileList(amyloidBinary);
	
	for (i = 0; i < imageFileList.length; i++) {
		run("Bio-Formats Importer", "open=[" + cleanedVesselProbMap + vessel_FL[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	    current = imageFileList[i];
	    sampleName = current.substring(0, 2);
	    print("Processing " + sampleName + "...");
	    print("Cleaned vessel image opened");
	    
	    for (j = 0; j < 10; j++) {
	    	run("Dilate");
	    }
	    
	    run("Create Selection");
	    roiManager("Add");
	    print("Dilated vessels and created selection");
	    
	    run("Bio-Formats Importer", "open=[" + amyloidBinary + allAmyloid_FL[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	    print("Opened comprehensive binary amyloid image");
	    
	    roiManager("Select", 0);
	    run("Clear");
	    print("Cleared vessels from binary image to isolate plaques");
		
	    saveAs("tiff", output + sampleName + "_isolated_plaque_signal.tif");
	    print("Saved image");
	    print("Done with " + sampleName + "\n");
	    
	    close("*");
	    close("ROI Manager");
	}
}

// run it!
probabilityMaps = getDirectory("Navigate to your thresholded vessel probability maps ");
binaryAmyloidImages = getDirectory("Navigate to the comprehensive binary IF amyloid images ");
outputDirectory = getDirectory("Navigate to the desired output folder ");

isolatingAmyloidPlaques(probabilityMaps, binaryAmyloidImages, outputDirectory);
print("Script finished");