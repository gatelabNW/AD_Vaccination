// ImageJ macro for cleaning manually thresholded amyloid vessel probability maps

function cleaningProbabilityMap(vesselProbMap, output) {
	print("Function running...");
	imageFileList = getFileList(vesselProbMap);
	
	for (i = 0; i < imageFileList.length; i++) {
		run("Bio-Formats Importer", "open=[" + vesselProbMap + imageFileList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	    current = imageFileList[i];
	    sampleName = current.substring(0, 2);
	    print("Processing " + sampleName + "...");
	    print("Image opened");
	    
	    for (j = 0; j < 3; j++) {
	    	run("Dilate");
	    }
	    for (j = 0; j < 3; j++) {
	    	run("Erode");
	    }
	    print("Completed dilations/erosions");
		
		run("Analyze Particles...", "size=500-infinity show=Masks clear");
		run("Invert LUT");
		print("Filtered");
		
	    saveAs("tiff", output + sampleName + "_cleaned_thresholded_blurred_prob_map.tif");
	    print("Saved image");
	    print("Done with " + sampleName + "\n");
	    
	    close("*");
	}
}

// run it!
probabilityMaps = getDirectory("Navigate to your thresholded vessel probability maps ");
outputDirectory = getDirectory("Navigate to the desired output folder ");

cleaningProbabilityMap(probabilityMaps, outputDirectory);
print("Script finished");
print("Now time to manually clean the rest :)!");