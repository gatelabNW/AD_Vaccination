// ImageJ macro for isolating amyloid vessel probability maps

// Assumes you have manually run the classifiers on the appropriate images,
// calculated the probability maps, show the probability maps in FIJI,
// converted them to 8-bit, and saved them to the same folder

function isolatingAndBlurringProbMap(vesselProbMap, output) {
	print("Function running...");
	imageFileList = getFileList(vesselProbMap);
	
	for (i = 0; i < imageFileList.length; i++) {
		run("Bio-Formats Importer", "open=[" + vesselProbMap + imageFileList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	    current = imageFileList[i];
	    sampleName = current.substring(0, 2);
	    print("Processing " + sampleName + "...");
	    print("Image opened");
	    
	    // isolate probability map of interest assuming channels in classifier are in the following order
	    	// background
	    	// foreground
	    	// amyloid_plaque
	    	// amyloid_vessel
	    	// neuron
	    	// microglia
	    	// dapi
	    	
	    run("Split Channels"); 
	    print("Split channels");
	    
	    selectImage("C1-" + current);
	    close();
	    selectImage("C2-" + current);
	    close();
	    selectImage("C3-" + current);
	    close();
	    selectImage("C4-" + current);
	    amyloid_vessel = getTitle();
	    selectImage("C5-" + current);
	    close();
	    selectImage("C6-" + current);
	    close();
	    selectImage("C7-" + current);
	    close();
	    print("Isolated probability map for amyloid in vessels");
	    
	    run("Despeckle");
		run("Gaussian Blur...", "sigma=4");
		print("Despeckle and Gaussian blur done");
		
	    saveAs("tiff", output + sampleName + "_blurred_vessel_prob_map.tif");
	    print("Saved image");
	    print("Done with " + sampleName + "\n");
	    
	    close("*");
	}
}

// run it!
probabilityMaps = getDirectory("Navigate to the probability maps");
outputFolder = getDirectory("Navigate to the desired output folder");
isolatingAndBlurringProbMap(probabilityMaps, outputFolder);

print("Script finished");