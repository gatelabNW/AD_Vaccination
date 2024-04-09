// Cleaning manually thresholded amyloid IF images

function cleaningAmyloidBinary(amyloidBinary, outputDir) {
	print("Function running...");
	amyloidFileList = getFileList(amyloidBinary);
	roiCounter = 0;
	
	for (i = 0; i < amyloidFileList.length; i++) {
		// dilation erosion
	    run("Dilate");
	    run("Dilate");
	    run("Erode");
	    run("Erode");
	    print("Ran dilation/erosion");
	    
	    run("Analyze Particles...", "size=82-infinity show=Masks clear");
	    	// adjust lower bound based on image resolution
	    run("Invert LUT");
	    print("Filtered");
	    
	    saveAs("tiff", outputDir + "_" + sampleName + "_cleaned_binary");
		print("Done with " + "_" + sampleName);
		close("*");
	}
}

// run it!
amyloidDirectory = getDirectory("Navigate to the images of interest (ex. cleaned binary amyloid plaque) ");
outputDirectory = getDirectory("Navigate to the output folder ");

cleaningAmyloidBinary(amyloidDirectory, outputDirectory);
print("Script finished");