// Deconvoluting original image & saving multipage components

function gettingMultipagePieces_eosin(originalTif, outputDir) {
	print("Function running...");
    fileList = getFileList(originalTif);
    for (i = 0; i < fileList.length; i++) {
        sampleName = fileList[i].substring(0, fileList[i].length - 23); // change if necessary 
        
        // Make output folders
        amyloidOutput = outputDir + "amyloid/";
        hematoxylinOutput = outputDir + "Decon-H/";
        dabOutput = outputDir + "Decon-DAB/";
        splitGreenOutput = outputDir + "splitTiff/"; // alignment images
		File.makeDirectory(amyloidOutput);
		File.makeDirectory(hematoxylinOutput);
		File.makeDirectory(dabOutput);
		File.makeDirectory(splitGreenOutput);
        
    	// Open and process the tif files
    	run("Bio-Formats Importer", "open=[" + originalTif + fileList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
    	full_tiff = getTitle();
    	print("Processing " + sampleName);
    	run("RGB Color");
    	run("Colour Deconvolution", "vectors=[H DAB]");
    	print("Deconvoluted");
    	
    	// Select irrelevant channel and close
    	selectImage(fileList[i] + " (RGB)-(Colour_3)");
    	close();
    	
    	// Hematoxylin channel
    	selectImage(fileList[i] + " (RGB)-(Colour_1)");
    	saveAs("tiff", hematoxylinOutput + sampleName + "-H.tif");
    	print("Saved hematoxylin channel");
    	
    	// DAB channel
    	selectImage(fileList[i] + " (RGB)-(Colour_2)");
    	saveAs("tiff", dabOutput + sampleName + "-DAB.tif");
    	print("Saved DAB channel");
    	
    	// Duplicate DAB channel to process amyloid image
    	run("Duplicate...", "title=[Amyloid Processing]");
    	saveAs("tiff", amyloidOutput + sampleName + "-amyloid.tif");
    	print("Saved duplicated DAB channel for thresholding amyloid signal");
    	close("*");
    	
    	// Split original tiff and save green channel for alignment image
    	run("Bio-Formats Importer", "open=[" + currentDir + fileList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
    	print("Re-opened: " + full_tiff);
    	run("Split Channels");
    	print("Ran split channels");
    	selectImage("C3-" + full_tiff);
    	close();
    	selectImage("C1-" + full_tiff);
    	close();
    	selectImage("C2-" + full_tiff);
    	saveAs("tiff", splitGreenOutput + sampleName + "-splitTiff.tif");
    	print("Saved alignment image for SpaceRanger");
	}
}

// run it!
original = getDirectory("Navigate to your original tifs ");
output = getDirectory("Navigate to your output folder ");
gettingMultipagePieces_eosin(original, output);

print("Script finished");