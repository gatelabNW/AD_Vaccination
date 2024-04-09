// Cleaning binary Iba1 images and isolating included selections on bleached, 
// background subtracted images

// Assumes images have been manually thresholded after undergoing the 
// bleach correction and background subtraction

function iba1Processing(bleach, binaryPath, output) {
	print("Function running...");
	fileList = getFileList(bleach);
	thresholdList = getFileList(binaryPath);
	roiCounter = 0;
	
	cleanedIba1Output = output + "bleached_rB_manualThresh_iba1_cleaned/";
	binaryIba1Output = output + "bleached_rB_manualThresh_binary_iba1_cleaned/";
	File.makeDirectory(cleanedIba1Output);
	File.makeDirectory(binaryIba1Output);
	print("Made output folders");
	
	for (i = 0; i < fileList.length; i++) {
		// open bleached iba1 image
		run("Bio-Formats Importer", "open=[" + bleach + fileList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		currentBleach = thresholdList[i];
		sampleName = currentBleach.substring(0, 2); // change if necessary
		bleachedImage = fileList[i];
		print("Processing " + sampleName + "...");
		
		// open binary iba1 image
		run("Bio-Formats Importer", "open=[" + binaryPath + thresholdList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		run("Invert");
		orig_binary = getTitle();
		run("Duplicate...", "title=[" + sampleName + "_binary_duplicate]");
		binary_duplicate = getTitle();
		selectImage(orig_binary);
		
		// clean up through a series of dilations, erosions, and filtering
		dilationsErosions();
		de1 = getTitle();
		run("Analyze Particles...", "size=10-Infinity pixel show=Masks clear");
		run("Invert LUT");
		
		dilationsErosions();
		de2 = getTitle();
		run("Analyze Particles...", "size=15-Infinity pixel show=Masks clear");
	    run("Invert LUT");
		
		dilationsErosions();
		de3 = getTitle();
		run("Analyze Particles...", "size=50-infinity pixels show=Masks clear");
			// adjust lower bound based on image resolution
		run("Invert LUT");
		mask = getTitle();
		print("Completed cleaning of binary image");
		
		selectImage(orig_binary);
		close();
		selectImage(de1);
		close();
		selectImage(de2);
		close();
		selectImage(de3);
		close();
		
		// create selection on binary image
		selectImage(mask);
		run("Create Selection");
		roiManager("Add");
		selectImage(mask);
		close();
		
		//restore selection on duplicated binary
		selectImage(binary_duplicate);
		roiManager("Select", roiCounter);
		run("Clear Outside"); // NEED TO MAKE SURE BACKGROUND = 0 AFTER THIS
		run("Select None");
		run("Analyze Particles...", "size=29-infinity pixels show=Masks clear");
			// adjust lower bound based on image resolution
		run("Invert LUT");
		saveAs("tiff", binaryIba1Output + sampleName + "_bleached_binary_cleaned.tif");
		cleanedMask = getTitle();
		print("Created cleaned binary image");
		roiCounter++;
		
		// create selection on saved binary image
		selectImage(cleanedMask);
		run("Create Selection");
		roiManager("Add");
		selectImage(sampleName + "_bleached_binary_cleaned.tif"); 
		close();
		
		//restore selection on duplicated iba1
		selectImage(bleachedImage);
		roiManager("Select", roiCounter);
		run("Clear Outside"); // NEED TO MAKE SURE BACKGROUND = 0 AFTER THIS
		run("Select None");
		saveAs("tiff", cleanedIba1Output + sampleName + "_bleached_cleaned.tif");
		print("Created cleaned bleached image");
		roiCounter++;
		
		close("*");
		
		print("Done with " + sampleName + "\n");
		
	}
}

function dilationsErosions() {
	for (j=0; j < 2; j++) {
			run("Dilate");
		}
		for (j=0; j < 2; j++) {
			run("Erode");
		}
}

// run it!
bleachedDir = getDirectory("Navigate to your bleached images ");
binaryDir = getDirectory("Navigate to your binary Iba1 images ");
outputDir = getDirectory("Navigate to your output folder ");

iba1Processing(bleachedDir, binaryDir, outputDir);
print("Script finished");