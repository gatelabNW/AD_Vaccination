// Cleaning binary Iba1 images and isolating included selections on bleached, 
// background subtracted images

// Assumes Iba1 images have undergone rollingBall background subtraction with
// a radius of 10 and then have been manually thresholded

function iba1Processing(rollingBallImages, binaryImages, outputDir) {
	print("Function running...");
	fileList = getFileList(rollingBallImages);
	thresholdList = getFileList(binaryImages);
	roiCounter = 0;
	
	for (i = 0; i < fileList.length; i++) {
		// open bleached iba1 image
		run("Bio-Formats Importer", "open=[" + rollingBallImages + fileList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		currentBleach = thresholdList[i];
		sampleName = currentBleach.substring(0, 2); // change if necessary
		bleachedImage = fileList[i];
		print("Processing " + sampleName + "...");
		
		// open binary iba1 image
		run("Bio-Formats Importer", "open=[" + binaryImages + thresholdList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		run("Invert");
		orig_binary = getTitle();
		run("Duplicate...", "title=[" + sampleName + "_binary_duplicate]");
		binary_duplicate = getTitle();
		selectImage(orig_binary);
		
		// clean up through a series of dilations, erosions, and filtering
		run("Dilate");
		run("Dilate");
		run("Erode");
		run("Erode");
		run("Analyze Particles...", "size=50-Infinity pixel show=Masks clear");
		run("Invert LUT");
		
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
		saveAs("tiff", binaryIba1Output + sampleName + "_rollingBall_binary_cleaned.tif");
		cleanedMask = getTitle();
		print("Created cleaned binary image");
		roiCounter++;
		
		// create selection on saved binary image
		selectImage(cleanedMask);
		run("Create Selection");
		roiManager("Add");
		selectImage(sampleName + "_rollingBall_binary_cleaned.tif"); 
		close();
		
		//restore selection on duplicated iba1
		selectImage(bleachedImage);
		roiManager("Select", roiCounter);
		run("Clear Outside"); // NEED TO MAKE SURE BACKGROUND = 0 AFTER THIS
		run("Select None");
		saveAs("tiff", cleanedIba1Output + sampleName + "_rollingBall_cleaned.tif");
		print("Created cleaned rollingBall image");
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
bleachedDir = getDirectory("Navigate to your rollingBall images ");
binaryDir = getDirectory("Navigate to your binary Iba1 images ");
outputDir = getDirectory("Navigate to your output folder ");

iba1Processing(bleachedDir, binaryDir, outputDir);
print("Script finished");