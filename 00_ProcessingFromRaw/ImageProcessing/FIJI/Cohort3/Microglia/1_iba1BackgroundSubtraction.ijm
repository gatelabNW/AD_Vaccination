// Performing rollingBall background subtraction on raw Iba1 images

function iba1BackgroundSubtraction(rawIba1Images, outputDir) {
	print("Function running...");
	fileList = getFileList(rawIba1Images);
	roiCounter = 0;
	
	for (i = 0; i < fileList.length; i++) {
		// open bleached iba1 image
		run("Bio-Formats Importer", "open=[" + rawIba1Images + fileList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		current = thresholdList[i];
		sampleName = current.substring(0, 2); // change if necessary
		print("Processing " + sampleName + "...");

		// subtract background
		run("Subtract Background...", "rolling=10");
		saveAs("tiff", backgroundSubImages + sampleName + "_rollingBall.tif");
		
		close("*");
		
		print("Done with " + sampleName + "\n");
		
	}
}

// run it!
rawDir = getDirectory("Navigate to your raw Iba1 images ");
outputDir = getDirectory("Navigate to your output folder ");

iba1BackgroundSubtraction(rawDir, outputDir);
print("Script finished");