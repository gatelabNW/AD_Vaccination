// Performs size correction on images to account for alteration in dimensions 
// of microglia images during bleach correction
// Ensures all multipage components will be the same size and are thus compatible
// to be combined into one file

function sizeCorrection(inputFolder) {
	print("Function running");
	images_FL = getFileList(inputFolder);
	
	for (i = 0; i < images_FL.length; i++) {
		run("Bio-Formats Importer", "open=[" + inputFolder + images_FL[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		current = images_FL[i];
		sampleName = current.substring(0, 2);
		original = getTitle();
		print("Processing " + sampleName);
		
		// make output folders
		outputFolder = inputFolder + "_un-re-stitched";
		File.makeDirectory(outputFolder);
		print("Made output folder");

		selectImage(original);
		run("Montage to Stack...", "columns=25 rows=25 border=0");
		print("Made stack");
		
		// re stitch into one image
		run("Make Montage...", "columns=25 rows=25 scale=1");
		
		saveAs("tiff", outputFolder + "/" + current + "_un-re-stitched");
		close();
		print("Saved image that has been unstitched and restitched");
		
		print("Done with " + sampleName + "\n");
		
		close("*");
	}
}

// run it!
imagesOfInterest = getDirectory("Navigate to the folder holding the images you wish to size correct ");
sizeCorrection(imagesOfInterest);

print("Script finished");