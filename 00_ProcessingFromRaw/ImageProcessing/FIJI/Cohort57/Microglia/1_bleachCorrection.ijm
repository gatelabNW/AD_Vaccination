// Standardizing microglia signal by performing bleach correction and background subtraction

function bleachCorrection(origDir, outputDir) {
	print("Function running");
	orig_FL = getFileList(origDir);
	
	for (i = 0; i < 1; i++) {
		run("Bio-Formats Importer", "open=[" + origDir + orig_FL[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		current = orig_FL[i]; // assumes files are name simply with their sampleName
		original = getTitle();
		print("Processing " + current);
		
		// make output folders
		bleachedImages = outputDir + "bleachedImages/";
		backgroundSubImages = outputDir + "backgroundSubtracted/";
		File.makeDirectory(bleachedImages);
		File.makeDirectory(backgroundSubImages);
		print("Made output folders");

		run("Montage to Stack...", "columns=25 rows=25 border=0");
		print("Made stack");
		
		numSlices = nSlices;
		
		// label slices in order
		for (j = 1; j <= numSlices; j++) {
			setSlice(j);
			setMetadata("Label",""+j);
		}
		print("Labeled slices");
		
		// make middle slides first to ensure the standard is representative of the tissue 
		run("Select All");
		run("Make Substack...", "slices=312-313 delete");
		run("Concatenate...", "  title=shuffled keep image1=[Substack (312-313)] image2=Stack");
		print("Reordered slices");
		
		// standardize signal to noise in all panels
		run("Bleach Correction", "correction=[Histogram Matching]");
		print("Bleaching done");
		
		// re order slices by label
		run("Make Substack...", "slices=3-313 delete");
		run("Concatenate...", "keep image1=[Substack (3-313)] image2=shuffled");
		print("Stack returned to original order");
		
		// re stitch into one image
		run("Make Montage...", "columns=25 rows=25 scale=1");
		saveAs("tiff", bleachedImages + current + "_bleached");
		print("Bleached image saved");
		
		// subtract background
		run("Subtract Background...", "rolling=3");
		saveAs("tiff", backgroundSubImages + sampleName + "_rollingBall_bleached.tif");
		
		close("*");
	}
}

iba1Orig = getDirectory("Navigate to the original Iba1 images ");
	// If necessary, isolate the iba1 channels from original microscope output files and move them all to this same folder
outputFolder = getDirectory("Navigate to the desired output folder ");

bleachCorrection(iba1Dir, outputFolder);

print("Script finished");