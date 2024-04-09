// Creating distance maps of cleaned binary amyloid images for neighboring Visium spot detection

function amyloidDistanceChannels(amyloidBinary, outputDir) {
	print("Function running...");
	amyloidFileList = getFileList(amyloidBinary);
	roiCounter = 0;
	
	for (i = 0; i < amyloidFileList.length; i++) {
		// open rollingBall iba1 image image
		run("Bio-Formats Importer", "open=[" + amyloidBinary + amyloidFileList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		original = getTitle();
		sampleName = original.substring(0, original.length - 12); // change if necessary
		print("Processing " + sampleName + "...");
	    
	    saveAs("tiff", outputDir + sampleName + "_cleaned");
	    print("Saved cleaned");
		
		dilated1st = outputDir + "_dilated1st_layer/";
		dilated2nd = outputDir + "_dilated2nd_layer/";
		dilated3rd = outputDir + "_dilated3rd_layer/";
		dilated4th = outputDir + "_dilated4th_layer/";
		dilated5th = outputDir + "_dilated5th_layer/";
		File.makeDirectory(dilated1st);
		File.makeDirectory(dilated2nd);
		File.makeDirectory(dilated3rd);
		File.makeDirectory(dilated4th);
		File.makeDirectory(dilated5th);
		
		// create all dilated images
		selectImage(original);
		nextDilation(roiCounter, original, "dilated_1st", 1);
		roiCounter++;
		saveAs("tiff", dilated1st + type + "_" + sampleName + "_dilated1st");
		d1st = getTitle();
		print("Created 1st dilation layer");
		
		nextDilation(roiCounter, d1st, "dilated_2nd", 2);
		roiCounter++;
		saveAs("tiff", dilated2nd + type + "_" + sampleName + "_dilated2nd");
		d2nd = getTitle();
		print("Created 2nd dilation layer");
		
		nextDilation(roiCounter, d2nd, "dilated_3rd", 3);
		roiCounter++;
		saveAs("tiff", dilated3rd + type + "_" + sampleName + "_dilated3rd");
		d3rd = getTitle();
		print("Created 3rd dilation layer");
		
		nextDilation(roiCounter, d3rd, "dilated_4th", 4);
		roiCounter++;
		saveAs("tiff", dilated4th + type + "_" + sampleName + "_dilated4th");
		d4th = getTitle();
		print("Created 4th dilation layer");
		
		nextDilation(roiCounter, d4th, "dilated_5th", 5);
		roiCounter++;
		saveAs("tiff", dilated5th + type + "_" + sampleName + "_dilated5th");
		d5th = getTitle();
		print("Created 5thdilation");
		
		// ensure no overlapping regions
		selectImage(d20);
		clearPreviousAndDivide(5*i, 1.2);
		d1st = getTitle();
		
		selectImage(d40);
		clearPreviousAndDivide((5*i)+1, 1.5);
		d2nd = getTitle();
		
		selectImage(d60);
		clearPreviousAndDivide((5*i)+2, 2);
		d3rd = getTitle();
		
		selectImage(d80);
		clearPreviousAndDivide((5*i)+3, 3);
		d4th = getTitle();
		
		selectImage(d100);
		clearPreviousAndDivide((5*i)+4, 6);
		d5th = getTitle();
		
		imageCalculator("Add create", original, d1st);
		result = "Result of " + current;
		imageCalculator("Add create", result, d2nd);
		result = "Result of " + result;
		imageCalculator("Add create", result, d3rd);
		result = "Result of " + result;
		imageCalculator("Add create", result, d4th);
		result = "Result of " + result;
		imageCalculator("Add create", result, d5th);
		saveAs("tiff", outputDir + "_" + sampleName + "_distance_map");
		print("Done with " + "_" + sampleName);
		close("*");
	}
}

function nextDilation(counter, latestImage, newestImage, iteration){
	run("Duplicate...", "title=[dilate_" + iteration*20 + "x]");
	duplicate = getTitle();
	print("Duplicated " + latestImage);
	
	selectImage(latestImage);
	run("Create Selection");
	roiManager("Add");
	print("Created selection");
	run("Select None");
	
	selectImage(duplicate);
	for (j=0; j < 20; j++) {
		run("Dilate");
		print(j);
	}
	print("Completed dilations");
	
	rename(newestImage);
}

function clearPreviousAndDivide(counter, n){
	roiManager("Select", counter)
	run("Clear");
	run("Select None");
	run("Divide...", "value=" + n);
}

// run it!
amyloidDirectory = getDirectory("Navigate to the images of interest (ex. cleaned binary amyloid plaque) ");
outputDirectory = getDirectory("Navigate to the output folder ");

amyloidDistanceChannels(amyloidDirectory, outputDirectory);
print("Script finished");