// Makes multipages for IF images

// NOTE
// Modifications to this code may be required depending on size of original images
// SpaceRanger will not properly process any multipage tifs that are bigger than 4GB
// May be necessary to split these images across two or three multipage tifs, depending
// on how many channels you want to sync with the spatial data
// Make sure that alignment image is always the first channel of all multipages

// You are at liberty to decide what channels you include in these multipage tifs
// and in what order

function ifMultipage(origDir, dapiDir, iba1BinDir, iba1RawDir, amyloidAllBinDir, amyloidPlaqueDistanceDir, amyloidVesselDistanceDir, outputDir) {
	print("Function running");
	orig_FL = getFileList(origDir);
	dapi_FL = getFileList(dapiDir);
	iba1Bin_FL = getFileList(iba1BinDir);
	iba1Raw_FL = getFileList(iba1RawDir);
	amyloidAllBin_FL = getFileList(amyloidAllBinDir);
	amyloidPlaqueDistance_FL = getFileList(amyloidPlaqueDistanceDir);
	amyloidVesselDistance_FL = getFileList(amyloidVesselDistanceDir);
	
	for (i = 0; i < orig_FL.length; i++) {
		run("Bio-Formats Importer", "open=[" + origDir + orig_FL[i] + "] autoscale color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		run("Flatten");
		orig = getTitle();
		sampleName = orig.substring(0, 2); // change if necessary
		print("Processing " + sampleName + "...");
		run("8-bit");
		print("Opened and flattened original");
		
		// open DAPI binary
		run("Bio-Formats Importer", "open=[" + dapiDir + dapi_FL[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		dapi = getTitle();
		run("8-bit");
		print("Opened DAPI binary");
	
		// open cleaned Iba1 binary
		run("Bio-Formats Importer", "open=[" + iba1BinDir + iba1Bin_FL[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		binaryIba1 = getTitle();
		run("8-bit");
		print("Opened cleaned Iba1 binary");
		
		// open raw Iba1 w binary selection
		run("Bio-Formats Importer", "open=[" + iba1RawDir + iba1Raw_FL[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		rawIba1 = getTitle();
		run("8-bit");
		print("Opened raw Iba1 w binary selection");
		
		// open cleaned amyloid binary
		run("Bio-Formats Importer", "open=[" + amyloidAllBinDir + amyloidAllBin_FL[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		binaryAllAmyloid = getTitle();
		run("8-bit");
		print("Opened cleaned amyloid binary");
		
		// open amyloid distance map
		run("Bio-Formats Importer", "open=[" + amyloidPlaqueDistanceDir + amyloidPlaqueDistance_FL[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		plaqueDist = getTitle();
		run("8-bit");
		print("Opened amyloid distance map");
		
		// open raw amyloid w binary selection 
		run("Bio-Formats Importer", "open=[" + amyloidVesselDistanceDir + amyloidVesselDistance_FL[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		vesselDist = getTitle();
		run("8-bit");
		print("Opened raw amyloid w binary selection");
		
		// when all images opened and processed
		// Set your grayscale settings here
	    run("Merge Channels...", "c1=[" + orig + "] c2=[" + dapi + "]  c3=[" + binaryIba1 + "] c4=[" + rawIba1 + "] c5=[" + binaryAllAmyloid + "] c6=[" + plaqueDist + "] c7=[" + vesselDist + "] create");
	    print("Ran merge channels");
	    Stack.setDisplayMode("grayscale");
	    Stack.setChannel(1);
	    Stack.setChannel(2);
	    Stack.setChannel(3);
	    Stack.setChannel(4);
	    Stack.setChannel(5);
	    Stack.setChannel(6);
	    Stack.setChannel(7);
	    print("Ran grayscale");
	    saveAs("tiff", outputDir + sampleName + "_IF_multipage_new");
	    print("Saved " + sampleName + "_IF_multipage_new.tif");
	    close("*");
	    print("Done with " + sampleName + "\n");
	}	
}

// run it!
original = getDirectory("Navigate to your original tifs ");
dapi_dir = getDirectory("Navigate to your binary DAPI images ");
iba1_binary = getDirectory("Navigate to your cleaned binary Iba1 images ");
iba1_raw = getDirectory("Navigate to your cleaned bleached and background subtracted Iba1 images ");
amyloid_all_binary = getDirectory("Navigate to your cleaned binary (plaque+vessel combined) amyloid images ");
amyloid_plaque_distance = getDirectory("Navigate to your amyloid plaque distance maps ");
amyloid_vessel_distance = getDirectory("Navigate to your amyloid vessel distance maps ");
output_directory = getDirectory("Navigate to your output directory ");
ifMultipage(original, dapi_dir, iba1_binary, iba1_raw, amyloid_all_binary, amyloid_plaque_distance, amyloid_vessel_distance, output_directory);

print("Script finished");