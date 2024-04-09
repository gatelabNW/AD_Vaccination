// Creating multipage and making all channels 8bit and grayscale
// Assumes that all pieces of the multipage are contained within one main folder that holds folders entitled
	// Decon-H (deconvoluted hematoxylin images)
	// Decon-DAB (deconvoluted DAB images)
	// binary_amyloid (uncleaned binary amyloid images; cleaning takes place in this script)
	// SplitTiff (first channel of original microscope image, only used for alignment on SpaceRanger)

function makeMultipage_eosin(mainFolder) {
	print("Function running...");
    fileList = getFileList(currentDir);
    roiCounter = 0;
    
    for (i = 0; i < fileList.length; i++) {
        sampleName = fileList[i].substring(0, fileList[i].length - 5); // change if necessary
        print("Processing " + sampleName + "...");
        
        //make new folders
        multipageOutput = mainFolder + "final_multipage/";
        deconDAB_output = mainFolder + "decon_DAB_cleaned/";
        rawCleanedDABOutput = mainFolder + "thresholded_binary/";
        File.makeDirectory(multipageOutput);
        File.makeDirectory(deconDAB_output);
        File.makeDirectory(rawCleanedDABOutput);
    	
    	//deconvoluted hematoxylin
    	run("Bio-Formats Importer", "open=[" + mainFolder + "Decon-H/" + sampleName + "-H.tif] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
    	run("8-bit");
    	run("Invert");
    	H = getTitle();
    	print("Opened H");
    	print("Ran 8-bit");
    	
    	//deconvoluted DAB
    	run("Bio-Formats Importer", "open=[" + mainFolder + "Decon-DAB/" + sampleName + "-DAB.tif] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
    	run("8-bit");
    	DAB = getTitle();
    	run("Duplicate...", "title=[" + sampleName + "-DAB_duplicate]");
    	run("Invert");
    	DAB_dupe = getTitle();
    	selectImage(DAB);
    	run("Invert");
    	print("Opened DAB");
    	print("Ran 8-bit");
    	
    	//clean binary amyloid
    	run("Bio-Formats Importer", "open=[" + mainFolder + "binary_amyloid/" + sampleName + "-uncleaned.tif] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
    	run("Invert");
    	uncleaned = getTitle();
    	run("Duplicate...", "title=[" + sampleName + "-uncleaned]");
    	uncleaned_dupe = getTitle();
    	selectImage(uncleaned);
    	run("Dilate");
    	run("Erode");
    	print("Opened binary amyloid");
    	print("Done with dilation/erosion");
    	run("Analyze Particles...", "size=15-infinity show=Masks clear");
    	print("Analyzed particles");
    	run("Create Selection");
    	roiManager("Add");
    	
    	// restore selection onto uncleaned binary amyloid
    	selectImage(uncleaned_dupe);
    	roiManager("Select", roiCounter);
    	run("Restore Selection");
    	run("Clear");
    	print("Restored selection on uncleaned binary + cleared outside");
    	saveAs("tiff", rawCleanedDABOutput + sampleName + "-cleaned.tif");
    	run("8-bit");
    	amyloid = getTitle();
    	run("Create Selection");
    	roiManager("Add");
    	roiCounter++; 
    	print("Cleaned binary amyloid");
    	print("Ran 8-bit");
    	
    	//add selection on Decon_DAB
    	selectImage(DAB_dupe);
    	roiManager("Select", roiCounter);
    	run("Restore Selection");
    	print("Restored selection");
    	run("Clear Outside");
    	print("Restored selection on Decon_DAB + cleared outside");
    	saveAs("tiff", deconDAB_output + sampleName + "-DAB_cleaned.tif");
    	run("8-bit");
    	DAB_cleaned = getTitle();
    	roiCounter++;
    	print("Made cleaned DAB image");
    	print("Ran 8-bit");
    	
    	//green channel of original tif
    	run("Bio-Formats Importer", "open=[" + mainFolder + "splitTiff/" + sampleName + "-splitTif.tif] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
    	run("8-bit");
    	run("Invert");
    	splitTiff = getTitle();
    	print("Opened raw image for alignment");
    	print("Ran 8-bit");
    	
    	// Set your grayscale settings here
        run("Merge Channels...", "c1=[" + splitTiff + "] c2=[" + amyloid + "] c3=[" + DAB_cleaned + "] c4=[" + DAB + "] c5=[" + H + "] create");
        print("Ran merge channels");
        Stack.setDisplayMode("grayscale");
        Stack.setChannel(1);
        Stack.setChannel(2);
        Stack.setChannel(3);
        Stack.setChannel(4);
        Stack.setChannel(5);
        print("Ran grayscale");
        saveAs("tiff", multipageOutput + sampleName + "_multi_final.tif");
        print("Saved " + sampleName + "_multi_final.tif");
        close("*");
        print("Done with " + sampleName + "\n");    
    }
}

// run it!
piecesMainDir = getDirectory("Navigate to the folder containing all of the multipage pieces ");
makeMultipage(piecesMainDir);

print("Script finished");