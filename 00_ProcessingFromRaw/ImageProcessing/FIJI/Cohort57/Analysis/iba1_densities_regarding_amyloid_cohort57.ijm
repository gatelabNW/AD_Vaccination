// determining Iba1+ density overlapping amyloid, neighboring amyloid, and not neighboring amyloid
// Cohort 3

function iba1DensitiesRegardingAmyloid(iba1Dir, greyMatterRoiDir, binaryPlaqueDir, binaryVesselDir, distanceLayer100PlaqueDir, distanceLayer100VesselDir, dapiDir, output, excelName, excelPath) {
	print("Function running...");
	iba1FileList = getFileList(iba1Dir);
	greyMatterRoiFileList = getFileList(greyMatterRoiDir);
	vesselFileList = getFileList(binaryVesselDir);
	plaqueFileList = getFileList(binaryPlaqueDir);
	distancePlaqueFileList = getFileList(distanceLayer100PlaqueDir);
	distanceVesselFileList = getFileList(distanceLayer100VesselDir)
	dapiFileList = getFileList(dapiDir);
		
	// also measure DAPI coverage in all of these regions
	
	for (i = 0; i < iba1FileList.length; i++) {
		// open cleaned binary iba1
		run("Bio-Formats Importer", "open=[" + iba1Dir + iba1FileList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	    currentIba1 = getTitle();
	    sample = iba1FileList[i];
	    sampleName = sample.substring(0, 2); // may need to be adjusted depending on naming convention
	    print("Processing " + sampleName + "...");
	
		// duplicate six times
		run("Duplicate...", "title=[iba1DirectPlaque]");
		run("Duplicate...", "title=[iba1NeighboringPlaque]");
		run("Duplicate...", "title=[iba1NotNeighboringPlaque]");
		
		run("Duplicate...", "title=[iba1DirectVessel]");
		run("Duplicate...", "title=[iba1NeighboringVessel]");
		run("Duplicate...", "title=[iba1NotNeighboringVessel]");
		
		print("Opened and duplicated cleaned binary Iba1 image");
		
		// open cleaned binary DAPI 
		run("Bio-Formats Importer", "open=[" + dapiDir + dapiFileList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	    run("Invert");
	    run("Invert LUT");
	    currentDAPI = getTitle();
	    
	    // duplicate six times		
		run("Duplicate...", "title=[dapiDirectPlaque]");
		run("Duplicate...", "title=[dapiNeighboringPlaque]");
		run("Duplicate...", "title=[dapiNotNeighboringPlaque]");
		
		run("Duplicate...", "title=[dapiDirectVessel]");
		run("Duplicate...", "title=[dapiNeighboringVessel]");
		run("Duplicate...", "title=[dapiNotNeighboringVessel]");
		
		print("Opened and duplicated cleaned binary DAPI image");
		
		// prepare excel file
		run("Read and Write Excel", "no_count_column file=[" + excelPath + excelName + "] sheet=" + sampleName + "");
		run("Read and Write Excel", "file=[" + excelPath + excelName + "] file_mode=read_and_open");
		
		// open grey matter ROI + add to manager (0)
		open(greyMatterRoiDir + greyMatterRoiFileList[i]);
		roiManager("Add");
		run("Select None");
		selectImage(currentIba1);
		rename("Grey Matter");
		roiManager("Select", 0);
		run("Clear Outside"); // isolating gray matter
		// measure value of total grey matter area
		run("Measure");
		print("Measured whole grey matter area");
		saveAs("tiff", output + sampleName + "_greyMatterSelection.tif");
		
		// measure total Iba1 coverage in grey matter
		rename("Iba1GreyMatter");
		run("Create Selection");
		run("Measure");
		print("Measured total area of Iba1 in grey matter");
		run("Select None");
		close();
		
		// measure total DAPI coverage in grey matter
		selectImage(currentDAPI);
		rename("DAPIGreyMatter");
		roiManager("Select", 0);
		run("Clear Outside"); // isolating gray matter
		run("Select None");
		run("Create Selection");
		run("Measure");
		print("Measured total area of DAPI in grey matter");
		run("Select None");
		close();
		
		// DIRECTLY OVERLAPPING
		// open cleaned binary plaque
		run("Bio-Formats Importer", "open=[" + binaryPlaqueDir + plaqueFileList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	    rename("binaryPlaque");
	    roiManager("Select", 0); // isolating plaques in gray matter
	    run("Clear Outside");
	    run("Select None");
	    run("Create Selection"); // make sure this does not need to be inverted
			// selection of amyloid plaque present in grey matter
		roiManager("Add");
		// measure value for total plaque area (1)
		run("Measure");
		print("Measured total area of plaque in grey matter");
		
		// open cleaned binary vessel
		run("Bio-Formats Importer", "open=[" + binaryVesselDir + vesselFileList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	    rename("binaryVessel");
	    roiManager("Select", 0);
	    run("Clear Outside"); // isolating vessels in gray matter
	    run("Select None");
	    run("Create Selection"); // make sure this does not need to be inverted
			// selection of amyloid vessel present in grey matter
		roiManager("Add");
		// measure value for total vessel area (2)
		run("Measure");
		print("Measured total area of vessel in grey matter");
		
		// select directPlaque cleaned binary iba1 image
		selectImage("iba1DirectPlaque");
		roiManager("Select", 1);
	    run("Clear Outside"); // isolating iba1 directly overlapping amyloid plaque
	    run("Select None");
	    run("Create Selection"); // make sure this does not need to be inverted
			// selection of iba1 directly overlapping plaque in grey matter
		run("Measure");
		print("Measured area of microglia directly on plaque");
		
		// select directVessel cleaned binary iba1 image
		selectImage("iba1DirectVessel");
		roiManager("Select", 2);
	    run("Clear Outside"); // isolating iba1 directly overlapping amyloid vessel
	    run("Select None");
	    run("Create Selection"); // make sure this does not need to be inverted
			// selection of iba1 directly overlapping vessel in grey matter
		run("Measure");
		print("Measured area of microglia directly on vessel");
		
		// select directPlaque cleaned binary DAPI image
		selectImage("dapiDirectPlaque");
		roiManager("Select", 1);
	    run("Clear Outside"); // isolating DAPI directly overlapping amyloid plaque
	    run("Select None");
	    run("Create Selection"); // make sure this does not need to be inverted
			// selection of DAPI directly overlapping plaque in grey matter
		run("Measure");
		print("Measured area of DAPI directly on plaque");
		
		// select directVessel cleaned binary DAPI image
		selectImage("dapiDirectVessel");
		roiManager("Select", 2);
	    run("Clear Outside"); // isolating DAPI directly overlapping amyloid vessel
	    run("Select None");
	    run("Create Selection"); // make sure this does not need to be inverted
			// selection of DAPI directly overlapping vessel in grey matter
		run("Measure");
		print("Measured area of DAPI directly on vessel");
		
		// NEIGHBORING PLAQUE
		// open plaque distance map layer 100 image
		run("Bio-Formats Importer", "open=[" + distanceLayer100PlaqueDir + distancePlaqueFileList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	    plaque_distance = getTitle();
	    run("Bio-Formats Importer", "open=[" + distanceLayer100VesselDir + distanceVesselFileList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	    run("Create Selection");
	    roiManager("Add"); // vessel distance map selection (3)
	    selectImage(plaque_distance);
	    roiManager("Select", 3);
	    run("Clear");
	    	// vessel distance maps erased from plaque distance maps
	    rename("neighboringPlaque");
	    roiManager("Select", 0);
	    run("Clear Outside"); // isolating neighboring plaque area in gray matter
	    run("Select None");
	    roiManager("Select", 1);
	    run("Clear"); // removing binary plaque signal
	    run("Select None");
	    run("Create Selection"); // make sure this does not need to be inverted
			// selection of areas neighboring plaque (4)
		roiManager("Add");
		// measure value for neighboring plaque area
		run("Measure");
		print("Measured area of gray matter neighboring plaque");
		
		// select neighboringPlaque cleaned binary iba1 image
		selectImage("iba1NeighboringPlaque");
		roiManager("Select", 1);
		run("Clear"); // removing binary plaque signal
		run("Select None");
		roiManager("Select", 4); 
		run("Clear Outside"); // isolating neighboring plaque area in gray matter
		run("Select None");
		run("Create Selection"); // make sure this does not need to be inverted
			// selection of iba1 neighboring plaque in grey matter
		run("Measure");
		print("Measured area of Iba1 neighboring plaque");
		
		// select neighboringPlaque cleaned binary DAPI image
		selectImage("dapiNeighboringPlaque");
		roiManager("Select", 1);
		run("Clear"); // removing binary plaque signal
		run("Select None");
		roiManager("Select", 4);
		run("Clear Outside"); // isolating neighboring plaque area in gray matter
		run("Select None");
		run("Create Selection"); // make sure this does not need to be inverted
			// selection of DAPI neighboring plaque in grey matter
		run("Measure");
		print("Measured area of DAPI neighboring plaque");
		
		// NEIGHBORING VESSEL
		// open vessel distance map layer 100 image
		run("Bio-Formats Importer", "open=[" + distanceLayer100VesselDir + distanceVesselFileList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	    rename("neighboringVessel");
	    roiManager("Select", 0);
	    run("Clear Outside"); // isolating gray matter
	    run("Select None");
	    roiManager("Select", 2);
	    run("Clear"); // clearing binary vessel signal
	    run("Select None");
	    run("Create Selection"); // make sure this does not need to be inverted
			// selection of areas neighboring vessel (5)
		roiManager("Add");
		// measure value for neighboring vessel area
		run("Measure");
		print("Measured area of gray matter neighboring vessel");
		
		// select neighboringVessel cleaned binary iba1 image
		selectImage("iba1NeighboringVessel");
		roiManager("Select", 2);
		run("Clear"); // clearing binary vessel signal
		run("Select None");
		roiManager("Select", 5);
		run("Clear Outside"); // isolating neighboring vessel area in gray matter
		run("Select None");
		run("Create Selection"); // make sure this does not need to be inverted
			// selection of iba1 neighboring vessel in grey matter
		run("Measure");
		print("Measured area of Iba1 neighboring vessel");
		
		// select neighboringVessel cleaned binary DAPI image
		selectImage("dapiNeighboringVessel");
		roiManager("Select", 2);
		run("Clear"); // clearing binary vessel signal
		run("Select None");
		roiManager("Select", 5);
		run("Clear Outside"); // isolating neighboring vessel area in gray matter
		run("Select None");
		run("Create Selection"); // make sure this does not need to be inverted
			// selection of DAPI neighboring vessel in grey matter
		run("Measure");
		print("Measured area of DAPI neighboring vessel");
		
		// NOT NEIGHBORING PLAQUE
		// select notNeighboringPlaque cleaned binary iba1 image
		selectImage("iba1NotNeighboringPlaque");
		// area not neighboring plaque = total grey matter area - amyloid area - neighboring amyloid area
		roiManager("Select", 0);
		run("Clear Outside"); // isolating gray matter
		run("Select None");
		roiManager("Select", 1);
		run("Clear"); // removing binary plaque signal
		run("Select None");
		roiManager("Select", 4);
		run("Clear"); // removing neighboring plaque signal
		run("Select None");
		run("Create Selection"); // make sure this does not need to be inverted
			// selection of iba1 not neighboring plaque in grey matter
		run("Measure");
		print("Measured area of Iba1 not neighboring plaque");
		
		// select notNeighboringPlaque cleaned binary DAPI image
		selectImage("dapiNotNeighboringPlaque");
		// area not neighboring plaque = total grey matter area - plaque area - neighboring plaque area
		roiManager("Select", 0);
		run("Clear Outside"); // isolating gray matter
		run("Select None");
		roiManager("Select", 1);
		run("Clear"); // removing binary plaque signal
		run("Select None");
		roiManager("Select", 4);
		run("Clear"); // removing neighboring plaque signal
		run("Select None");
		run("Create Selection"); // make sure this does not need to be inverted
			// selection of DAPI not neighboring plaque in grey matter
		run("Measure");
		print("Measured area of DAPI not neighboring plaque");
		
		// NOT NEIGHBORING VESSEL
		// select notNeighboringVessel cleaned binary iba1 image
		selectImage("iba1NotNeighboringVessel");
		// area not neighboring vessel = total grey matter area - vessel area - neighboring vessel area
		roiManager("Select", 0);
		run("Clear Outside"); // isolating gray matter
		run("Select None");
		roiManager("Select", 2);
		run("Clear"); // removing binary vessel signal
		run("Select None");
		roiManager("Select", 5);
		run("Clear"); // removing neighboring vessel signal
		run("Select None");
		run("Create Selection"); // make sure this does not need to be inverted
			// selection of iba1 not neighboring vessel in grey matter
		run("Measure");
		print("Measured area of Iba1 not neighboring vessel");
		
		// select notNeighboringVessel cleaned binary DAPI image
		selectImage("dapiNotNeighboringVessel");
		// area not neighboring vessel = total grey matter area - vessel area - neighboring vessel area
		roiManager("Select", 0);
		run("Clear Outside"); // isolating gray matter
		run("Select None");
		roiManager("Select", 2);
		run("Clear"); // removing binary vessel signal
		run("Select None");
		roiManager("Select", 5);
		run("Clear"); // removing neighboring vessel signal
		run("Select None");
		run("Create Selection"); // make sure this does not need to be inverted
			// selection of DAPI not neighboring vessel in grey matter
		run("Measure");
		print("Measured area of DAPI neighboring vessel");
		
		// save all output to excel sheet
		run("Read and Write Excel", "no_count_column dataset_label=[" + sampleName + "] file_mode=queue_write");
		run("Read and Write Excel", "file_mode=write_and_close");
		print("Saved data in Excel");
		
		close("*");
		close("Results");
		close("ROI Manager");
		print("\n");
	}
}

// run it!
binaryIba1 = getDirectory("Navigate to folder containing the cleaned binary Iba1 images");
watershedDAPI = getDirectory("Navigate to folder containing the cleaned binary watershed DAPI images");
greyMatterROIs = getDirectory("Navigate to gray matter traces (saved as ROIs) ");
binaryPlaque = getDirectory("Navigate to folder containing the cleaned binary amyloid PLAQUE images");
binaryVessel = getDirectory("Navigate to folder containing the cleaned binary amyloid VESSEL images")
layer100PlaqueDistance = getDirectory("Navigate to folder containing layer 100 from the distance maps of the cleaned binary amyloid PLAQUE images");
layer100VesselDistance = getDirectory("Navigate to folder containing layer 100 from the distance maps of the cleaned binary amyloid VESSEL images");
outputDir = getDirectory("Navigate to the desired output folder ");
excelOutputName = getString("Enter the name you would like the Excel spreadsheet containing your output to have ", "default");
excelOutputPath = getDirectory("Navigate to where you would like the Excel spreadsheet containing your output to be stored ");

iba1DensitiesRegardingAmyloid(binaryIba1, greyMatterROIs, binaryPlaque, binaryVessel, layer100PlaqueDistance, layer100VesselDistance, watershedDAPI, outputDir, excelOutputName, excelOutputPath);
print("Script Finished");