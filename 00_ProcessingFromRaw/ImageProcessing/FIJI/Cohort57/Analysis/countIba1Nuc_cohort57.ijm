// determining Iba1+ nuclei for cohorts 5/7

function countIba1Nuc(iba1Dir, dapiWatershedDir, greyMatterRoiDir, binaryPlaqueDir, binaryVesselDir, distanceLayer100PlaqueDir, distanceLayer100VesselDir, output, excelName, excelPath) {
	print("Function running...");
	iba1FileList = getFileList(iba1Dir);
	dapiFileList = getFileList(dapiWatershedDir);
	binaryPlaqueFileList = getFileList(binaryPlaqueDir);
	binaryVesselFileList = getFileList(binaryVesselDir);
	distancePlaqueFileList = getFileList(distanceLayer100PlaqueDir);
	distanceVesselFileList = getFileList(distanceLayer100VesselDir);
	roiFileList = getFileList(greyMatterRoiDir);
	
	for (i = 0; i < iba1FileList.length; i++) {
		roiCounter = 0;
	    // open Iba1 cleaned binary image
	    run("Bio-Formats Importer", "open=[" + iba1Dir + iba1FileList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	    current = getTitle();
	    sample = iba1FileList[i];
	    sampleName = sample.substring(0, 2); // // may need to be adjusted depending on naming convention
	    print("Processing " + sampleName + "...");
	    // isolating neughoring amyloid area
	    // prepare excel sheel for data output
	    run("Read and Write Excel", "file_mode=write_and_close");
		run("Read and Write Excel", "file=[" + excelPath + excelName + "] sheet=" + sampleName);
		run("Read and Write Excel", "file=[" + excelPath + excelName + "] file_mode=read_and_open");
	
	    // process Iba1 image
	    run("Erode");
	    run("Erode");
	    print("Erosions done");
	    
	    // isolate Iba1 signal in the gray matter
	    open(greyMatterRoiDir + roiFileList[i]);
		roiManager("Add"); // 0 gray matter signal
		run("Select None");
		selectImage(current);
		roiManager("Select", 0);
		run("Clear Outside");
		run("Select None");
		roiCounter++;
		print("Gray matter isolated");
		
		// create selection of remaining Iba1 cells
	    run("Create Selection"); 
		roiManager("Add"); // 1 Iba1 cells in gray matter
		print("Selection made");
	    
		// open DAPI cleaned WATERSHED image
		run("Bio-Formats Importer", "open=[" + dapiWatershedDir + dapiFileList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		run("Invert");
		run("Invert LUT");
		run("Fill Holes");
		print("Ran fill holes");
		roiManager("Select", 0); // isolated dapi signal in gray matter
		run("Clear Outside");
		run("Select None");
		
		// duplicate three times
		dapiMain = getTitle();
		run("Duplicate...", "title=[dapi_count_directPlaque]");
		run("Duplicate...", "title=[dapi_count_neighboringPlaque]");
		run("Duplicate...", "title=[dapi_count_notNeighboringPlaque]");
		run("Duplicate...", "title=[dapi_count_directVessel]");
		run("Duplicate...", "title=[dapi_count_neighboringVessel]");
		run("Duplicate...", "title=[dapi_count_notNeighboringVessel]");
		print("Duplicated DAPI image");
		
		// also count DAPI cells here
		selectImage(dapiMain);
		rename("dapi_count");
		run("Find Maxima...", "prominence=10 exclude output=[Count]");
			// prominence doesn't matter bc images are binary
		run("Read and Write Excel", "no_count_column dataset_label=[Total_DAPI_Count] file_mode=queue_write");
		print("Wrote baseline gray matter DAPI count");
		
		// restore Iba1+ cells selection
		roiManager("Select", 1);
		run("Clear Outside"); // isolate Iba1+DAPI+ cells in gray matter
		run("Select None");
		print("Isolated Iba1+ nuclei");
		
		run("Analyze Particles...", "size=10-infinity pixel show=Masks clear");
		run("Invert LUT");
		saveAs("tiff", output + sampleName + "_Iba1Nuc.tif");
		print("Analyzed particles and saved mask");
		
		// duplicate three times
		isolatedIba1Nuc = getTitle();
		run("Duplicate...", "title=[iba1Nuc_count_directPlaque]");
		run("Duplicate...", "title=[iba1Nuc_count_neighboringPlaque]");
		run("Duplicate...", "title=[iba1Nuc_count_notNeighboringPlaque]");
		run("Duplicate...", "title=[iba1Nuc_count_directVessel]");
		run("Duplicate...", "title=[iba1Nuc_count_neighboringVessel]");
		run("Duplicate...", "title=[iba1Nuc_count_notNeighboringVessel]");
		print("Duplicated Iba1+ nuclei image");
		
		// count Iba1+ nuclei
		selectImage(isolatedIba1Nuc);
		rename("iba1Nuc_count_gray_matter");
		run("Find Maxima...", "prominence=10 exclude output=[Count]");
		run("Read and Write Excel", "no_count_column stack_results dataset_label=[Iba1+_Nucleus_Count] file_mode=queue_write");
		print("Wrote Iba1+ nuclei count");
		
		// determine Iba1+ nuclei overlapping amyloid PLAQUES
		run("Bio-Formats Importer", "open=[" + binaryPlaqueDir + binaryPlaqueFileList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	    rename("binaryPlaque");
	    roiManager("Select", 0);
	    run("Clear Outside"); // isolated gray matter
	    run("Select None");
	    run("Create Selection"); // make sure this does not need to be inverted
			// selection of amyloid plaques present in grey matter (2)
		roiManager("Add");
		
		// counts overlapping amyloid PLAQUES
		selectImage("iba1Nuc_count_directPlaque");
		roiManager("Select", 2); 
		run("Clear Outside"); // isolating binary amyloid plaque in gray matter
		run("Find Maxima...", "prominence=10 exclude output=[Count]");
		run("Read and Write Excel", "no_count_column stack_results dataset_label=[Iba1+_Nucleus_Count_Direct_Plaque] file_mode=queue_write");
		print("Wrote Iba1+ nuclei overlapping amyloid plaque count");
		
		selectImage("dapi_count_directPlaque");
		roiManager("Select", 2);
		run("Clear Outside"); // isolating binary amyloid plaque in gray matter
		run("Find Maxima...", "prominence=10 exclude output=[Count]");
		run("Read and Write Excel", "no_count_column stack_results dataset_label=[DAPI_Count_Direct_Plaque] file_mode=queue_write");
		print("Wrote DAPI overlapping amyloid plaque count");
		
		// determine Iba1+ nuclei overlapping amyloid VESSELS
		run("Bio-Formats Importer", "open=[" + binaryVesselDir + binaryVesselFileList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	    rename("binaryVessel");
	    roiManager("Select", 0);
	    run("Clear Outside"); // isolated gray matter
	    run("Select None");
	    run("Create Selection"); // make sure this does not need to be inverted
			// selection of amyloid vessels present in grey matter (3)
		roiManager("Add");
		
		// counts overlapping amyloid VESSELS
		selectImage("iba1Nuc_count_directVessel");
		roiManager("Select", 3);
		run("Clear Outside"); // isolating binary amyloid vessels in gray matter
		run("Find Maxima...", "prominence=10 exclude output=[Count]");
		run("Read and Write Excel", "no_count_column stack_results dataset_label=[Iba1+_Nucleus_Count_Direct_Vessel] file_mode=queue_write");
		print("Wrote Iba1+ nuclei overlapping amyloid vessel count");
		
		selectImage("dapi_count_directVessel");
		roiManager("Select", 3);
		run("Clear Outside"); // isolating binary amyloid vessels in gray matter
		run("Find Maxima...", "prominence=10 exclude output=[Count]");
		run("Read and Write Excel", "no_count_column stack_results dataset_label=[DAPI_Count_Direct_Vessel] file_mode=queue_write");
		print("Wrote DAPI overlapping amyloid vessel count");
		
		// determine Iba1+ nuclei neighboring amyloid VESSELS
		run("Bio-Formats Importer", "open=[" + distanceLayer100VesselDir + distanceVesselFileList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	    rename("neighboringVessel");
	    roiManager("Select", 0);
	    run("Clear Outside"); // isolated gray matter
	    run("Select None");
	    roiManager("Select", 3);
	    run("Clear"); // removing binary amyloid vessels in gray matter
	    run("Select None");
	    run("Create Selection"); // make sure this does not need to be inverted
			// selection of areas neighboring amyloid vessels (4)
		roiManager("Add");
		
		// counts neighboring amyloid VESSELS
		selectImage("iba1Nuc_count_neighboringVessel");
		roiManager("Select", 4);
		run("Clear Outside"); // isolating neighboring vessel area
		run("Find Maxima...", "prominence=10 exclude output=[Count]");
		run("Read and Write Excel", "no_count_column stack_results dataset_label=[Iba1+_Nucleus_Count_Neighbor_Vessel] file_mode=queue_write");
		print("Wrote Iba1+ nuclei neighboring amyloid vessel count");
		
		selectImage("dapi_count_neighboringVessel");
		roiManager("Select", 4);
		run("Clear Outside"); // isolating neighboring vessel area
		run("Find Maxima...", "prominence=10 exclude output=[Count]");
		run("Read and Write Excel", "no_count_column stack_results dataset_label=[Iba1+_Nucleus_Count_Neighbor_Vessel] file_mode=queue_write");
		print("Wrote DAPI neighboring amyloid vessel count");
		
		// determine Iba1+ nuclei neighboring amyloid PLAQUES
		run("Bio-Formats Importer", "open=[" + distanceLayer100PlaqueDir + distancePlaqueFileList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	    rename("neighboringPlaque");
	    roiManager("Select", 0);
	    run("Clear Outside"); // isolated gray matter
	    run("Select None");
	    roiManager("Select", 2);
	    run("Clear"); // removing binary amyloid plaque in gray matter
	    run("Select None");
	    roiManager("Select", 4);
	    run("Clear"); // removing neighboring vessel area
	    run("Select None");
	    run("Create Selection"); // make sure this does not need to be inverted
			// selection of areas neighboring amyloid plaque (5)
		roiManager("Add");
		run("Select None");
		saveAs("tiff", "//fsmresfiles.fsm.northwestern.edu/fsmresfiles/Neurology/Gate_Lab/projects/AN1792-vacc/04 raw-data/07 Image analysis/cohort_5-7/amyloid/distance_maps/plaques_layer100_minus_vessels_layer100/" + sampleName + "_plaqueLayer100WithoutVesselDM");
		
		// counts neighboring amyloid PLAQUES
		selectImage("iba1Nuc_count_neighboringPlaque");
		roiManager("Select", 5);
		run("Clear Outside"); // isolating neighboring plaque area
		run("Find Maxima...", "prominence=10 exclude output=[Count]");
		run("Read and Write Excel", "no_count_column stack_results dataset_label=[Iba1+_Nucleus_Count_Neighbor_Plaque] file_mode=queue_write");
		print("Wrote Iba1+ nuclei neighboring amyloid plaque count");
		
		selectImage("dapi_count_neighboringPlaque");
		roiManager("Select", 5);
		run("Clear Outside"); // isolating neighboring plaque area
		run("Find Maxima...", "prominence=10 exclude output=[Count]");
		run("Read and Write Excel", "no_count_column stack_results dataset_label=[Iba1+_Nucleus_Count_Neighbor_Plaque] file_mode=queue_write");
		print("Wrote DAPI neighboring amyloid plaque count");
		
		
		// determine Iba1+ nuclei not neighboring amyloid PLAQUE
		selectImage("iba1Nuc_count_notNeighboringPlaque");
		// area not neighboring amyloid = total grey matter area - amyloid area - neighboring amyloid area
		roiManager("Select", 0);
		run("Clear Outside"); // isolated gray matter
		run("Select None");
		roiManager("Select", 2);
		run("Clear"); // isolating binary amyloid plaque in gray matter
		run("Select None");
		roiManager("Select", 5);
		run("Clear"); // removing neighboring plaque area
		run("Select None");
		run("Find Maxima...", "prominence=10 exclude output=[Count]");
		run("Read and Write Excel", "no_count_column stack_results dataset_label=[Iba1+_Nucleus_Count_Not_Neighbor_Plaque] file_mode=queue_write");
		print("Wrote Iba1+ nuclei not neighboring amyloid plaque count");
		
		// determine Iba1+ nuclei not neighboring amyloid VESSEL
		selectImage("iba1Nuc_count_notNeighboringPlaque");
		// area not neighboring amyloid = total grey matter area - amyloid area - neighboring amyloid area
		roiManager("Select", 0);
		run("Clear Outside"); // isolated gray matter
		run("Select None");
		roiManager("Select", 3);
		run("Clear"); // isolating binary amyloid vessels in gray matter
		run("Select None");
		roiManager("Select", 4);
		run("Clear"); // removing neighboring vessel area
		run("Select None");
		run("Find Maxima...", "prominence=10 exclude output=[Count]");
		run("Read and Write Excel", "no_count_column stack_results dataset_label=[Iba1+_Nucleus_Count_Not_Neighbor_Vessel] file_mode=queue_write");
		print("Wrote Iba1+ nuclei not neighboring amyloid vessel count");
		
		selectImage("dapi_count_notNeighboringPlaque");
		// area not neighboring amyloid = total grey matter area - amyloid area - neighboring amyloid area
		roiManager("Select", 0);
		run("Clear Outside"); // isolated gray matter
		run("Select None");
		roiManager("Select", 2);
		run("Clear"); // isolating binary amyloid plaque in gray matter
		run("Select None");
		roiManager("Select", 5);
		run("Clear"); // removing neighboring plaque area
		run("Select None");
		run("Find Maxima...", "prominence=10 exclude output=[Count]");
		run("Read and Write Excel", "no_count_column stack_results dataset_label=[Iba1+_Nucleus_Count_Not_Neighbor_Plaque] file_mode=queue_write");
		print("Wrote DAPI not neighboring amyloid plaque count");
		
		selectImage("dapi_count_notNeighboringVessel");
		// area not neighboring amyloid = total grey matter area - amyloid area - neighboring amyloid area
		roiManager("Select", 0);
		run("Clear Outside"); // isolated gray matter
		run("Select None");
		roiManager("Select", 3);
		run("Clear"); // isolating binary amyloid vessels in gray matter
		run("Select None");
		roiManager("Select", 4);
		run("Clear"); // removing neighboring vessel area
		run("Select None");
		run("Find Maxima...", "prominence=10 exclude output=[Count]");
		run("Read and Write Excel", "no_count_column stack_results dataset_label=[Iba1+_Nucleus_Count_Not_Neighbor_Vessel] file_mode=queue_write");
		print("Wrote DAPI not neighboring amyloid vessel count");
		
		// save all output to excel sheet
		run("Read and Write Excel", "file_mode=write_and_close");
		print("Saved data in Excel\n");
		
		close("*");
		close("Results");
		close("ROI Manager");
	}
}

// run it!
binaryIba1 = getDirectory("Navigate to the folder containing the cleaned binary Iba1 images ");
watershedDAPI = getDirectory("Navigate to folder containing the cleaned binary watershed DAPI images ");
greyMatterROIs = getDirectory("Navigate to folder containing the gray matter traces (saved as ROIs) ");
binaryPlaque = getDirectory("Navigate to folder containing the cleaned binary amyloid PLAQUE images ");
binaryVessel = getDirectory("Navigate to folder containing the cleaned binary amyloid VESSEL images ");
layer100PlaqueDistance = getDirectory("Navigate to folder containing layer 100 of the distance maps for the cleaned binary amyloid PLAQUE images ");
layer100VesselDistance = getDirectory("Navigate to folder containing layer 100 of the distance maps for the cleaned binary amyloid VESSEL images ");
outputDirectory = getDirectory("Navigate to the desired output folder ");
excelOutputName = getString("Enter the name you would like the Excel spreadsheet containing your output to have ", "default");
excelOutputPath = getDirectory("Navigate to where you would like the Excel spreadsheet containing your output to be stored ");

countIba1Nuc(binaryIba1, watershedDAPI, greyMatterROIs, binaryPlaque, binaryVessel, layer100PlaqueDistance, layer100VesselDistance, outputDirectory, excelOutputName, excelOutputPath);
print("Script Finished");