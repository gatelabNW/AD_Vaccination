// determining Iba1+ and DAPI+ density in regions overlapping amyloid, 
// neighboring amyloid, and not neighboring amyloid
// Cohort 3

function iba1DensitiesRegardingAmyloid(iba1Dir, greyMatterRoiDir, binaryAmyloidDir, distanceLayer100AmyloidDir, dapiDir, output, excelFileName, excelOutput) {
	print("Function running...");
	iba1FileList = getFileList(iba1Dir);
	greyMatterRoiFileList = getFileList(greyMatterRoiDir);
	amyloidFileList = getFileList(binaryAmyloidDir);
	distanceAmyloidFileList = getFileList(distanceLayer100AmyloidDir);
	dapiFileList = getFileList(dapiDir);
		
	// also measure DAPI coverage in all of these regions
	
	for (i = 0; i < iba1FileList.length; i++) {
		// open cleaned binary iba1
		run("Bio-Formats Importer", "open=[" + iba1Dir + iba1FileList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	    currentIba1 = getTitle();
	    sample = iba1FileList[i];
	    sampleName = sample.substring(0, 2); // may need to be adjusted
	    print("Processing " + sampleName + "...");
	
		// duplicate three times
		run("Duplicate...", "title=[iba1DirectAmyloid]");
		run("Duplicate...", "title=[iba1NeighboringAmyloid]");
		run("Duplicate...", "title=[iba1NotNeighboringAmyloid]");
		
		print("Opened and duplicated cleaned binary Iba1 image");
		
		// open cleaned binary DAPI 
		run("Bio-Formats Importer", "open=[" + dapiDir + dapiFileList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	    run("Invert");
	    run("Invert LUT");
	    currentDAPI = getTitle();
	    
	    // duplicate three times		
		run("Duplicate...", "title=[dapiDirectAmyloid]");
		run("Duplicate...", "title=[dapiNeighboringAmyloid]");
		run("Duplicate...", "title=[dapiNotNeighboringAmyloid]");
		
		print("Opened and duplicated cleaned binary DAPI image");
		
		// prepare excel file
		run("Read and Write Excel", "no_count_column file=[" + excelOutput + excelFileName + "] sheet=" + sampleName + "");
		run("Read and Write Excel", "file=[" + excelOutput + excelFileName + " file_mode=read_and_open");
		
		// open grey matter ROI + add to manager (0)
		open(greyMatterRoiDir + greyMatterRoiFileList[i]);
		roiManager("Add");
		run("Select None");
		selectImage(currentIba1);
		rename("Grey Matter");
		roiManager("Select", 0);
		run("Clear Outside"); // isolating gray matter area
		// measure value of total grey matter area
		run("Measure");
		print("Measured whole grey matter area");
		saveAs("tiff", output + sampleName + "_greyMatterSelection.tif");
		
		// measure total Iba1 and DAPI coverage in grey matter
		rename("Iba1GreyMatter");
		run("Create Selection");
		run("Measure");
		print("Measured total area of Iba1 in grey matter");
		run("Select None");
		close();
		
		selectImage(currentDAPI);
		rename("DAPIGreyMatter");
		roiManager("Select", 0);
		run("Clear Outside"); // isolating gray matter area
		run("Select None");
		run("Create Selection");
		run("Measure");
		print("Measured total area of DAPI in grey matter");
		run("Select None");
		close();
		
		// DIRECTLY OVERLAPPING
		// open cleaned binary amyloid
		run("Bio-Formats Importer", "open=[" + binaryAmyloidDir + amyloidFileList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	    rename("binaryAmyloid");
	    roiManager("Select", 0);
	    run("Clear Outside"); // isolating gray matter area
	    run("Select None");
	    run("Create Selection"); // make sure this does not need to be inverted
			// selection of amyloid present in grey matter
		roiManager("Add");
		// measure value for total amyloid area (1)
		run("Measure");
		print("Measured total area of amyloid in grey matter");
		
		// NEIGHBORING AMYLOID
		// open plaque distance map layer 100 image
		run("Bio-Formats Importer", "open=[" + distanceLayer100AmyloidDir + distanceAmyloidFileList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	    amyloid_distance = getTitle();
	    rename("neighboringAmyloid");
	    roiManager("Select", 0);
	    run("Clear Outside"); // isolating gray matter area
	    run("Select None");
	    roiManager("Select", 1); // clearing binary amyloid signal
	    run("Clear");
	    run("Select None");
	    run("Create Selection"); // make sure this does not need to be inverted
			// selection of areas neighboring amyloid (2)
		roiManager("Add");
		// measure value for neighboring amyloid area
		run("Measure");
		print("Measured area of gray matter neighboring amyloid");
		
		// select neighboringAmyloid cleaned binary iba1 image
		selectImage("iba1NeighboringAmyloid");
		roiManager("Select", 1);
		run("Clear"); // clearing binary amyloid signal
		run("Select None");
		roiManager("Select", 2);
		run("Clear Outside"); // isolating neighboring amyloid signal
		run("Select None");
		run("Create Selection"); // make sure this does not need to be inverted
			// selection of iba1 neighboring amyloid in grey matter
		run("Measure");
		print("Measured area of Iba1 neighboring amyloid");
		
		// select neighboringAmyloid cleaned binary DAPI image
		selectImage("dapiNeighboringAmyloid");
		roiManager("Select", 1);
		run("Clear"); // clearing binary amyloid signal
		run("Select None");
		roiManager("Select", 2);
		run("Clear Outside"); // isolating neighboring amyloid signal
		run("Select None");
		run("Create Selection"); // make sure this does not need to be inverted
			// selection of DAPI neighboring amyloid in grey matter
		run("Measure");
		print("Measured area of DAPI neighboring amyloid");
	
		// NOT NEIGHBORING AMYLOID
		// select notNeighboringAmyloid cleaned binary iba1 image
		selectImage("iba1NotNeighboringAmyloid");
		// area not neighboring amyloid = total grey matter area - amyloid area - neighboring amyloid area
		roiManager("Select", 0);
		run("Clear Outside"); // isolating gray matter area
		run("Select None");
		roiManager("Select", 1);
		run("Clear"); // clearing binary amyloid signal
		run("Select None");
		roiManager("Select", 2);
		run("Clear"); // clearing neighboring amyloid signal
		run("Select None");
		run("Create Selection"); // make sure this does not need to be inverted
			// selection of iba1 not neighboring amyloid in grey matter
		run("Measure");
		print("Measured area of Iba1 not neighboring amyloid");
		
		// select notNeighboringAmyloid cleaned binary DAPI image
		selectImage("dapiNotNeighboringAmyloid");
		// area not neighboring amyloid = total grey matter area - plaque area - neighboring plaque area
		roiManager("Select", 0);
		run("Clear Outside"); // isolating gray matter area
		run("Select None");
		roiManager("Select", 1);
		run("Clear"); // clearing binary amyloid signal
		run("Select None");
		roiManager("Select", 2);
		run("Clear"); // clearing neighboring amyloid signal
		run("Select None");
		run("Create Selection"); // make sure this does not need to be inverted
			// selection of DAPI not neighboring amyloid in grey matter
		run("Measure");
		print("Measured area of DAPI not neighboring amyloid");
		
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
binaryAmyloid = getDirectory("Navigate to folder containing the cleaned binary amyloid images");
layer100AmyloidDistance = getDirectory("Navigate to folder containing layer 100 from the distance maps of the cleaned binary amyloid images");
outputDir = getDirectory("Navigate to the desired output folder ");
excelOutputName = getString("Enter the name you would like the Excel spreadsheet containing your output to have ", "default");
excelOutputPath = getDirectory("Navigate to where you would like the Excel spreadsheet containing your output to be stored ");

iba1DensitiesRegardingAmyloid(binaryIba1, greyMatterROIs, binaryAmyloid, layer100AmyloidDistance, watershedDAPI, outputDir, excelOutputName, excelOutputPath);
print("Script Finished");