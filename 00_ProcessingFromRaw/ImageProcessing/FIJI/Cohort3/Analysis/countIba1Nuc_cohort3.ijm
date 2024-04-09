// determining Iba1+ nuclei for cohort 3

function countIba1Nuc(iba1Dir, dapiWatershedDir, greyMatterRoiDir, binaryAmyloidDir, distanceLayer100AmyloidDir, output, excelName, excelPath) {
	print("Function running...");
	iba1FileList = getFileList(iba1Dir);
	dapiFileList = getFileList(dapiWatershedDir);
	binaryAmyloidFileList = getFileList(binaryAmyloidDir);
	distanceAmyloidFileList = getFileList(distanceLayer100AmyloidDir);
	roiFileList = getFileList(greyMatterRoiDir);
	
	for (i = 0; i < iba1FileList.length; i++) {
	    // open Iba1 cleaned binary image
	    run("Bio-Formats Importer", "open=[" + iba1Dir + iba1FileList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	    run("Invert LUT");
	    current = getTitle();
	    numCharsBeforeSpace = 0;
		for (c = current.length - 1; c > 0; c--) {
		    if (fromCharCode(charCodeAt(current, c)) == "/") { // may need to be adjusted depending on naming convention
		        break;
		    }
		    numCharsBeforeSpace++;
		}
		sampleName = current.substring(current.length - numCharsBeforeSpace, current.length - 4);
	    print("Processing " + sampleName + "...");
	    
	    // prepare excel sheel for data output
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
		run("Clear Outside"); // isolating gray matter signal
		run("Select None");
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
		roiManager("Select", 0); // isolating gray matter signal
		run("Clear Outside");
		run("Select None");
		
		// duplicate three times
		dapiMain = getTitle();
		run("Duplicate...", "title=[dapi_count_directAmyloid]");
		run("Duplicate...", "title=[dapi_count_neighboringAmyloid]");
		run("Duplicate...", "title=[dapi_count_notNeighboringAmyloid]");
		print("Duplicated DAPI image");
		
		// also count DAPI cells here
		selectImage(dapiMain);
		rename("dapi_count");
		run("Find Maxima...", "prominence=10 exclude output=[Count]");
			// prominence doesn't matter bc images are binary
		run("Read and Write Excel", "no_count_column dataset_label=[Total_DAPI_Count] file_mode=queue_write");
		print("Wrote baseline gray matter DAPI count");
		
		// restore Iba1+ cell selection
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
		run("Duplicate...", "title=[iba1Nuc_count_directAmyloid]");
		run("Duplicate...", "title=[iba1Nuc_count_neighboringAmyloid]");
		run("Duplicate...", "title=[iba1Nuc_count_notNeighboringAmyloid]");
		print("Duplicated Iba1+ nuclei image");
		
		// count Iba1+ nuclei
		selectImage(isolatedIba1Nuc);
		rename("iba1Nuc_count_gray_matter");
		run("Find Maxima...", "prominence=10 exclude output=[Count]");
		run("Read and Write Excel", "no_count_column stack_results dataset_label=[Iba1+_Nucleus_Count] file_mode=queue_write");
		print("Wrote Iba1+ nuclei count");
		
		// determine Iba1+ nuclei overlapping amyloid
		run("Bio-Formats Importer", "open=[" + binaryAmyloidDir + binaryAmyloidFileList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	    run("Invert LUT");
	    rename("binaryPlaque");
	    roiManager("Select", 0);
	    run("Clear Outside"); // isolating gray matter signal
	    run("Select None");
	    run("Create Selection"); // make sure this does not need to be inverted
			// selection of amyloid present in grey matter (2)
		roiManager("Add");
		
		// counts overlapping amyloid
		selectImage("iba1Nuc_count_directAmyloid");
		roiManager("Select", 2);
		run("Clear Outside"); // isolating overlapping binary amyloid region
		run("Find Maxima...", "prominence=10 exclude output=[Count]");
		run("Read and Write Excel", "no_count_column stack_results dataset_label=[Iba1+_Nucleus_Count_Direct_Amyloid] file_mode=queue_write");
		print("Wrote Iba1+ nuclei overlapping amyloid count");
		
		selectImage("dapi_count_directAmyloid");
		roiManager("Select", 2);
		run("Clear Outside"); // isolating overlapping binary amyloid region
		run("Find Maxima...", "prominence=10 exclude output=[Count]");
		run("Read and Write Excel", "no_count_column stack_results dataset_label=[DAPI_Count_Direct_Amyloid] file_mode=queue_write");
		print("Wrote DAPI overlapping amyloid count");
		
		// determine Iba1+ nuclei neighboring amyloid
		run("Bio-Formats Importer", "open=[" + distanceLayer100AmyloidDir + distanceAmyloidFileList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	    rename("neighboringAmyloid");
	    roiManager("Select", 0);
	    run("Clear Outside");
	    run("Select None");
	    roiManager("Select", 2);
	    run("Clear"); // removing overlapping binary amyloid region
	    run("Select None");
	    run("Create Selection"); // make sure this does not need to be inverted
			// selection of neighboring amyloid area in grey matter (3)
		roiManager("Add");
	    
		// counts neighboring amyloid
		selectImage("iba1Nuc_count_neighboringAmyloid");
		roiManager("Select", 3);
		run("Clear Outside"); // isolating neughoring amyloid area
		run("Find Maxima...", "prominence=10 exclude output=[Count]");
		run("Read and Write Excel", "no_count_column stack_results dataset_label=[Iba1+_Nucleus_Count_Neighbor_Amyloid] file_mode=queue_write");
		print("Wrote Iba1+ nuclei neighboring amyloid count");
		
		selectImage("dapi_count_neighboringAmyloid");
		roiManager("Select", 3);
		run("Clear Outside"); // isolating neughoring amyloid area
		run("Find Maxima...", "prominence=10 exclude output=[Count]");
		run("Read and Write Excel", "no_count_column stack_results dataset_label=[Iba1+_Nucleus_Count_Neighbor_Amyloid] file_mode=queue_write");
		print("Wrote DAPI neighboring amyloid count");
		
		
		// determine Iba1+ nuclei not neighboring amyloid
		selectImage("iba1Nuc_count_notNeighboringAmyloid");
		// area not neighboring amyloid = total grey matter area - amyloid area - neighboring amyloid area
		roiManager("Select", 0);
		run("Clear Outside");
		run("Select None");
		roiManager("Select", 2);
		run("Clear"); // removing overlapping binary amyloid region
		run("Select None");
		roiManager("Select", 3);
		run("Clear"); // removing neughoring amyloid area
		run("Select None");
		run("Find Maxima...", "prominence=10 exclude output=[Count]");
		run("Read and Write Excel", "no_count_column stack_results dataset_label=[Iba1+_Nucleus_Count_Not_Neighbor_Amyloid] file_mode=queue_write");
		print("Wrote Iba1+ nuclei not neighboring amyloid count");
		
		selectImage("dapi_count_notNeighboringAmyloid");
		// area not neighboring amyloid = total grey matter area - amyloid area - neighboring amyloid area
		roiManager("Select", 0);
		run("Clear Outside");
		run("Select None");
		roiManager("Select", 2);
		run("Clear"); // removing overlapping binary amyloid region
		run("Select None"); 
		roiManager("Select", 3);
		run("Clear"); // removing neughoring amyloid area
		run("Select None");
		run("Find Maxima...", "prominence=10 exclude output=[Count]");
		run("Read and Write Excel", "no_count_column stack_results dataset_label=[Iba1+_Nucleus_Count_Not_Neighbor_Amyloid] file_mode=queue_write");
		print("Wrote DAPI not neighboring amyloid count");
		
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
binaryAmyloid = getDirectory("Navigate to folder containing the cleaned binary amyloid images ");
layer100AmyloidDistance = getDirectory("Navigate to folder containing layer 100 of the distance maps for the cleaned binary amyloid images ");
outputDirectory = getDirectory("Navigate to the desired output folder ");
excelOutputName = getString("Enter the name you would like the Excel spreadsheet containing your output to have ", "default");
excelOutputPath = getDirectory("Navigate to where you would like the Excel spreadsheet containing your output to be stored ");

countIba1Nuc(binaryIba1, watershedDAPI, greyMatterROIs, binaryAmyloid, layer100AmyloidDistance, outputDirectory, excelOutputName, excelOutputPath);
print("Script Finished");