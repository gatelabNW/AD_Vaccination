// Makes Z-projections using Maximum Intensity for microglia-neuron interaction analysis

// Assumes multiple images per sample
// Output stores all images per sample in same folder, which is created in the script

function maxProjections(originalDir, outputDir) {
	print("Function running...");
	all_images = getFileList(originalDir);
	
	for (i = 0; i < all_images.length; i++) {
		folder = all_images[i];
		sampleName = folder.substring(0, folder.length-1);
		imagesPath = originalDir + all_images[i] + "20x/";
		indiv_images = getFileList(imagesPath);
		print("Processing " + sampleName);
		
		sampleOutput = outputDir + sampleName + "/";
		File.makeDirectory(sampleOutput);
		print("Made output folder");
		
		for (j = 0; j < indiv_images.length; j++) {
			current = indiv_images[j];
			run("Bio-Formats Importer", "open=[" + imagesPath + indiv_images[j] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
			run("Z Project...", "projection=[Max Intensity]");
			saveAs("tiff", sampleOutput + sampleName + "_20x_image" + j + "_maxProject");
			print("Made max projection of " + current);
			close("*");
		}
		print("Done with " + sampleName);
	}
}

// run it!

originalFiles = getDirectory("Navigate to the original images you wish to make max projections of ");
outputFolder = getDirectory("Navigate to the desired output folder ");

maxProjections(originalFiles, outputFolder);

print("Script finished");