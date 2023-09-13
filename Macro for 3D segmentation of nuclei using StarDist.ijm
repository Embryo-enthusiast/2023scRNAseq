// Alberio Lab image analysis macro for 3D objects

// Before you start you need to install Fiji from this link: https://imagej.net/software/fiji/
//Then go to Help -> Update... -> Manage update sites -> Select 3D ImageJ Suite, Bio-Formats, CSBDeep, StarDist, TensorFlow, TrackMate-StarDist -> Close -> Apply changes


// Base_dir is the path to the directory in which your files are held
// image_folder is the path to the directory in which your TIF image files are held
// DAPI_channel_number=the order of the DAPI channel in a multichannel image e.g. C1 (channel 1), C2 (channel 2)
// Number_of_channels=the total number of channels in the image (currently written for a 4 channel image)
// This Macro assumes that the /in/, /out/ and /csv/ directories are all within the base dir (ex: C/USER/SEGMENTATION/PROJECT/in,  C/USER/SEGMENTATION/PROJECT/out and C/USER/SEGMENTATION/PROJECT/csv)
//this pipe line will out put CSVs for each object and each channel. These can be further processed using the R code entitled "Downstream 3d Segmentation.R"
//Files should be .czi on input  

base_dir="USER_INPUT_HERE";
image_folder=base_dir + "/in/";
segmentation_output_folder=base_dir + "/out/";
csvs_folder=base_dir + "/csv/";

DAPI_channel_number="C1"

Number_of_channels=4



///////////////// NO USER INPUT BEYOND HERE ///////////////////////////////////////

image_list=getFileList(image_folder);

for (i=0; i < image_list.length; i++) {
	image_name=image_list[i];
	image_path = image_folder + image_name;
	run("Bio-Formats Importer", "open=["+image_path+"] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	rename(image_name);
	image_title=getTitle();
	//image_title = File.nameWithoutExtension;
	run("Split Channels");

	selectWindow("C1-" + image_title);
	run("Grays");
	//saveAs("Tiff", segmentation_output_folder + "C1-" + image_title + ".tif");
	rename("C1-" + image_title + ".tif");
	
	selectWindow("C2-" + image_title);
	run("Green");
	//saveAs("Tiff", segmentation_output_folder + "C2-" + image_title + ".tif");
	rename("C2-" + image_title + ".tif");

	selectWindow("C3-" + image_title);
	run("Cyan");
	//saveAs("Tiff", segmentation_output_folder + "C3-" + image_title + ".tif");
	rename("C3-" + image_title + ".tif");

	selectWindow("C4-" + image_title);
	run("Magenta");
	//saveAs("Tiff", segmentation_output_folder + "C4-" + image_title + ".tif");
	rename("C4-" + image_title + ".tif");

	
	selectWindow(DAPI_channel_number + "-" + image_title + ".tif");
	run("Duplicate...", "duplicate");
	run("Unsharp Mask...", "radius=15 mask=0.80 stack");
	rename("Image");
	
	run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
	run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'Image', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'99.8', 'probThresh':'0.479071', 'nmsThresh':'0.3', 'outputType':'Both', 'nTiles':'1', 'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
	
	run("TrackMate");
	waitForUser("Manual running of the Trackmate-StarDist pipeline", "DO NOT PRESS OKAY UNTIL MANUAL RUNNING OF THE TRACKMATE-STARDIST PIPELINE IS COMPLETE \n Follow these selection instructions: \n 1. Click Next, there is no need to adjust the crop settings \n 2. Select 'Label image detector' then click Next \n 3. Click 'Preview' to see how well the 2D StarDist detector is performing on an individual slice and when satisifed click Next \n 4. Wait for the Detection bar at the top to complete and then click Next (this may take a few minutes) \n 5. If the number of spots is > the number of slices* the approximate number of cells then consider adjusting the Quality bar to reduce the number of spots (or objects) \n taken forward, otherwise click Next \n 6. Wait for the calculation to complete then add filters to the spots using the green plus button e.g. filters for Area/Average Ch1 (DAPI) intensity can help to filter out erroneous segmentation \n (this can be checked by seeing which spots disappear from the visible slice), then click Next \n 7. Select the LAP Tracker then click Next \n 8. Set the Frame-to-frame linking to 4uM, tick the Allow gap closure box and set the Max Track segment gap closing to 4uM and Max frame gap to 3, then click Next \n 9. Wait for the tracks to be calculated, then click Next \n 10. Add some filters to the tracks using the green plus button e.g. Track duration set to Above 3.45, then click Next \n 11. Do not change any Display options, click Next \n 12. Click Next \n 13. Under Select an action chose 'Export label image' and click Execute \n 14. Select 'Export only spots in tracks' \n Close the TrackMate window and click okay on this message");

	run("Duplicate...", "duplicate");
	rename("Segmentation");
	run("glasbey inverted");
	run("Edit LUT...");

	if(Number_of_channels==4) {
	run("Merge Channels...", "c1=[C1-" + image_title + ".tif" + "] c2=[C2-" + image_title + ".tif" + "] c3=[C3-" + image_title + ".tif" + "] c4=[C4-" + image_title + ".tif" + "] c5=Segmentation create");
	} else if (Number_of_channels==3) {
	run("Merge Channels...", "c1=[C1-" + image_title + ".tif" + "] c2=[C2-" + image_title + ".tif" + "] c3=[C3-" + image_title + ".tif" + "]  c4=[" + Segmentation + "] create");
	} else if (Number_of_channels==5) {
		run("Merge Channels...", "c1=[C1-" + image_title + ".tif" + "] c2=[C2-" + image_title + ".tif" + "] c3=[C3-" + image_title + ".tif" + "] c4=[C4-" + image_title + ".tif" + "] c5=[C5-" + image_title + ".tif" + "] c6=[" + Segmentation + "] create");
	} else {
		print("Error image must have between 3 and 5 channels");
	}

	saveAs("Tiff", segmentation_output_folder + "Segmented_" + image_title + ".tif");

	selectImage("Segmented_" + image_title + ".tif");
	segmented_image_title = getTitle();
	run("Split Channels");
	
	run("3D Manager Options", "volume surface compactness fit_ellipse integrated_density mean_grey_value std_dev_grey_value minimum_grey_value maximum_grey_value centroid_(pix) centroid_(unit) distance_to_surface centre_of_mass_(pix) centre_of_mass_(unit) bounding_box radial_distance closest distance_between_centers=10 distance_max_contact=1.80 drawing=Contour");
	

	
	if(Number_of_channels==4) {
			selectWindow("C5-" + segmented_image_title);
			run("3D Manager");
			Ext.Manager3D_AddImage();
			Ext.Manager3D_SelectAll();
			Ext.Manager3D_Measure();
			Ext.Manager3D_SaveResult("M", csvs_folder + "Morphological_measurements_" + image_title + ".csv");
			Ext.Manager3D_CloseResult("M");

			selectWindow("C1-" + segmented_image_title);
			Ext.Manager3D_Quantif();
			Ext.Manager3D_SaveResult("Q", csvs_folder + "C1_" + image_title + ".csv");
			Ext.Manager3D_CloseResult("Q");
			
			selectWindow("C2-" + segmented_image_title);
			Ext.Manager3D_Quantif();
			Ext.Manager3D_SaveResult("Q", csvs_folder + "C2_" + image_title + ".csv");
			Ext.Manager3D_CloseResult("Q");
		
			selectWindow("C3-" + segmented_image_title);
			Ext.Manager3D_Quantif();
			Ext.Manager3D_SaveResult("Q", csvs_folder + "C3_" + image_title + ".csv");
			Ext.Manager3D_CloseResult("Q");
		
			selectWindow("C4-" + segmented_image_title);
			Ext.Manager3D_Quantif();
			Ext.Manager3D_SaveResult("Q", csvs_folder + "C4_" + image_title + ".csv");
			Ext.Manager3D_CloseResult("Q");
		
			Ext.Manager3D_SelectAll();
			Ext.Manager3D_Delete();
			
			run("Clear Results");
			close();
			run("Close All");
	} else if (Number_of_channels==3) {
			selectWindow("C4-" + segmented_image_title);
			run("3D Manager");
			Ext.Manager3D_AddImage();
			Ext.Manager3D_SelectAll();
			Ext.Manager3D_Measure();
			Ext.Manager3D_SaveResult("M", csvs_folder + "Morphological_measurements_" + image_title + ".csv");
			Ext.Manager3D_CloseResult("M");
		
			selectWindow("C1-" + segmented_image_title);
			Ext.Manager3D_Quantif();
			Ext.Manager3D_SaveResult("Q", csvs_folder + "C1_" + image_title + ".csv");
			Ext.Manager3D_CloseResult("Q");
		
			selectWindow("C2-" + segmented_image_title);
			Ext.Manager3D_Quantif();
			Ext.Manager3D_SaveResult("Q", csvs_folder + "C2_" + image_title + ".csv");
			Ext.Manager3D_CloseResult("Q");
		
			selectWindow("C3-" + segmented_image_title);
			Ext.Manager3D_Quantif();
			Ext.Manager3D_SaveResult("Q", csvs_folder + "C3_" + image_title + ".csv");
			Ext.Manager3D_CloseResult("Q");
		
			Ext.Manager3D_SelectAll();
			Ext.Manager3D_Delete();
			
			run("Clear Results");
			close();
			run("Close All");
	} else if (Number_of_channels==5) {
		selectWindow("C6-" + segmented_image_title);
			run("3D Manager");
			Ext.Manager3D_AddImage();
			Ext.Manager3D_SelectAll();
			Ext.Manager3D_Measure();
			Ext.Manager3D_SaveResult("M", csvs_folder + "Morphological_measurements_" + image_title + ".csv");
			Ext.Manager3D_CloseResult("M");

			selectWindow("C1-" + segmented_image_title);
			Ext.Manager3D_Quantif();
			Ext.Manager3D_SaveResult("Q", csvs_folder + "C1_" + image_title + ".csv");
			Ext.Manager3D_CloseResult("Q");
			
			selectWindow("C2-" + segmented_image_title);
			Ext.Manager3D_Quantif();
			Ext.Manager3D_SaveResult("Q", csvs_folder + "C2_" + image_title + ".csv");
			Ext.Manager3D_CloseResult("Q");
		
			selectWindow("C3-" + segmented_image_title);
			Ext.Manager3D_Quantif();
			Ext.Manager3D_SaveResult("Q", csvs_folder + "C3_" + image_title + ".csv");
			Ext.Manager3D_CloseResult("Q");
		
			selectWindow("C4-" + segmented_image_title);
			Ext.Manager3D_Quantif();
			Ext.Manager3D_SaveResult("Q", csvs_folder + "C4_" + image_title + ".csv");
			Ext.Manager3D_CloseResult("Q");
			
			selectWindow("C5-" + segmented_image_title);
			Ext.Manager3D_Quantif();
			Ext.Manager3D_SaveResult("Q", csvs_folder + "C5_" + image_title + ".csv");
			Ext.Manager3D_CloseResult("Q");
		
			Ext.Manager3D_SelectAll();
			Ext.Manager3D_Delete();
			
			run("Clear Results");
			close();
			run("Close All");
	} else {
		print("Error image must have between 3 and 5 channels")
	}
	
}



	