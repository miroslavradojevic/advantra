package advantra.commands;

import ij.ImagePlus;
import ij.io.FileSaver;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.Calendar;

import advantra.general.CreateDirectory;
import advantra.general.DebugExport;
import advantra.tools.MeanShift;

public class DemoMeanShiftAll {
	
public static void main(String[] args){
	
	String image_path 	= "";
	int h_spatial		= 4;
	int max_iter		= 100;
	double epsilon 		= 0.0001;
	
	if(args.length == 4){
		image_path		= 								args[0];
		h_spatial		= (int)		Integer.parseInt(	args[1]);
		max_iter		= (int)		Integer.parseInt(	args[2]);
		epsilon			= (double)	Double.parseDouble(	args[3]);
	}	
	else{
		System.err.println(
				"Takes the image and does mean-shift (noblurring) for predefined number of iterations. \n" +
				"Output (convergence points) are stored in a .csv file:\n" +
				"each row is one convergence point 1st element is row, 2nd is col, 3rd is how many points reached that conv. point \n\n"+
				"Arguments:\n" +
				"1-image path\n" +
				"2-h_spatial\n" +
				"3-max-iter\n" +
				"4-epsilon\n" +
				"Enter again...");
		System.exit(1);
	}
		
	// current time/date for the export file/folder
	String current_moment 	= (new SimpleDateFormat("dd-MM-yyyy-HH-mm-ss")).format(Calendar.getInstance().getTime()); // used for export filenames & export directory name
	String log_dir_name 	= System.getProperty("user.home")+File.separator+"DemoMeanShift."+current_moment+File.separator;
		
	// create dir & file
	CreateDirectory.createOneDir(log_dir_name);
	
	// check input path
	image_path = (new File(image_path)).getAbsolutePath();
	if( (!(new File(image_path)).exists()) ){
		System.err.println(image_path+"  does not exist!");
		System.exit(1);
	}
	// load image		
	ImagePlus img = new ImagePlus(image_path);
	
	if(img.getStack().getSize()>1){
		System.out.println("Image stack was at the input.");
		System.exit(1);
	}
	
	// export the image
	(new FileSaver(img)).saveAsTiff(log_dir_name+"input_image.tif");
	
	MeanShift 	ms 		= new MeanShift(img, h_spatial);
	
	for (int itr = 5; itr <= max_iter; itr+=1) {
			
		ms.run_WrapY(itr, epsilon); 
		
		int[][] clust = ms.extractClust(1);
		
		String output_filename = String.format("%siter%05d.meanshift", log_dir_name, itr);
		DebugExport metadata_file 	= new DebugExport(output_filename);
		for (int i = 0; i < clust.length; i++) {
			metadata_file.writeLine(String.format("%d, %d, %d", clust[i][0], clust[i][1], clust[i][2]));
		}
		metadata_file.closeDebug();
		System.out.println("\nlog exported to:    "+output_filename);
		
	}
	
	}

}
