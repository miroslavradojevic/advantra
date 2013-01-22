package advantra.commands;

import ij.ImagePlus;
import ij.io.FileSaver;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;

import advantra.general.CreateDirectory;
import advantra.general.DebugExport;
import advantra.general.ImageConversions;
import advantra.shapes.Sphere;
import advantra.shapes.Sphere.Planisphere_Extr_Mode;

public class ExtractPlanisphereAtPoss {
	
	public static void main(String[] args){
		
		String	pos_file_path 		= ""; // file containing the coordinates & radiuses
		String 	image_path 			= ""; // source image
		double 	radius_scale 		= 2.5;
		int 	resolution 			= 32;
		double r_start				= 0.5;
		
		switch(args.length){
		case 5:
			pos_file_path 				= args[0];
			image_path					= args[1];
			radius_scale				= (double) Double.parseDouble(args[2]);
			resolution					= (int) Integer.parseInt(args[3]);
			r_start						= (double) Double.parseDouble(args[4]);
			break;
		default:
			System.err.println(
					"Takes the file with positions & extracts the planisphere images at given position.\n"+
					"Arguments: \n" +
					"1 - file with positions \n" +
					"2 - source image\n" +
					"3 - k...  planisphere_radius = k*neurite radius \n" +
					"4 - resolution \n" +
					"5 - r_start \n" +
					"enter again...");
			System.exit(1);
			break;
			
		}
		
		pos_file_path 	= (new File(pos_file_path)).getAbsolutePath();
		if(!(new File(pos_file_path)).exists()){
			System.err.println(""+pos_file_path+" file does not exist!");
			System.exit(1);
		}
		
		image_path 	= (new File(image_path)).getAbsolutePath();
		if(!(new File(image_path)).exists()){
			System.err.println(""+image_path+" file does not exist!");
			System.exit(1);
		}
		
		String export_dir 	= 
				System.getProperty("user.dir")									+ File.separator +
				(new File(pos_file_path)).getName() +	"_extract_planisphere"	+ File.separator ;
		CreateDirectory.createOneDir(export_dir);
		
		ImagePlus input_image = new ImagePlus(image_path);
		ImagePlus img_rgb = ImageConversions.ImagePlusToRGB(input_image);
		
		// read the input file (same reading principle as the one for swc)
		/*
		 *  SCAN THE LINES
		 */
//		int file_length	= 0;
//		try {
//					FileInputStream fstream 	= new FileInputStream(pos_file_path);
//					BufferedReader br 			= new BufferedReader(new InputStreamReader(new DataInputStream(fstream)));
//					String read_line;
//		  			// check the length of the reconstruction first
//			  		System.out.println("loading positions...");
//			  		while ( (read_line = br.readLine()) != null ) {
//			  			if(!read_line.trim().startsWith("#")) file_length++; // # are comments
//			  		}
//			  		br.close();
//				    fstream.close();
//		}
//		catch (Exception e){
//			System.err.println("Error: " + e.getMessage());
//		}
//		System.out.println(file_length+" lines (nodes) found.");
				
		// create log file
		DebugExport planisph_log = new DebugExport(export_dir+"planispheres.log");
		planisph_log.writeLine(String.format("input: %s ", image_path));
		
		/*
		 * READ THE LINES, POSITIONS
		 */

		try {
			FileInputStream fstream 	= new FileInputStream(pos_file_path);
			BufferedReader br 			= new BufferedReader(new InputStreamReader(new DataInputStream(fstream)));
			// should loop again the same amount of times it did at the initialization
			
			String read_line;
			int read_line_number = 0;
			double poss_x, poss_y, poss_z, poss_r;
			Sphere sp = new Sphere();
			
			while ((read_line = br.readLine()) != null) {
				
				read_line = read_line.trim();
			  	
				if(!read_line.startsWith("#")){
					String[] tokens = read_line.split("\\s+");
				    if(tokens.length == 4){
				    	poss_x 	= Double.valueOf(tokens[0].trim()).doubleValue();
				    	poss_y 	= Double.valueOf(tokens[1].trim()).doubleValue();
				    	poss_z 	= Double.valueOf(tokens[2].trim()).doubleValue();
				    	poss_r 	= Double.valueOf(tokens[3].trim()).doubleValue();
				    	
				    	sp.setSphere(poss_x, poss_y, poss_z, radius_scale*poss_r);
				    	// extract
				    	ImagePlus planisphere_output = sp.extractPlanisphereView(input_image, resolution, r_start, Planisphere_Extr_Mode.LOOP_SPHERICAL);
				    	// save it
				    	(new FileSaver(planisphere_output)).saveAsTiff(String.format("%s_planisphere_%d.tif", export_dir, read_line_number));
				    	// draw it
				    	sp.drawOverColorImage(img_rgb, 255);
				    	
						planisph_log.writeLine(String.format(
						"%d: planisphere at (%f, %f, %f), neurite radius %f, planisphere radius range (%f-%f (x%f)), res:%dx%d, extracted %d voxels",
						read_line_number, poss_x, poss_y, poss_z, poss_r, r_start*sp.getR(), sp.getR(), radius_scale, resolution, 2*resolution, sp.getExtractedPtsNr()
						));
				    	
				    	read_line_number++; // counts the lines that were read
				    			
				    }
				    else{
System.out.format(
		"Managed to read %d values for line: %s\n. Extracted more than 4 elements for position... skipping this row: \n%s\n", read_line);
				    }
				    		
				  }
			  } // end looping the file
		}
		catch (Exception e){
			System.err.println("Error: " + e.getMessage());
		}
		
		(new FileSaver(img_rgb)).saveAsTiffStack(export_dir+"where_are_planispheres.tif");
		System.out.println("data exported to: "+export_dir);

		planisph_log.closeDebug();
//		/*
//		 * MAIN lOOP
//		 */
//		
//		
//		for (int i = 0; i < file_length; i++) {
//			if(!(sp.getExtractedPtsNr()>0)){
//				planisph_log.writeLine(String.format("\tsphere hits the borders or is out of the image... exported image is blank! \n"));
//			}
//		}
//		
		
		
	}

}



//public class ExtractPlanispheresSWCPoints {
//
//	public static void main(String[] args){
//		
//		// params
//		String 	file_path 			= "";
//		String 	image_path 			= "";
//		int 	res 				= 0;
//		int extraction_mode_index	= 1;
//		
//
//		// make it absolute for later (folder will be referenced wrt. to this file - needs to be full address)
//		file_path 	= (new File(file_path)).getAbsolutePath();
//		image_path 	= (new File(image_path)).getAbsolutePath();
//		
//		// check paths
//		if(!(new File(file_path)).exists()){
//			System.err.println(""+file_path+" file does not exist!");
//			System.exit(1);
//		}
//		
//		if(!(new File(image_path)).exists()){
//			System.err.println(""+image_path+" file does not exist!");
//			System.exit(1);
//		}
//		
//		// check resolution
//		if(res<16){
//			System.err.println("Resolution cannot be less than 16");
//			System.exit(1);
//		}
//		
//		/*
//		 * analyze swc into nodes[]
//		 */
//		
//		AnalyzeSWC alyzer_swc = new AnalyzeSWC(file_path, image_path);
//		
//		alyzer_swc.load();
//		
//  		System.out.format("undefined: \t\t %d \n", 			alyzer_swc.getNumberOfUndefined());
//  		System.out.format("end-node: \t\t %d \n", 			alyzer_swc.getNumberOfEndpoints());
//  		System.out.format("body-node: \t\t%d \n", 			alyzer_swc.getNumberOfBodypoints());
//  		System.out.format("bifurcation-node: \t%d \n", 		alyzer_swc.getNumberOfBifurcations());
//  		int total = 
//  				alyzer_swc.getNumberOfUndefined()+
//  				alyzer_swc.getNumberOfEndpoints()+
//  				alyzer_swc.getNumberOfBodypoints()+
//  				alyzer_swc.getNumberOfBifurcations();
//  		System.out.format("total: \t\t\t%d \n", total);
//  		
//  		/*
//  		 * extract planispheres
//  		 */
//  		//ImagePlus img = new ImagePlus(image_path);
//  		alyzer_swc.extractPlanispheres(res, extraction_mode); // should make a new folder right where .swc is & export planipshere images there
//  		
//	}
//
//}
