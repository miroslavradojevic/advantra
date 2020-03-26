package advantra.commands;

import ij.ImagePlus;
import ij.io.FileSaver;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.Calendar;

import advantra.general.CreateDirectory;
import advantra.general.DebugExport;
import advantra.general.ImageConversions;
import advantra.shapes.Sphere;
import advantra.shapes.Sphere.Planisphere_Extr_Mode;

public class ExtractPlanisphereAtPos {

	public static void main(String[] args){
		
		double pos_x = 0;
		double pos_y = 0;
		double pos_z = 0;
		double pos_r = 0;
		String 	image_path 			= "";
		
		double r_ratio 				= 0.5;
		int 	resolution			= 32;
		int extraction_mode_index	= 1;
		
		switch(args.length){
		case 8:
			pos_x = (double) Double.parseDouble(args[0]);
			pos_y = (double) Double.parseDouble(args[1]);
			pos_z = (double) Double.parseDouble(args[2]);
			pos_r = (double) Double.parseDouble(args[3]);
			image_path 				= args[4];
			
			r_ratio					= Double.parseDouble(args[5]);
			resolution				= (int)	Integer.parseInt(args[6]);
			extraction_mode_index	= (int)	Integer.parseInt(args[7]);
			break;
		case 7:
			pos_x = (double) Double.parseDouble(args[0]);
			pos_y = (double) Double.parseDouble(args[1]);
			pos_z = (double) Double.parseDouble(args[2]);
			pos_r = (double) Double.parseDouble(args[3]);
			image_path 				= args[4];
			
			r_ratio					= Double.parseDouble(args[5]);
			resolution				= (int)	Integer.parseInt(args[6]);
			break;
		case 6:
			pos_x = (double) Double.parseDouble(args[0]);
			pos_y = (double) Double.parseDouble(args[1]);
			pos_z = (double) Double.parseDouble(args[2]);
			pos_r = (double) Double.parseDouble(args[3]);
			image_path 				= args[4];
			
			r_ratio					= Double.parseDouble(args[5]);
			break;
		case 5:
			pos_x = (double) Double.parseDouble(args[0]);
			pos_y = (double) Double.parseDouble(args[1]);
			pos_z = (double) Double.parseDouble(args[2]);
			pos_r = (double) Double.parseDouble(args[3]);
			image_path 				= args[4];
			break;
		default:
			System.err.println(
					"Takes image file, extracts the planisphere image at given position.\n"+
					"Arguments: \n" +
					"1 - x \n" +
					"2 - y \n" +
					"3 - z \n" +
					"4 - r \n" +
					"5 - image file path\n" +
					"6(OPTIONAL) - start_r    	0.1<start_r<0.9	\n" +
					"7(OPTIONAL) - resolution 			\n" +
					"8(OPTIONAL) - extraction mode(0-LOOP_CARTESIAN, 1-LOOP_SPHERICAL)\n" +
					"enter again...");
			System.exit(1);
			break;
			
		}
		
		// check extraction mode
		Planisphere_Extr_Mode extraction_mode  = Planisphere_Extr_Mode.LOOP_SPHERICAL; // default
		switch (extraction_mode_index) {
				case 0:
					extraction_mode = Planisphere_Extr_Mode.LOOP_CARTESIAN;
					break;
				case 1:
					extraction_mode = Planisphere_Extr_Mode.LOOP_SPHERICAL;
					break;	
				default:
					System.err.println("planisphere extraction mode index can be 0 or 1 only");
					System.exit(1);
					break;
		}
		
		image_path 	= (new File(image_path)).getAbsolutePath();
		if(!(new File(image_path)).exists()){
			System.err.println(""+image_path+" file does not exist!");
			System.exit(1);
		}
		
		// check resolution
		if(resolution<16){
					System.err.println("Resolution cannot be less than 16");
					System.exit(1);
		}
		
		// create export folder "extract_planisphere_dd-MM-yyyy-HH_mm_ss" in current folder
  		String folder_name = (new SimpleDateFormat("dd-MM-yyyy-HH_mm_ss")).format(Calendar.getInstance().getTime());
		String export_dir 	= 
				System.getProperty("user.dir")			+ File.separator +
				"extract_planisphere_"+ folder_name		+ File.separator;
		CreateDirectory.createOneDir(export_dir);
		
		// create log file
		DebugExport planisph_log = new DebugExport(export_dir+"planisphere.log");
		
		/*
		 * MAIN
		 */
		
		Sphere sp = new Sphere(pos_x, pos_y, pos_z, pos_r);
		
		System.out.format("extracting planisphere \nat (%f, %f, %f), radius %f \n", pos_x, pos_y, pos_z, pos_r);
		System.out.format("input img %s \n", image_path);
		System.out.format("staring from radius %f to %f \n", r_ratio*sp.getR(), sp.getR());
		System.out.format("resolution: %d \n", resolution);
		System.out.format("extraction mode: %s \n", extraction_mode);
		
		ImagePlus input_image = new ImagePlus(image_path);
		ImagePlus planisphere_output = sp.extractPlanisphereView(input_image, resolution, r_ratio, extraction_mode);		
		
		System.out.format("data exported to: %s\n", export_dir);
		
		(new FileSaver(planisphere_output)).saveAsTiff(export_dir+"planisphere.tif");
		
		planisph_log.writeLine(String.format("extracting planisphere \nat (%f, %f, %f), radius %f \n", pos_x, pos_y, pos_z, pos_r));
		planisph_log.writeLine(String.format("input img %s \n", image_path));
		planisph_log.writeLine(String.format("starting from radius %f to %f \n", r_ratio*sp.getR(), sp.getR()));
		planisph_log.writeLine(String.format("resolution: %d \n", resolution));
		planisph_log.writeLine(String.format("extraction mode: %s \n", extraction_mode));
		planisph_log.writeLine(String.format("extracted: %d voxels, distributed on %d x %d grid\n", sp.getExtractedPtsNr(), resolution, 2*resolution ));
		if(!(sp.getExtractedPtsNr()>0)){
			planisph_log.writeLine(String.format("sphere hits the borders or is out of the image... exported image is blank! \n"));
		}
		
		
		planisph_log.closeDebug();
		
		// save the color visualization
		ImagePlus img_rgb = ImageConversions.ImagePlusToRGB(input_image);
		sp.drawOverColorImage(img_rgb, 255);
		(new FileSaver(img_rgb)).saveAsTiffStack(export_dir+"where_is_planisphere.tif");
		
	}
	
}
