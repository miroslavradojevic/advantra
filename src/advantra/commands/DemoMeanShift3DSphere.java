package advantra.commands;

import ij.ImagePlus;
import ij.io.FileSaver;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.Calendar;

import advantra.general.ArrayHandling;
import advantra.general.CreateDirectory;
import advantra.general.ImageConversions;
import advantra.shapes.Sphere;
import advantra.tools.MeanShift3DSphere;

public class DemoMeanShift3DSphere {

	public static void main(String[] args){
		
		System.out.println("Uses extracted sphere ImagePlus object...");
		
		double pos_x 					= 0;
		double pos_y 					= 0;
		double pos_z 					= 0;
		double pos_r 					= 0;
		String 	image_path 				= "";
		
		int 	resolution				= 32;	// resolution x resolution x resolution
		
		if(args.length>=5 && args.length<=6){
			pos_x = (double) Double.parseDouble(	args[0]);
			pos_y = (double) Double.parseDouble(	args[1]);
			pos_z = (double) Double.parseDouble(	args[2]);
			pos_r = (double) Double.parseDouble(	args[3]);
			image_path 	= 							args[4];
			if(args.length>=6){
				resolution = (int)Integer.parseInt(	args[5]);
			}
		}
		else{
			System.err.println(
					"Takes image file, extracts the sphere at given position.\n"+
					"Arguments: \n" +
					"1 - x \n" +
					"2 - y \n" +
					"3 - z \n" +
					"4 - r \n" +
					"5 - image file path\n" +
					"6(OPTIONAL) - resolution (32 default)\n" +
					"enter again...");
			System.exit(1);
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
				System.getProperty("user.dir")					+ File.separator +
				"mean_shift_bifurcation_"+ folder_name			+ File.separator;
		CreateDirectory.createOneDir(export_dir);
		
		/*
		 * MAIN
		 */
		
		Sphere sp = new Sphere(pos_x, pos_y, pos_z, pos_r);
		
		System.out.format("extracting planisphere \nat (%f, %f, %f), radius %f \n", pos_x, pos_y, pos_z, pos_r);
		System.out.format("input img %s \n", image_path);
		System.out.format("resolution: %d \n", resolution);
		
		ImagePlus input_image 	= new ImagePlus(image_path);
		
		ImagePlus output 		= sp.extract(input_image, resolution);	
		(new FileSaver(output)).saveAsTiffStack(export_dir+"sphere_stack.tif");	

		double angle_range_deg 	= 20;
		double angle_range_rad 	= (angle_range_deg/180)*Math.PI;
		int 	N = 200;
		MeanShift3DSphere ms3dSph = new MeanShift3DSphere(output, sp, angle_range_rad, N);
		int max_iter = 100;
		double epsilon = 0.0001;
		ms3dSph.run(max_iter, epsilon);
		ms3dSph.extractClusters(0.05, 25);
		double[][] out_dirs = ms3dSph.getClusterDirs();
		if(out_dirs!=null){
			System.out.format("### DONE! %d directions! :\n", out_dirs.length);
			ArrayHandling.print2DArray(out_dirs);
		}
		else{
			System.out.format("no directions!\n");
		}
		
		// save the color visualization
		ImagePlus img_rgb = ImageConversions.ImagePlusToRGB(input_image);
		sp.drawOverColorImage(img_rgb, 255, 0, 0);
		(new FileSaver(img_rgb)).saveAsTiffStack(export_dir+"where_is_sphere.tif");

		System.out.println("files exported in :\n" +
				export_dir);
	}
	
}
