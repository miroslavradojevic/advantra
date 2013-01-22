package advantra.trials;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.Calendar;

import advantra.general.CreateDirectory;
import advantra.general.DebugExport;
import advantra.shapes.Sphere.Planisphere_Extr_Mode;
import advantra.tools.BranchModel;

import ij.IJ;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

public class CriticalPointDetection_ implements PlugIn {

	public static void main(String[] args){
		
		// default parameter valuess
		int width=128, height=128, length=32;
		double radius_std_dev=2.0;
		
		double planisphere_radius 			= 3*radius_std_dev;
		int planisphere_resolution 			= 32;
		int planisphere_extraction_mode 	= 1;
		
		if(args.length == 7){
			// params for model
			width					= (int)			Integer.parseInt(	args[0]);
			height					= (int)			Integer.parseInt(	args[1]);
			length					= (int)			Integer.parseInt(	args[2]);
			radius_std_dev			= (double)		Double.parseDouble(	args[3]);
			// params for planisphere transformation
			planisphere_radius		= (double)		Double.parseDouble(	args[4]);
			planisphere_resolution	= (int)			Integer.parseInt(	args[5]);
			planisphere_extraction_mode	= (int)		Integer.parseInt(	args[6]);
		}	
		else{
			System.err.println("Arguments: \n" +
					"1-width  2-height  3-length  4-radius_std_dev 5-planisphere_radius 6-planisphere_resolution 7-planisphere_extraction_mode \n" +
					"enter again...");
			System.exit(1);
		}
		
		/*
		 * check extraction mode parameter
		 */
		Planisphere_Extr_Mode mode = Planisphere_Extr_Mode.LOOP_SPHERICAL;
		switch(planisphere_extraction_mode){
		case 0: 
			mode = Planisphere_Extr_Mode.LOOP_CARTESIAN;
			break;
		case 1:
			mode = Planisphere_Extr_Mode.LOOP_SPHERICAL;
			break;
		default:
			System.err.println("planisphere_extraction_mode can be 0 or 1 only");
			System.exit(1);
			break;
		}
		
		/*
		 *  example
		 */
		
		String folder_name = (new SimpleDateFormat("dd-MM-yyyy-HH_mm_ss")).format(Calendar.getInstance().getTime());
		String export_dir 	= 
				System.getProperty("user.home")	+ File.separator +
				"gen_"+folder_name				+ File.separator;
		CreateDirectory.createOneDir(export_dir);
		
		/*
		 *  make the branch model
		 */
		BranchModel bm = new BranchModel(height, width, length, radius_std_dev);
		bm.drawHorizontalModel();
		bm.saveModel(export_dir);
		bm.exportModelParams(export_dir);
		
		/*
		 *  export planispheres 
		 */
		bm.exportPlanispheres(planisphere_radius, planisphere_resolution, mode, export_dir);
		//ImagePlus[] planisphere_images = bm.extractPlanispheres(planisphere_radius, planisphere_resolution, mode);
		
		/*
		// set of random points
		Random generator = new Random();
		double[][] set_of_points = new double[5][3];
		for (int i = 0; i < set_of_points.length; i++) {
			set_of_points[i][0] = generator.nextDouble() * height; 	// x
			set_of_points[i][1] = generator.nextDouble() * width; 	// y
			set_of_points[i][2] = generator.nextDouble() * length; 	// z
		}
		bm.exportPlanispheresAtPoints(set_of_points, planisphere_radius, planisphere_resolution, mode, export_dir);
		*/
		
		/*
		 *  mean-shift on planispheres 
		 */
		
		DebugExport mean_shift_log = new DebugExport(export_dir+"ms_export.m");  
		//int h_s = 3; 
		//int h_r = 20;
		//MeanShift ms = new MeanShift(planisphere_images[0], h_s, h_r, 5);
//		ms_log.writeLine(String.format("rows = []; cols = [];  \n"));
//		for (int i = 0; i < 10; i++) {
//			for (int j = 0; j < 20; j++) {
//				int[]res = ms.runAtPositionIntensityWeighted(i, j);
//				ms_log.writeLine(
//						String.format("rows = [rows; %d]; cols = [cols; %d];  \n", 
//								res[0], res[1]));//, 
//								//i+1, j+1, planispheres[pt].getProcessor().getPixelValue(i, j)));
//			}
//		}
//		ms_log.writeLine("im = imread('p2_planisphere.tif');figure; imshow(im); hold on;plot(cols,rows,'r.');");
		mean_shift_log.closeDebug();

		
		//BranchModel bm = new BranchModel;
		//bm.generate
		
	}
	
	public void run(String arg0) {
		
		// default param vals
		int width=128, height=128, length=32;
		double radius_std_dev=2.0;
		
		GenericDialog gd = new GenericDialog("Branch model formation", IJ.getInstance());
		// label, defaultVal, digits, columns, units
		gd.centerDialog(true);	
		gd.addMessage("Parameters ");
		gd.addNumericField( "width	:", 			width, 			0, 10, 		"" );
		gd.addNumericField( "height	:",				height, 		0, 10, 		""); 
		gd.addNumericField( "length	:",				length, 		0, 10, 		"");
		gd.addNumericField( "radius std. dev.:",	radius_std_dev, 1, 10, 		"" );
		
		
		gd.showDialog();
		if (gd.wasCanceled()) {
			return;
		}
		
		width 				= (int)gd.getNextNumber();
		height				= (int)gd.getNextNumber();
		length				= (int)gd.getNextNumber();
		radius_std_dev		= gd.getNextNumber();

		/*
		 *  example
		 */
		
		String folder_name = (new SimpleDateFormat("dd-MM-yyyy-HH_mm_ss")).format(Calendar.getInstance().getTime());
		String export_dir 	= 
				System.getProperty("user.home")	+ File.separator +
				"gen_"+folder_name				+ File.separator;
		CreateDirectory.createOneDir(export_dir);
		
		//BranchModel bm = new BranchModel;
		//bm.generate
		
	}
	
}