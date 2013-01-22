package advantra.commands;

import java.io.File;

import advantra.shapes.Sphere.Planisphere_Extr_Mode;
import advantra.swc.AnalyzeSWC;

public class ExtractPlanispheresSWCPoints {

	public static void main(String[] args){
		
		// params
		String 	file_path 			= "";
		String 	image_path 			= "";
		int 	res 				= 0;
		int extraction_mode_index	= 1;
		
		// take params
		if(args.length == 4){
			file_path 				= args[0];
			image_path 				= args[1];
			res						= (int)		Integer.parseInt(	args[2]);
			extraction_mode_index	= (int)		Integer.parseInt(	args[3]);
		}
		else{
			System.err.println(
					"Takes .swc file, extracts the points, and planispheres at those locations.\n"+
					"Arguments: \n" +
					"1-.swc file path \n" +
					"2-.tif image path \n" +
					"3-resolution \n" +
					"4-planisphere extraction (0:loop cartesian, 1: loop spherical) \n" +
					"enter again...");
			System.exit(1);
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
		
		// make it absolute for later (folder will be referenced wrt. to this file - needs to be full address)
		file_path 	= (new File(file_path)).getAbsolutePath();
		image_path 	= (new File(image_path)).getAbsolutePath();
		
		// check paths
		if(!(new File(file_path)).exists()){
			System.err.println(""+file_path+" file does not exist!");
			System.exit(1);
		}
		
		if(!(new File(image_path)).exists()){
			System.err.println(""+image_path+" file does not exist!");
			System.exit(1);
		}
		
		// check resolution
		if(res<16){
			System.err.println("Resolution cannot be less than 16");
			System.exit(1);
		}
		
		/*
		 * analyze swc into nodes[]
		 */
		
		AnalyzeSWC alyzer_swc = new AnalyzeSWC(file_path, image_path);
		
		alyzer_swc.load();
		
  		System.out.format("undefined: \t\t %d \n", 			alyzer_swc.getNumberOfUndefined());
  		System.out.format("end-node: \t\t %d \n", 			alyzer_swc.getNumberOfEndpoints());
  		System.out.format("body-node: \t\t%d \n", 			alyzer_swc.getNumberOfBodypoints());
  		System.out.format("bifurcation-node: \t%d \n", 		alyzer_swc.getNumberOfBifurcations());
  		int total = 
  				alyzer_swc.getNumberOfUndefined()+
  				alyzer_swc.getNumberOfEndpoints()+
  				alyzer_swc.getNumberOfBodypoints()+
  				alyzer_swc.getNumberOfBifurcations();
  		System.out.format("total: \t\t\t%d \n", total);
  		
  		/*
  		 * extract planispheres
  		 */
  		//ImagePlus img = new ImagePlus(image_path);
  		alyzer_swc.extractPlanispheres(res, extraction_mode); // should make a new folder right where .swc is & export planipshere images there
  		
	}

}
