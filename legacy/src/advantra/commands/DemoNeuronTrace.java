package advantra.commands;

import ij.ImagePlus;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.Calendar;

import advantra.general.CreateDirectory;
import advantra.trace.NeuronTrace;

public class DemoNeuronTrace {
	
	public static void main(String[] args){
		
		String image_stack_path 		= "";
		int start_pos_x					= 0;
		int start_pos_y					= 0;
		int start_pos_z					= 0;
		
		int branch_nr					= 20;
		
		if(args.length >= 4 && args.length <= 5){
			image_stack_path		= 								args[0];
			image_stack_path = (new File(image_stack_path)).getAbsolutePath();
			if( (!(new File(image_stack_path)).exists()) ){
				System.err.println(image_stack_path+"  does not exist!");
				System.exit(1);
			}
			
			start_pos_x				= (int)		Integer.parseInt(	args[1]);
			start_pos_y				= (int)		Integer.parseInt(	args[2]);
			start_pos_z				= (int)		Integer.parseInt(	args[3]);
			
			if(args.length==5){
				branch_nr 			= (int)		Integer.parseInt(	args[4]);
			}
		
		}	
		else{
			System.err.println(
					"Takes the image stack and does neuron tracing starting from the defined 3d point. \n" +
					"Arguments:\n" +
					"1-image stack path\n" +
					"2-x \n" +
					"3-y \n" +
					"4-z \n" +
					"5 (optional)-branch number \n" +
					"Enter again...");
			System.exit(1);
		}
		
		// load image stack		
		ImagePlus img 		= new ImagePlus(image_stack_path);
		
		if(!(img.getStack().getSize()>1) ){
			System.out.println("Image was not stack.");
			System.exit(1);
		}
		
		if(!(img.getType()==ImagePlus.GRAY8)){
			System.out.println("Image was not gray8.");
			System.exit(1);
		}
		
		// current time/date for the export file/folder
		String current_moment 	= (new SimpleDateFormat("dd-MM-yyyy-HH-mm-ss")).format(Calendar.getInstance().getTime()); 
		// used for export filenames & export directory name
		String log_dir_name 		= 
				System.getProperty("user.home")+File.separator+
				"DemoNeuronTrace."+current_moment+File.separator;
		
		// create dir & file
		CreateDirectory.createOneDir(log_dir_name);
		
		NeuronTrace neuron_tr = new NeuronTrace(img, branch_nr); 
		long start_time = System.currentTimeMillis();
		neuron_tr.trace(start_pos_x, start_pos_y, start_pos_z);
		System.out.format("elasped time: %f sec.\n", (double)(System.currentTimeMillis()-start_time)/1000 );
		
		neuron_tr.export_swc(log_dir_name+"recon.swc");
		
//		ImagePlus img_recon = neuron_tr.showReconstruction(RoiType.CYLINDER, 255, 0, 0); // RoiType.SPHERE
//		(new FileSaver(img)).saveAsTiffStack(log_dir_name+"input_image.tif");
//		(new FileSaver(img_recon)).saveAsTiffStack(log_dir_name+"reconstruction.tif");
	
	}
	

}
