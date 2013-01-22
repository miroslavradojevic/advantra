package advantra.swc;

import ij.ImagePlus;
import ij.io.FileSaver;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.text.SimpleDateFormat;
import java.util.Calendar;

import advantra.general.CreateDirectory;
import advantra.general.DebugExport;
import advantra.general.ImageConversions;
import advantra.shapes.Sphere;
import advantra.shapes.Sphere.Planisphere_Extr_Mode;
import advantra.trace.component.NeuronNode;
import advantra.trace.component.NeuronNode.Type;

public class AnalyzeSWC {
	
	// .swc file path
	private String  			file_path;
	private	String				folder_path;
	// source image path
	private ImagePlus			image;
	// some counters used
	private int					file_length;
	// structural elements
	private NeuronNode[]		nodes;	// each element is one node of the reconstruction
	
	private	int nr_bifurcation;
	private	int nr_end;
	private	int nr_body;
	private	int nr_undefined;
	
	public AnalyzeSWC(String swc_file_path, String corresponding_image_path){
		
		// make them absolute
		swc_file_path  = (new File(swc_file_path)).getAbsolutePath();
		corresponding_image_path = (new File(corresponding_image_path)).getAbsolutePath();
		
		// check whether files exist
		if(!(new File(swc_file_path).exists())){
			System.err.println(""+swc_file_path+" file does not exist!");
			return;
		}
		if(!(new File(corresponding_image_path).exists())){
			System.err.println(""+corresponding_image_path+" file does not exist!");
			return;
		}
		
		System.out.println("swc: "+swc_file_path);
		System.out.println("img: "+corresponding_image_path);
		
		this.file_path		= swc_file_path;
		this.folder_path	= (new File(file_path)).getParent(); // maybe add separator at the end
		this.image			= new ImagePlus(corresponding_image_path);
		
		// image has to be stack
		if(!(image.getStack().getSize()>1)){
			System.err.println(""+corresponding_image_path+" has to be stack!");
			System.exit(1);
		}
		
		this.file_length	= 0;
		
		this.nr_bifurcation	= 0;
		this.nr_end			= 0;
		this.nr_body		= 0;
		this.nr_undefined	= 0;
		
		// scan the file
		try {
			FileInputStream fstream 	= new FileInputStream(file_path);
			BufferedReader br 			= new BufferedReader(new InputStreamReader(new DataInputStream(fstream)));
			String read_line;
  			// check the length of the reconstruction first
	  		System.out.println("loading SWC...");
	  		while ( (read_line = br.readLine()) != null ) {
	  			if(!read_line.trim().startsWith("#")) file_length++; // # are comments
	  		}
	  		br.close();
		    fstream.close();
		}
		catch (Exception e){
  			System.err.println("Error: " + e.getMessage());
  		}
		System.out.println(file_length+" lines (nodes) found.");
		
		// allocate the "nodes"
		nodes 				= new NeuronNode[file_length];
		for (int i = 0; i < nodes.length; i++) {
			nodes[i]	= new NeuronNode();	// set nodes to zero at the initialization
		}
		
		System.out.format("Current folder is :  \t\t %s \n", folder_path);
		System.out.format("Current file is :    \t\t %s \n", file_path);
		System.out.println("Initialization done.");
		
	}
	
	public void load(){
		
		System.out.println("processing SWC...");
		
		// variables - data to be read 
  		int 		current_label 		= 0;
  		int 		current_type 		= 0;
  		double[] 	current_coords 		= new double[3];
  		double		current_radius 		= 0; 
  		int 		current_parent 		= 0;
		
  		String 		read_line;
  		int 		read_line_number	= 0;
  		
		// scan the file
		try {
			FileInputStream fstream 	= new FileInputStream(file_path);
			//fstream.getChannel().position(0);	// reset to the beginning of the file
			BufferedReader br 			= new BufferedReader(new InputStreamReader(new DataInputStream(fstream)));
			
			// should loop again the same amount of times it did at the initialization
			while ((read_line = br.readLine()) != null) {
				
				read_line = read_line.trim();
	  			if(!read_line.startsWith("#")){
		    		String[] tokens = read_line.split("\\s+");
		    		
		    		if(tokens.length == 7){
		    			
		    			current_label		= Integer.valueOf(tokens[0].trim()).intValue();
		    			current_type		= Integer.valueOf(tokens[1].trim()).intValue();
		    			
		    			current_coords[0] 	= Double.valueOf(tokens[2].trim()).doubleValue();
		    			current_coords[1] 	= Double.valueOf(tokens[3].trim()).doubleValue();
		    			current_coords[2] 	= Double.valueOf(tokens[4].trim()).doubleValue();
		    			current_radius		= Double.valueOf(tokens[5].trim()).doubleValue();
		    			current_parent 		= Integer.valueOf(tokens[6].trim()).intValue();
		    			
		    			read_line_number++; // counts the lines that were read
		    			
		    			if(read_line_number!=current_label){
		    				System.out.format("Warning: index didn't match with iteration number!\nrow: %s\n",read_line);
		    			}
		    			if(current_label>file_length){
		    				System.out.format("Warning: label was higher than the number of lines!\nrow: %s\n",read_line);
		    			}
		    			
		    			// use current_label to set 
		    			nodes[current_label-1].set(current_coords[0], current_coords[1], current_coords[2], current_radius);
		    			System.out.format(">>[%s]\n", read_line);
		    			System.out.format("\t>>Set node %d (idx %d) as %s \n", current_label, current_label-1, nodes[current_label-1].getType());
		    			if(current_label>1){
		    				nodes[current_parent-1].addNeighbour(current_coords[0], current_coords[1], current_coords[2]);
		    				System.out.format("\t\t>>Set node %d (idx %d) as %s \n", current_parent, current_parent-1, 
			    					nodes[current_parent-1].getType());
		    				
		    				nodes[current_label-1].addNeighbour(nodes[current_parent-1].getX(), nodes[current_parent-1].getY(), nodes[current_parent-1].getZ());
		    				System.out.format("\t\t\t>>Set node %d (idx %d) as %s \n", current_label, current_label-1, 
			    					nodes[current_label-1].getType());
		    			}
		    		}
		    		else{
		    			System.out.format("Exracted more than 7 elements... skipping row: \n%s\n", read_line);
		    		}
		    		
		    	}
	  		} // end looping the file
	  		
	  	br.close();
		fstream.close();	  		
  		
  		// all the values are read now
  		System.out.println("all nodes loaded.");
  		
		for (int i = 0; i < file_length; i++) {
			switch(nodes[i].getType()){
			case UNDEFINED_POINT:
				nr_undefined++;
				break;
			case END_POINT:
				nr_end++;
				break;
			case BODY_POINT:
				nr_body++;
				break;
			case BIFURCATION_POINT:
				nr_bifurcation++;
				break;
			default:
					System.err.println("AnalyzeSWC:load():This type of NeuronNode desn't exist");
					System.exit(1);
					break;
			}
		}
  		
  		}
		catch (Exception e){
		  			System.err.println("Error: " + e.getMessage());
		}
		
		for (int i = 0; i < current_coords.length; i++) {
			
		}
	}
	
	public int 	getNumberOfBifurcations(){
		return this.nr_bifurcation;
	}
	
	public int 	getNumberOfEndpoints(){
		return this.nr_end;
	}
	
	public int 	getNumberOfBodypoints(){
		return this.nr_body;
	}
	
	public int 	getNumberOfUndefined(){
		return this.nr_undefined;
	}
	
	public void extractBifurcations(int color){ // output will be in current folder, extract_bifurcations_swc_dd-MM-yyyy-HH_mm_ss
		
		String name_of_the_swc = (new File(file_path)).getName();
		String export_dir 	= 
				System.getProperty("user.dir")	      			+ File.separator +
				 name_of_the_swc+".bifurcations"				+ File.separator ;
		CreateDirectory.createOneDir(export_dir);
		
		// names of the 2 output files
		String bifur_log_path = export_dir+name_of_the_swc+".bif";
		String bifur_viz_path = export_dir+name_of_the_swc+".bif.tif";
		
		// create output rgb image and initialize it with source image stack
		ImagePlus img_rgb = ImageConversions.ImagePlusToRGB(image);
		Sphere sp = new Sphere();
		// create file with positions
		DebugExport bifur_log = new DebugExport(bifur_log_path);
		bifur_log.writeLine(String.format("# bifurcation_x bifurcation_y bifurcation_z neurite_r \n" +
				"# FORMAT TO STORE BIFURCATIONS"));
		for (int i = 0; i < file_length; i++) {
			if(nodes[i].getType()==Type.BIFURCATION_POINT){
				double row = nodes[i].getY();
				double col = nodes[i].getX();
				double lay = nodes[i].getZ();
				double rad = nodes[i].getR();
				bifur_log.writeLine(String.format("%f  %f  %f  %f ", row, col, lay, rad));
				sp.setSphere(row, col, lay, rad);
				sp.drawOverColorImage(img_rgb, color);
			}
		}
		bifur_log.closeDebug();
		new FileSaver(img_rgb).saveAsTiffStack(bifur_viz_path);
		System.out.println("exported:   \n"+bifur_log_path+"\n"+bifur_viz_path);
		
	}
	
	public void extractEndpoints(int color){ // output will be in current folder, extract_bifurcations_swc_dd-MM-yyyy-HH_mm_ss
		
		String name_of_the_swc = (new File(file_path)).getName();
		String export_dir 	= 
				System.getProperty("user.dir")	      			+ File.separator +
				 name_of_the_swc+".endpoints"				+ File.separator ;
		CreateDirectory.createOneDir(export_dir);
		
		// names of the 2 output files
		String endpt_log_path = export_dir+name_of_the_swc+".end";
		String endpt_viz_path = export_dir+name_of_the_swc+".end.tif";

		
		// create output rgb image and initialize it with source image stack
		ImagePlus img_rgb = ImageConversions.ImagePlusToRGB(image);
		Sphere sp = new Sphere();
		// create file with positions
		DebugExport endpt_log = new DebugExport(endpt_log_path);
		endpt_log.writeLine(String.format("# endpoint_x endpoint_y endpoint_z neurite_r \n" +
				"# FORMAT TO STORE ENDPOINTS"));
		for (int i = 0; i < file_length; i++) {
			if(nodes[i].getType()==Type.END_POINT){
				double row = nodes[i].getY();
				double col = nodes[i].getX();
				double lay = nodes[i].getZ();
				double rad = nodes[i].getR();
				endpt_log.writeLine(String.format("%f  %f  %f  %f ", row, col, lay, rad));
				sp.setSphere(row, col, lay, rad);
				sp.drawOverColorImage(img_rgb, 255);
			}
		}
		endpt_log.closeDebug();
		new FileSaver(img_rgb).saveAsTiffStack(endpt_viz_path);
		System.out.println("exported:   \n"+endpt_log_path+"\n"+endpt_viz_path);
	}
	
	public void extractPlanispheres(int resolution, Planisphere_Extr_Mode mode){ // will try to avoid direct planisphere extraction
		
		byte[][] output_image_values = new byte[image.getStack().getSize()][image.getHeight()*image.getWidth()];
		
		for (int i = 0; i < image.getStack().getSize(); i++) {
			// converted to comply with byte8 
			output_image_values[i] = (byte[])(ImageConversions.ImagePlusToGray8(image)).getStack().getPixels(i+1);
		}
		
		// create export folder
  		String folder_name = (new SimpleDateFormat("dd-MM-yyyy-HH_mm_ss")).format(Calendar.getInstance().getTime());
		String export_dir 	= 
				System.getProperty("user.dir")	+ File.separator +
				"extract_planispheres_swc_"+ folder_name		+ File.separator;
		CreateDirectory.createOneDir(export_dir);
		
		// create log file
		DebugExport planisph_log = new DebugExport(export_dir+"planispheres.log");
		
		for (int i = 0; i < file_length; i++) { // loops all nodes[] - i is node index
			
			String file_name = export_dir+String.format("%d_%s.tif", i, nodes[i].getType()); // name of the planisphere file
			
			(new FileSaver(extractPlanisphereAtNode(i, output_image_values, image.getWidth(), resolution, mode))).saveAsTiff(file_name); 
			System.out.format("node %d saved to %s  \n", i, file_name);
			
			double[] v1 = nodes[i].getV1();
			double[] v2 = nodes[i].getV2();
			double[] v3 = nodes[i].getV3();
			
			switch (nodes[i].getType()) {
				
				case UNDEFINED_POINT:
					planisph_log.writeLine(String.format("%d", i));
				break;
					
				case END_POINT:
					planisph_log.writeLine(String.format("%d %.3f %.3f %.3f", 
							i, v1[0], v1[1], v1[2]));
				break;	
				
				case BODY_POINT:
					planisph_log.writeLine(String.format("%d %.3f %.3f %.3f %.3f %.3f %.3f", 
							i, v1[0], v1[1], v1[2], v2[0], v2[1], v2[2]));
				break;
				
				case BIFURCATION_POINT:
					planisph_log.writeLine(String.format("%d %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f", 
							i, v1[0], v1[1], v1[2], v2[0], v2[1], v2[2], v3[0], v3[1], v3[2]));
				break;
				
				default:
					System.err.println("AnalizeSWC:extractPlanispheresOnNodes():this option for node type does not exist.");
					System.exit(1);
				break;
			}
			
			
		}
		
		planisph_log.closeDebug();
		
	}
	
	private ImagePlus 	extractPlanisphereAtNode(
			int node_index,//double node_x, double node_y, double node_z,
			byte[][] input_image, 
			int input_image_width,
			int resolution, 
			Planisphere_Extr_Mode mode){
		double planisphere_x = nodes[node_index].getY();
		double planisphere_y = nodes[node_index].getX();
		double planisphere_z = nodes[node_index].getZ();
		double planisphere_r = 3*nodes[node_index].getR();
		planisphere_r = (planisphere_r<2.0)?2.0:planisphere_r; // doesn't really make sense to take small spheres
		Sphere sph = new Sphere(planisphere_x, planisphere_y, planisphere_z, planisphere_r);
		//System.out.format("%f  %f, %f, %f \n", planisphere_x, planisphere_y, planisphere_z, planisphere_r);
		
		return sph.extractPlanisphereView(ImageConversions.toGray8(input_image, input_image_width), input_image_width, resolution, mode); 
	}
	
}
