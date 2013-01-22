package advantra.trace;

//import flanagan.math.ArrayMaths;
import flanagan.math.ArrayMaths;
import flanagan.math.Matrix;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.NewImage;

import advantra.general.DebugExport;
import advantra.processing.Functions;
import advantra.processing.IntensityCalc;
import advantra.shapes.Cylinder;
import advantra.shapes.Sphere;

public class TubeConfiguration{ // TODO: check whether it is used at all
	
	// is actually another representation of the cylinder 
	// contains image data, coordinates cross-section coordinates
	// and radial distances and cylinder instance used to allocate all of these
	
	int 			real_size;
	int[][] 		coord;
	int[] 			val;
	double[][] 		cross_section_x_y;
	double[]		radial_dist;	
	Cylinder		associated_cylinder;
	
	double 			estimatedRadius;
	boolean 		radiusEstConverged;
	
	public TubeConfiguration(Cylinder cylinder_shape){
		
		// main purpose is to allocate the memory for 
		// coord, val, cross_section_x_y, radial_dist
		// hence cylinder_shape should be biggest expected
		// this one does not fill the variables, just allocates
		
		associated_cylinder		= cylinder_shape;
		real_size 				= 0;
		
		int alloc_size          = Sphere.numberOfVoxInSphere(
				(int) Math.ceil(Math.sqrt(
						Math.pow(cylinder_shape.getR(), 2) + 
						Math.pow(cylinder_shape.getH(), 2)))
						);
		
		coord 					= new int		[3]	[alloc_size];
		val   					= new int			[alloc_size];
		cross_section_x_y 		= new double	[2]	[alloc_size];
		radial_dist				= new double		[alloc_size];
		
		estimatedRadius 	= 1.0;
		radiusEstConverged 	= false;
	}
	
	public int setTubeConfiguration(
			Cylinder cyl, 
			int[][] image_array,
			int image_height,
			int image_width,
			int image_len
			){
		
		// posible to happen that allocated size of the arrays is not big enough to support 
		// the cylinder real_size - that might be reported as a problem within the extractVox() method
		// hapens in case cyl is bigger than the cylinder_shape used in initialization
		
		associated_cylinder.setCylinder(cyl);
		
		// will set coord, val,	cross_section_x_y, radial_dist and real_size
		real_size = cyl.extractVox(
				image_array, 
				image_height, 
				image_width, 
				image_len, 
				coord, 
				val, 
				cross_section_x_y,
				radial_dist);
		
		return real_size;
		
	}
	
	public double tubeLikelihood(double Ic, double s, double Iprev, double Rprev, double std){
		
		double likelihood 		= 0;
		double max_diff_in_out 	= 0;
		double diff_in_out 		= 0;
		double in				= 0;
		
		if(
				associated_cylinder.getR()>=1 &&
				real_size>0
		){
			
			//for (int r = 1; r <= associated_cylinder.getR(); r++) {
				
				
				double[] avgInOut = IntensityCalc.averageIn_Out(
						cross_section_x_y,
						radial_dist,
						real_size,
						0, 0, Rprev,
						val
						);
				diff_in_out = Math.pow((avgInOut[0]-avgInOut[1])/Ic, s);
				//if(diff_in_out > max_diff_in_out){
					max_diff_in_out = diff_in_out;
				//	in = avgInOut[0];
				//}
				
			//}
			
			likelihood = max_diff_in_out * Functions.gaussianFn(in, Iprev, std);
			
		}
		
		return likelihood;
		
	}
	
	public double[] avgIntensityInOut(double radius){
		
		radius = (radius>associated_cylinder.getR())? associated_cylinder.getR() : radius ;
		radius = (radius<0)? 0 : radius ;
		
		double[] in_out = new double[2];
		double sum_in = 0, sum_out = 0;
		double cnt_in = 0, cnt_out = 0;
		
		for (int i = 0; i < real_size; i++) {
			if(radial_dist[i]<=radius){
				sum_in += val[i];
				cnt_in ++;
			}
			else{
				sum_out += val[i];
				cnt_out ++;
			}
		}
		
		in_out[0] = sum_in	/	cnt_in;
		in_out[1] = sum_out /	cnt_out;
		
		return in_out;
		
	}
	
	public void exportMatlab(String matlabFilename){
		
		System.out.println("writing configuration variables to "+matlabFilename+" file...");
		
		DebugExport dbg_exp = new DebugExport(matlabFilename);
		
		dbg_exp.write(String.format("coord = zeros( 3 , %d );", real_size));
		dbg_exp.write(String.format("val = zeros( 1 , %d );", real_size));
		dbg_exp.write(String.format("cs_x_y = zeros( 2 , %d );", real_size));
		dbg_exp.write(String.format("rad_dist = zeros( 1 , %d );", real_size));
		
		
		for (int i = 0; i < real_size; i++) {
			dbg_exp.writeLine(String.format("disp('load line %d / %d')", i, real_size));
			// export coord
			dbg_exp.write(String.format("coord( :, %d ) = [ %d ; %d ; %d ];", i+1, coord[0][i], coord[1][i], coord[2][i]));
			// export val
			dbg_exp.write(String.format("val( 1, %d ) = %d;", i+1, val[i]));
			// export cs_x_y
			dbg_exp.write(String.format("cs_x_y( :, %d ) = [ %f ; %f];", i+1, cross_section_x_y[0][i], cross_section_x_y[1][i]));
			// export radial_dist
			dbg_exp.write(String.format("rad_dist( 1, %d ) = %f;", i+1, radial_dist[i]));
			
			dbg_exp.writeLine("");
			
		}
		
		dbg_exp.closeDebug();
	}
	
	public double tubeConfigurationRadiusEst(double th, int buffer_size){
		                           
		// create radiuses
		double r_start 	= 1.0;
		double r_values = r_start;
		double r_step  	= 0.2;
		int    cnt 		= 1;
		
		// cnt will count how many of them exist
		while(r_values<associated_cylinder.getR()){
			cnt++;
			r_values+=r_step;
		}
		cnt--;
		
		double[] 	radiuses = new double[cnt]; // will contain radius values to be tested
		
		for (int i = 0; i < cnt; i++) {
			radiuses[i] = r_start + i * r_step;
		}

		// parameters threshold and filter width
//		double th 		= 5;
//		int buffer_size = 6;
		
		double[] 	energy_diff 	= new double[buffer_size];
		
		// int[]		s				= new int	[real_size]; 
		// this one will be defined for each radius
		// weighting each value - there are real_size values
		
		double energy = 0, energy_prev = 0;
		
		for (int i = 0; i < cnt; i++) {
			
			// radiuses[i] energy & energy_diff
			energy = 0;
			
			if(i==0){
				
				energy_prev = energy;
				
			}
			
			// calculate snakuscules energy
			for (int j = 0; j < real_size; j++) {
				energy += ((radial_dist[j]<radiuses[i])? -1 : +1 ) * val[j] ;
			}
			
			energy 		/= 8*radiuses[i]*radiuses[i];
			// shift buffer values
			for (int j = 0; j <= buffer_size-2; j++) {
				
				energy_diff[j] 	= energy_diff[j+1];
				
			}
			// fill the last buffer value
			energy_diff[buffer_size-1] 	= Math.abs(energy - energy_prev);
			
			
			if(i>=buffer_size){// buffer is full now
				
				// check whether all the energy_diff values from the buffer are small enough
				//boolean allSmallEnough = true;
				ArrayMaths am = new ArrayMaths(energy_diff);
				double max_diff = Math.abs(am.maximumDifference());
				
				//for (int k = 0; k < energy_diff.length; k++) {
				//	if(energy_diff[k]>=th){
				//		allSmallEnough = false;
				//		break;
				//	}
				//}
				
				if(max_diff<=th){
					estimatedRadius		= radiuses[i];
					radiusEstConverged	= true;
					System.out.println("converged!");
					return estimatedRadius;
					
				}
			}
			
			energy_prev = energy;
			
		}
		
		System.out.println("returning default...");
		estimatedRadius		= 1.0;
		radiusEstConverged 	= false;
		
		return estimatedRadius;

	}
	
	public void saveAsTiff(String path){
		
		int configuration_height = (int)(associated_cylinder.getH()+1);
		double[] max_dist = Matrix.rowMatrix(radial_dist).maximumElement();
		int configuration_size = (int) (Math.ceil(max_dist[0]*2)+1);
		
		ImagePlus img_conf = NewImage.createByteImage(
				"Configuration stack", 
				configuration_size, 
				configuration_size, 
				configuration_height, 
				NewImage.FILL_BLACK);
		
		//fill the configuration
		double[] unit_orientation = associated_cylinder.getV();
		
		double root_x = 
				associated_cylinder.getCenterX() - 
				(associated_cylinder.getH()/2) * 
				unit_orientation[0];
		double root_y = 
				associated_cylinder.getCenterY() - 
				(associated_cylinder.getH()/2) *
				unit_orientation[1];
		double root_z = 
				associated_cylinder.getCenterZ() - 
				(associated_cylinder.getH()/2) * 
				unit_orientation[2];
		
		for (int i = 0; i < real_size; i++) {
			
			int row = (int) Math.round(cross_section_x_y[0][i]+max_dist[0]); 
			row = (row >= configuration_size)? configuration_size-1 : row; 
			int col = (int) Math.round(cross_section_x_y[1][i]+max_dist[0]); 
			col = (col >= configuration_size)? configuration_size-1 : col; 
			// layer will be dot product of the vector origined in root_i
			int layer = (int)
					Math.round((coord[0][i]-root_x)*unit_orientation[0] +
					(coord[1][i]-root_y)*unit_orientation[1] +
					(coord[2][i]-root_z)*unit_orientation[2]);
			layer = (layer<0)? 0 : layer;
			layer = (layer>=configuration_height)? configuration_height-1 : layer;
			
			// now set the value taking one from the image
			img_conf.getStack().setVoxel(row, col, layer, (double)val[i]);
			
		}
		
		IJ.saveAs(img_conf, "Tiff", path);
		
	}

	public ImageStack toStack(){
		
		// this stack will be aligned along 
		
		
		int slices = (int)(associated_cylinder.getH());
		double[] max_radial_distance = Matrix.rowMatrix(radial_dist).maximumElement();
		int base_width = (int) (Math.ceil(max_radial_distance[0]*2)+1);
		
		ImagePlus 	conf_plus 	= NewImage.createByteImage("Blank Image", base_width, base_width, slices, NewImage.FILL_BLACK);
		ImageStack 	conf_stack 	= conf_plus.getStack();
		
		//fill the configuration
		double[] unit_orientation = associated_cylinder.getV();
		
		double root_x = 
				associated_cylinder.getCenterX() - 
				(associated_cylinder.getH()/2) * 
				unit_orientation[0];
		double root_y = 
				associated_cylinder.getCenterY() - 
				(associated_cylinder.getH()/2) *
				unit_orientation[1];
		double root_z = 
				associated_cylinder.getCenterZ() - 
				(associated_cylinder.getH()/2) * 
				unit_orientation[2];
		
		for (int i = 0; i < real_size; i++) {
			
			int row = (int) Math.round(cross_section_x_y[0][i]+max_radial_distance[0]); 
			row = (row >= base_width)? base_width-1 : row; 
			int col = (int) Math.round(cross_section_x_y[1][i]+max_radial_distance[0]); 
			col = (col >= base_width)? base_width-1 : col; 
			// layer will be dot product of the vector origined in root_i
			int layer = (int)
					Math.round((coord[0][i]-root_x)*unit_orientation[0] +
					(coord[1][i]-root_y)*unit_orientation[1] +
					(coord[2][i]-root_z)*unit_orientation[2]);
			layer = (layer<0)? 0 : layer;
			layer = (layer>=slices)? slices-1 : layer;
			
			// now set the value taking one from the image
			conf_stack.setVoxel(row, col, layer, (double)val[i]);

		}
		
		return conf_stack;
		
	}
	
 	public void saveCrossSectionAvgAsTiff(String path){

		double[] max_dist = Matrix.rowMatrix(radial_dist).maximumElement();
		int configuration_size = (int) (Math.ceil(max_dist[0]*2)+1);
		
		ImagePlus img_cs_avg = NewImage.createByteImage(
				"Configuration stack", 
				configuration_size, 
				configuration_size, 
				1, 
				NewImage.FILL_BLACK);
		
		byte[] img_cs_byte = (byte[])img_cs_avg.getStack().getPixels(1);
		
		// create average image 
		
		int[] img_cs_sum 		= new int[configuration_size*configuration_size];
		int[] img_cw_sum_count 	= new int[configuration_size*configuration_size];
		
		for (int i = 0; i < real_size; i++) {
			// extract the layer coordinates
			int row = (int) Math.round(cross_section_x_y[0][i] + max_dist[0]); 
			// first  row are x-coords of projections -> rows in image
			
			row = (row >= configuration_size)? configuration_size-1 : row; 
			
			int col = (int) Math.round(cross_section_x_y[1][i] + max_dist[0]); 
			// second row are y-coords of projections -> cols in image
			
			col = (col >= configuration_size)? configuration_size-1 : col; 
			// col can be [0, configuration_size-1]
			
			img_cs_sum[col+row*configuration_size] 			+= val[i];
			img_cw_sum_count[col+row*configuration_size] 	+= 1;
			
		}
		
		for (int i = 0; i < configuration_size*configuration_size; i++) {
			if(img_cw_sum_count[i]>0){
				
				img_cs_byte[i] = (byte)((img_cs_sum[i]/img_cw_sum_count[i]) & 0xff);
			
			}
		}
		
		// save it
		IJ.saveAs(img_cs_avg, "Tiff", path);
		
		
	}

	public int getRealLength() {
		return real_size;
	}
	
	public int[][] getCoord(){
		return coord;
	}
	
	public double[][] getCrossSectionXY(){
		return cross_section_x_y;
	}
	
	public int[] getVal(){
		return val;
	}
	
	public double[] getRadialDist(){
		return radial_dist;
	}
	
	public Cylinder getAssociatedCyl(){
		return associated_cylinder;
	}
}
