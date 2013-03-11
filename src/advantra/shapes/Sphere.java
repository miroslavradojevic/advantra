package advantra.shapes;

import ij.ImagePlus;
import ij.gui.NewImage;
import ij.process.ColorProcessor;
import advantra.general.ArrayHandling;
import advantra.general.ArrayHandling.IdxMode;
import advantra.general.ColourTransf;
import advantra.general.ImageConversions;
import advantra.general.Transf;
import advantra.processing.IntensityCalc;

public class Sphere extends RegionOfInterest {

	double   	r;
	
	int 				numberOfExtractedPoints;
	
	static 	int			R_MAX = 20;
	static  int			R_MIN = 1;
	
	public static enum Planisphere_Extr_Mode {LOOP_SPHERICAL, LOOP_CARTESIAN};
	
	public Sphere(){ 	// dummy construction
		super(0, 0, 0);
		this.roi_type = RoiType.SPHERE;
		this.r = 1;
		this.numberOfExtractedPoints = -1; // eachtime some extract... method is called - to obtain the number of points taken
		
	}
	
	public Sphere(double x, double y, double z, double r){
		super(x, y, z);
		this.roi_type = RoiType.SPHERE;
		this.r = (r<R_MIN)? R_MIN : ((r>R_MAX)?R_MAX:r); 
		this.numberOfExtractedPoints = -1; 
	}
	
	public Sphere(double[] xyz, double r) {
		super(xyz[0], xyz[1], xyz[2]);
		this.roi_type = RoiType.SPHERE;
		this.r = (r<R_MIN)? R_MIN : ((r>R_MAX)?R_MAX:r);
		this.numberOfExtractedPoints = -1;
	}
	
	public Sphere(Sphere s){
		super(s.getCenterX(), s.getCenterY(), s.getCenterZ());
		this.roi_type = RoiType.SPHERE;
		this.r = (s.getR()<R_MIN)? R_MIN : ((s.getR()>R_MAX)?R_MAX:s.getR()); 
		this.numberOfExtractedPoints = -1; 
	}
	
	public void setSphere(double x, double y, double z, double r){
		this.x = x;
		this.y = y;
		this.z = z;
		this.r = (r<R_MIN)? R_MIN : ((r>R_MAX)?R_MAX:r);
		this.numberOfExtractedPoints = -1; 	// reset it
	}
	
	public void setSphere(double[] xyz, double r){
		this.x = xyz[0];
		this.y = xyz[1];
		this.z = xyz[2];
		this.r = (r<R_MIN)? R_MIN : ((r>R_MAX)?R_MAX:r);
		this.numberOfExtractedPoints = -1; 	// reset it
	}

	public void setR(double r){
		this.r = (r<R_MIN)? R_MIN : ((r>R_MAX)?R_MAX:r);
		this.numberOfExtractedPoints = -1;	// reset it
	}
	
	public double getR(){
		return this.r;
	}
	
	public int 	getExtractedPtsNr(){
		return this.numberOfExtractedPoints;
	}
	
	/*	##################################
	 *  
	 * 	##################################
	 */	
	
	public int extractVox(
			final int[][] image_stack, 
			final int height, 
			final int width, 
			final int len, 
			final int[][] 	roi_coord,
			final int[] 	roi_vals
			){  
		
		// index = y + x*width 
		// y ~ columns 	~ width
		// x ~ rows 	~ height
		
		// extracts Sphere region of interest from particular image stack given as int[][] image_stack
		// images are gray8 but converted to int for processing
		// sphere defines center and radius
		// roi_coord is 3 x count integer array where rows 1-3 contain voxel coordinates
		// roi_vals  is 1 x count with voxel intensities
		
		// method will extract sphere from the image_stack
		// around the shpere center, with sphere's radius r
		
		int radius 			= (int)Math.round(this.r);
		int vox_in_sphere 	= Sphere.numberOfVoxInSphere((int) Math.ceil(r)); // expected number of voxels for given radius

		if(roi_coord[0].length<vox_in_sphere){
			System.err.println("Sphere:extractVox():\n" +
					"roi needs to have size 3 x"+Sphere.numberOfVoxInSphere((int) Math.ceil(r))+"+ columns. Allocate some more space!");
			System.exit(1);
		}
		
		if(roi_coord[0].length != roi_vals.length){			
			System.err.println("Sphere:extractVox():\n roi_coord needs to have the same size as roi_vals.");
			System.exit(1);
		}
		
		if(radius>R_MAX){
			System.err.println("Sphere:extractVox():\n radius cannot be more than "+R_MAX+" pixels...");
			System.exit(1);
		}
		
		// center in pixel coordinates
		int[] centerPix = {(int)Math.round(this.x), (int)Math.round(this.y), (int)Math.round(this.z)};
		
		// define the range
		int startX = 0; if (centerPix[0]-radius > 0) startX = centerPix[0]-radius;
		int startY = 0; if (centerPix[1]-radius > 0) startY = centerPix[1]-radius;
		int startZ = 0; if (centerPix[2]-radius > 0) startZ = centerPix[2]-radius;
			
		int endX   = height;if (centerPix[0]+radius < height) 	endX = centerPix[0]+radius;
		int endY   = width; if (centerPix[1]+radius < width) 	endY = centerPix[1]+radius;
		int endZ   = len;  	if (centerPix[2]+radius < len) 		endZ = centerPix[2]+radius;

		int count=0;
		
		for(int x=startX; x<endX; x++){
			for(int y=startY; y<endY; y++){
				for(int z=startZ; z<endZ; z++){
					// pixels = (byte[])stack.getPixels(z_layer);
					if((z-centerPix[2])*(z-centerPix[2])+(x-centerPix[0])*(x-centerPix[0])+(y-centerPix[1])*(y-centerPix[1])<=radius*radius){
						
						roi_coord[0][count] = x; 
						roi_coord[1][count] = y; 
						roi_coord[2][count] = z;
						
						roi_vals[count]  	= image_stack[z][y + width * x];
						
						count ++;
					}
				}
			}
		}	
		
		if(count<=0){
			System.err.println("Sphere:extractVox():\n Extracted number of voxels from the sphere was <=0 !!!");
			System.exit(1);
		}
		
		// border condition if sphere reached the image border
		if(count<=0.5*vox_in_sphere){
			System.err.println("Sphere:extractVox():\n Sphere used to extact voxels reached the end of the stack! Ending program...");
			System.exit(1);
		}
		
		this.numberOfExtractedPoints = count;
		
		return count;
	}

	/*	##################################
	 *  EXTRACTION METHODS
	 * 	##################################
	 */	
	
	public double[][] 	extractSphericalCoords( // TODO: replace this one 
			final ImagePlus 			input_image 
			){
		
		//outputs double[][] spherical_coords 	where 1.r, 2.phi, 3.theta 
		// referenced by the sphere's center
		
		
		double r_range = 0.9;
		// define the range
		int startX = 0; if (Math.floor(x-r) > 0) startX = (int)Math.floor(x-r);
		int startY = 0; if (Math.floor(y-r) > 0) startY = (int)Math.floor(y-r);
		int startZ = 0; if (Math.floor(z-r) > 0) startZ = (int)Math.floor(z-r);
				
		int height = input_image.getStack().getHeight();
		int width  = input_image.getStack().getWidth();
		int lenth  = input_image.getStack().getSize();
						
		int endX   = height; if (Math.ceil(x+r) < height) 	endX = (int)Math.ceil(x+r);
		int endY   = width;  if (Math.ceil(y+r) < width) 	endY = (int)Math.ceil(y+r);
		int endZ   = lenth;  if (Math.ceil(z+r) < lenth) 	endZ = (int)Math.ceil(z+r);
		
		// allocate output - see how many there are
		int count=0;
		for(int x=startX; x<endX; x++){
			for(int y=startY; y<endY; y++){
				for(int z=startZ; z<endZ; z++){
					
					double pos_x = x-this.x;
					double pos_y = y-this.y;
					double pos_z = z-this.z;
					
					if(pos_x*pos_x+pos_y*pos_y+pos_z*pos_z<=r*r &&
							pos_x*pos_x+pos_y*pos_y+pos_z*pos_z>=r_range*r*r_range*r){
						count ++;
					}
				}
			}
		}
		double[][] spherical_coords = new double[count][3];
		
		int idx = 0;
		for(int x=startX; x<endX; x++){
			for(int y=startY; y<endY; y++){
				for(int z=startZ; z<endZ; z++){
					
					double pos_x = x-this.x;
					double pos_y = y-this.y;
					double pos_z = z-this.z;
					
					if(pos_x*pos_x+pos_y*pos_y+pos_z*pos_z<=r*r && 
							pos_x*pos_x+pos_y*pos_y+pos_z*pos_z>=r_range*r*r_range*r){
						
						Transf.cart2sph(pos_x, pos_y, pos_z, spherical_coords[idx]);
						//System.out.format("%d, %d, %d, %f, %f, %f \n", 
						//		x, y, z, spherical_coords[idx][0], spherical_coords[idx][1], spherical_coords[idx][2]);
						idx++;
					}
				}
			}
		}	
		if(count<=0){
			System.err.println("Sphere:extractSphericalCoords():\n" +
					"Extracted number of voxels from the sphere was <=0 !!!");
		}
		
		this.numberOfExtractedPoints = count;
		
		return spherical_coords;
		
	}

	public ImagePlus extract( 
			final IntensityCalc 		calc, 
			final int 					resolution
			){  

		// method will extract sphere from the image_stack around the sphere center, with sphere's radius r
		// images are gray8, sphere defines center and radius

		int H = calc.getImgHeight();//input_image.getHeight();
		int W = calc.getImgWidth();
		int L = calc.getImgLength();//input_image.getStack().getSize();
		
		//IntensityCalc calc = new IntensityCalc(input_image.getStack());
		byte[][] pix_layers = new byte[resolution][resolution*resolution];
		
		// 'resolution' corresponds to sphere radius
		for (int r = 0; r < resolution; r++) {
			for (int c = 0; c < resolution; c++) {
				for (int l = 0; l < resolution; l++) {
					
					float sc = (float)resolution/2;
					float relative_row = (r-(sc-0.5f))/sc;
					float relative_col = (c-(sc-0.5f))/sc;
					float relative_lay = (l-(sc-0.5f))/sc;
					float relative_dist = (float)Math.sqrt(Math.pow(relative_row, 2)+Math.pow(relative_col, 2)+Math.pow(relative_lay, 2));	
					
					if(relative_dist<=1.00 ){ //  && relative_dist>=0.0
						// belongs to sphere shell at predefined radius
						// take the interpolated value from the image
						float coord_row = relative_row*(float)this.r + (float)this.x;
						float coord_col = relative_col*(float)this.r + (float)this.y; 
						float coord_lay = relative_lay*(float)this.r + (float)this.z;
//						System.out.format("\n%f, %f, %f (%d, %d, %d)...", 
//								coord_row, coord_col, coord_lay, r, c, l);
						if(coord_row>=0 && coord_row<=(H-1) && coord_col>=0 && coord_col<=(W-1) && coord_lay>=0 && coord_lay<=(L-1)){
							// take the interpolated value
							float value = calc.interpolateAt_new(coord_row, coord_col, coord_lay);
							pix_layers[l][r*resolution+c] = (byte)((int)Math.round(value));
//							System.out.format("... added %f (byte %d)", 
//									value, pix_layers[l][r*resolution+c]);
						}
					}
				}
			}
		}

		ImagePlus output 		= NewImage.createByteImage(
				"exported_sphere", resolution, resolution, resolution, NewImage.FILL_BLACK);
		for (int i = 0; i < resolution; i++) {
			output.getStack().setPixels(pix_layers[i], (i+1));
		}
		
		return output;
		
	}	
	
	public ImagePlus extractPlanisphereView( 	
			final ImagePlus 			source_image, 
			final int 					resolution,
			final double				r_range,
			final Planisphere_Extr_Mode how_to_extract
			){  
		
		double r_ratio = (r_range>=0.9)? 0.9 : r_range;
		r_ratio = (r_range<=0.1)? 0.1 : r_range;

		int radius 			= (int)Math.round(this.r);
		//int vox_in_sphere 	= Sphere.numberOfVoxInSphere((int) Math.ceil(r)); // expected number of voxels for given radius
		
		// center in pixel coordinates
		int[] centerPix = {(int)Math.round(this.x), (int)Math.round(this.y), (int)Math.round(this.z)};
		
		// allocate for the output
		double[] planisphere 			= 	new double[resolution*resolution*2]; // output will be Resolution x (2*Resolution)	 
		double[] planisphere_count 		= 	new double[resolution*resolution*2]; // sums
		int[] output_vals 				= 	new int[planisphere.length];
		
		// define the sizes
		int height 	= source_image.getHeight();
		int width	= source_image.getWidth();
		int length	= source_image.getStack().getSize();
		
		if (!(length>1)){
			System.err.println("Sphere:extractPlanisphereView():'source_image' has to be stack");
			System.exit(1);
		}
		
		boolean isInImage = 
				(centerPix[0]+radius < height)  &&
				(centerPix[1]+radius < width)   &&
				(centerPix[2]+radius < length)  &&
				(centerPix[0]-radius >= 0)      &&
				(centerPix[1]-radius >= 0)      &&
				(centerPix[2]-radius >= 0);
				
		// some auxiliary variables
		int 	index_theta=0, 	index_phi=0;
		// for conversions between spherical/cartesian
		double[] local_r_phi_theta 	= new double[3];
		double[] local_x_y_z 		= new double[3];
		
		int count=0;
		
		/*
		 * it is possible to loop&sample in cartesian space & assign values in spherical 
		 * or loop in spherical & sample values from corresponding cartesian coordinate  
		 */
		
		if(isInImage){
		switch(how_to_extract){
		case LOOP_CARTESIAN: 
			
			
			// TODO: border conditions have to be revised! what if the sphere is completely out of the range?!
			// here it is a bit special because sphere is independent from image
			int startX = centerPix[0]-radius;
			int startY = centerPix[1]-radius;
			int startZ = centerPix[2]-radius;
				
			int endX = centerPix[0]+radius;
			int endY = centerPix[1]+radius;
			int endZ = centerPix[2]+radius;
			
			// for each cartesian coordinate
			for(int x=startX; x<=endX; x++){
				for(int y=startY; y<=endY; y++){
					for(int z=startZ; z<=endZ; z++){
						
						local_x_y_z[0] = x-centerPix[0];
						local_x_y_z[1] = y-centerPix[1];
						local_x_y_z[2] = z-centerPix[2];
						
						Transf.cart2sph(
								local_x_y_z[0], 
								local_x_y_z[1], 
								local_x_y_z[2], 
								local_r_phi_theta);
						
						if((local_r_phi_theta[0]<=1.0*r)&&(local_r_phi_theta[0]>=r_ratio*r)){ 
							index_phi 		= ArrayHandling.value2index(local_r_phi_theta[1], 	IdxMode.LAST_INCLUDED, 0.0, Math.PI, 	resolution);
							index_theta 	= ArrayHandling.value2index(local_r_phi_theta[2],	IdxMode.LAST_EXCLUDED, 0.0, 2*Math.PI, 	2*resolution);
							// output width is 2*Resolution
							int idx = ArrayHandling.sub2index_2d(index_phi, index_theta, 2*resolution);
							planisphere[idx] 
									+= (int)source_image.getStack().getProcessor(z+1).getPixel(y, x);
							planisphere_count[ArrayHandling.sub2index_2d(index_phi, index_theta, 2*resolution)] ++;
							count ++;
						}
						
					}
				}
			}
			break;
			
			
		case LOOP_SPHERICAL:
			
			IntensityCalc im_calc = new IntensityCalc(source_image.getStack());
			
			for (double local_r = r_ratio*r; local_r <= 1.0*r; local_r+=0.1*r) {
				for (index_phi = 0; index_phi < resolution; index_phi++) {
					for (index_theta = 0; index_theta < 2*resolution; index_theta++) {
						
						local_r_phi_theta[0] = local_r;
						local_r_phi_theta[1] = ArrayHandling.index2value(index_phi, 	IdxMode.LAST_INCLUDED, 0.0, Math.PI, 	resolution);
						local_r_phi_theta[2] = ArrayHandling.index2value(index_theta,	IdxMode.LAST_EXCLUDED, 0.0, 2*Math.PI, 	2*resolution);
						
						Transf.sph2cart(
								local_r_phi_theta[0], 
								local_r_phi_theta[1], 
								local_r_phi_theta[2], 
								local_x_y_z);
						
						double x = this.x + local_x_y_z[0];
						double y = this.y + local_x_y_z[1];
						double z = this.z + local_x_y_z[2];
						
						if((x>=0 && x<=height-1) && (y>=0 && y<=width-1) && (z>=0 && z<=length-1)){
							
							int idx = ArrayHandling.sub2index_2d(index_phi, index_theta, 2*resolution);
							
							planisphere[idx] += 
									im_calc.interpolateAt_new((float)x, (float)y, (float)z);
									 
							planisphere_count[idx] += 1; 												 
							count ++;
						}
					}
				}
			}
			
			break;
		default:
			System.err.println("Sphere:extractPlanisphereView():this extraction mode is not possible.");
			System.exit(1);
			break;
		} // switch
		
		this.numberOfExtractedPoints = count;
		for (int i = 0; i < planisphere.length; i++) {
			if(planisphere_count[i]>0){
				output_vals[i] = (int)(planisphere[i] / planisphere_count[i] + 0.5);
			}
		}
		
		}// if isInImage otherwise give zeros
		
		ImagePlus planisph_imp = ImageConversions.toGray8(output_vals, 2*resolution);
		
		return planisph_imp;
	}	
	
	public int extractVox(
			final ImagePlus input_image, 
			final int[][] 	roi_coord,
			final int[] 	roi_vals
			){  
		
		// index = y + x*width 
		// y ~ columns 	~ width
		// x ~ rows 	~ height
		
		final int MAX_R = 20;
		//final int MIN_SPHERE_VOX_NR = 7;
		
		// extracts Sphere region of interest from particular input image
		// images are gray8 ImagePlus
		// sphere defines center and radius
		// roi_coord is 3 x count integer array where rows 1-3 contain voxel coordinates
		// roi_vals  is 1 x count with voxel intensities
		
		// method will extract sphere from the image_stack
		// around the shpere center, with sphere's radius r
		
		int radius = (int)Math.round(this.r);
		
		// checkings...
		if(roi_coord[0].length<Sphere.numberOfVoxInSphere((int) Math.ceil(r)) ){
			System.err.println("Sphere:extractVox():\n" +
					"roi_coord needs to have size 3 x"+Sphere.numberOfVoxInSphere((int) Math.ceil(r))+" at least! allocate some more space!");
			System.exit(1);
		}
		
		if(roi_coord.length!=3){
			System.err.println("Sphere:extractVox():\n roi_coord needs to have 3 rows.");
			System.exit(1);
		}
		
		if(roi_coord[0].length != roi_vals.length){			
			System.err.println("Sphere:extractVox():\n roi_coord needs to have the same size as roi_vals.");
			System.exit(1);
		}
		
		if(radius>MAX_R){
			System.err.println("Sphere:extractVox():\n radius cannot be more than "+MAX_R+" pixels...");
			System.exit(1);
		}
		
		// center in pixel coordinates
		int[] centerPix = {(int)Math.round(this.x), (int)Math.round(this.y), (int)Math.round(this.z)};
		
		// define the range
		int startX = 0; if (centerPix[0]-radius > 0) startX = centerPix[0]-radius;
		int startY = 0; if (centerPix[1]-radius > 0) startY = centerPix[1]-radius;
		int startZ = 0; if (centerPix[2]-radius > 0) startZ = centerPix[2]-radius;
		
		int height = input_image.getStack().getHeight();
		int width  = input_image.getStack().getWidth();
		int   len  = input_image.getStack().getSize();
				
		int endX   = height;if (centerPix[0]+radius < height) 	endX = centerPix[0]+radius;
		int endY   = width; if (centerPix[1]+radius < width) 	endY = centerPix[1]+radius;
		int endZ   = len;  	if (centerPix[2]+radius < len) 		endZ = centerPix[2]+radius;

		int count=0;
		
		for(int x=startX; x<endX; x++){
			for(int y=startY; y<endY; y++){
				for(int z=startZ; z<endZ; z++){
					if((z-centerPix[2])*(z-centerPix[2])+(x-centerPix[0])*(x-centerPix[0])+(y-centerPix[1])*(y-centerPix[1])<=radius*radius){
						
						roi_coord[0][count] = x; 
						roi_coord[1][count] = y; 
						roi_coord[2][count] = z;
						
						roi_vals[count] = (int)Math.round(input_image.getStack().getVoxel(y, x, z));
						
						count ++;
					}
				}
			}
		}	
		if(count<=0){
			System.err.println("Sphere:extractVox():\n Extracted number of voxels from the sphere was <=0 !!!");
		}
		
		this.numberOfExtractedPoints = count;
		
		return count;
	}

	public int extractCoords(
			final ImagePlus input_image, 
			final int[][] roi_coord
			){  
		
		final int MAX_R = 20;
		//final int MIN_SPHERE_VOX_NR = 7;
		
		// extracts coordinates of the Sphere region of interest from particular image stack given as input
		// input image is gray8 ImagePlus
		// sphere defines center and radius
		// output is 3 x count integer array where rows 1-3 contain voxel coordinates
		
		// method will extract sphere coords from the image_stack
		// around the shpere center, with sphere's radius r
		
		int radius = (int)Math.round(this.r);
		
		// checkings...
		if(roi_coord[0].length<Sphere.numberOfVoxInSphere(radius)){
			System.err.println("Sphere:extractCoords():\n output needs to have size 3x"+Sphere.numberOfVoxInSphere(radius)+" at least! allocate some more space!");
			System.exit(1);
		}
		
		if(roi_coord.length!=3){
			System.err.println("Sphere:extractCoords():\n output needs to have 3 rows.");
			System.exit(1);
		}
		
		if(radius>MAX_R){
			System.err.println("Sphere:extractCoords():\n radius cannot be more than "+MAX_R+" pixels...\n" +
					"setting it to "+MAX_R+" \n");
			radius = MAX_R;
		}
		
		// center in pixel coordinates
		int[] centerPix = {(int)Math.round(this.x), (int)Math.round(this.y), (int)Math.round(this.z)};
		
		// define the range
		int startX = 0; if (centerPix[0]-radius > 0) startX = centerPix[0]-radius;
		int startY = 0; if (centerPix[1]-radius > 0) startY = centerPix[1]-radius;
		int startZ = 0; if (centerPix[2]-radius > 0) startZ = centerPix[2]-radius;
		
		int height = input_image.getStack().getHeight();
		int width  = input_image.getStack().getWidth();
		int len    = input_image.getStack().getSize();
			
		int endX   = height;if (centerPix[0]+radius < height) 	endX = centerPix[0]+radius;
		int endY   = width; if (centerPix[1]+radius < width) 	endY = centerPix[1]+radius;
		int endZ   = len;  	if (centerPix[2]+radius < len) 		endZ = centerPix[2]+radius;

		int count=0;
		
		for(int x=startX; x<endX; x++){
			for(int y=startY; y<endY; y++){
				for(int z=startZ; z<endZ; z++){
					if((z-centerPix[2])*(z-centerPix[2])+(x-centerPix[0])*(x-centerPix[0])+(y-centerPix[1])*(y-centerPix[1])<=radius*radius){
						
						roi_coord[0][count] = x; 
						roi_coord[1][count] = y; 
						roi_coord[2][count] = z;
						
						count ++;
					}
				}
			}
		}	
		if(count<=0){
			System.out.println("Sphere:extractSphereCoords():\n Extracted number of voxels from the sphere was <=0 !!!");
			//System.exit(1);
		}
		
		this.numberOfExtractedPoints = count;
		
		return count;
	}
	
	/*	##################################
	 *  GENERATE pts
	 * 	##################################
	 */	
	
	public void 		generate3DSemiSpherePts(
			int N, 
			double ax, 
			double ay,
			double az,
			final double[][] 	points,
			final double[][] 	point_orients
			) {
		
		// checking...
		if(points.length!=3 || points[0].length!=N){
			System.err.println("Sphere:generate3DSemiSpherePts(): argument for storing 3d points has to have 3 rows and "+N+" columns.");
			System.exit(-1);
		}
		
		if(point_orients.length!=3 || point_orients[0].length!=N){
			System.err.println("Sphere:generate3DSemiSpherePts(): argument for storing 3d points has to have 3 rows and "+N+" columns.");
			System.exit(-1);
		}
		
		double h_k, theta_k, phi_k, phi_k_1 = 0;
		
		for (int k = 0; k < N; k++) {
			
			h_k = (double)k/(N-1); // 0 : 1
			
			theta_k = Math.acos(h_k);
			
			if(k==0 || k==(N-1)){
				
				phi_k   = 0;
				phi_k_1 = 0;
			
			}
			else{
				
				phi_k = phi_k_1 + 3.6 / ( Math.sqrt(N) * Math.sqrt(1-h_k*h_k));
				phi_k_1 = phi_k;
				
			}
			
			// cartesian coordinates
			point_orients[0][k] = points[0][k] = Math.sin(theta_k) * Math.cos(phi_k);
			point_orients[1][k] = points[1][k] = Math.sin(theta_k) * Math.sin(phi_k);
			point_orients[2][k] = points[2][k] = Math.cos(theta_k);
			
			// scale it 
			points[0][k] = this.r * points[0][k];
			points[1][k] = this.r * points[1][k];
			points[2][k] = this.r * points[2][k];
			
		}
		// set the orientation according to ax, ay, az - for all the points
		double norm = Math.sqrt(Math.pow(ax, 2)+Math.pow(ay, 2)+Math.pow(az, 2));
		
		Transf.rotate3xN(ax/norm, ay/norm, az/norm, points);  // rotation result is stored in points	
		
		Transf.rotate3xN(ax/norm, ay/norm, az/norm, point_orients); // have to add this one so that they're properly oriented
		
		for (int k = 0; k < N; k++) {
			// set sphere center as origin - translation
			points[0][k] += this.x;
			points[1][k] += this.y;
			points[2][k] += this.z;
		}
		
	}

	public double[][] 	generate3DSemiSpherePts(
			int N, 
			double ax, 
			double ay,
			double az
			) {
		
		double h_k, theta_k, phi_k, phi_k_1 = 0;
		
		double[][] points = new double[3][N];
		
		for (int k = 0; k < N; k++) {
			
			h_k = (double)k/(N-1); // 0 : 1
			
			theta_k = Math.acos(h_k);
			
			if(k==0 || k==(N-1)){
				
				phi_k   = 0;
				phi_k_1 = 0;
			
			}
			else{
				
				phi_k = phi_k_1 + 3.6 / ( Math.sqrt(N) * Math.sqrt(1-h_k*h_k));
				phi_k_1 = phi_k;
				
			}
			
			// points in cartesian coordinates
			points[0][k] = Math.sin(theta_k) * Math.cos(phi_k);
			points[1][k] = Math.sin(theta_k) * Math.sin(phi_k);
			points[2][k] = Math.cos(theta_k);
			
			// scale it 
			points[0][k] = this.r * points[0][k];
			points[1][k] = this.r * points[1][k];
			points[2][k] = this.r * points[2][k];
			
		}
		// set the orientation according to ax, ay, az - for all the points
		double norm = Math.sqrt(Math.pow(ax, 2)+Math.pow(ay, 2)+Math.pow(az, 2));
		
		Transf.rotate3xN(ax/norm, ay/norm, az/norm, points);  // rotation result is stored in points	
		
		for (int k = 0; k < N; k++) {
			// set sphere center as origin - translation
			points[0][k] += x;
			points[1][k] += y;
			points[2][k] += z;
		}
		
		return points;
		
	}
	
	public void 		generate3DFullSpherePts(
			int N, 
			final double[][] 	points,
			final double[][] 	point_orients
	){
		
		// checking...
		if(points.length!=3 || points[0].length!=N){
			System.err.println("Sphere:generate3DSemiSpherePts(): argument for storing 3d points has to have 3 rows and "+N+" columns.");
			System.exit(-1);
		}
		
		if(point_orients.length!=3 || point_orients[0].length!=N){
			System.err.println("Sphere:generate3DSemiSpherePts(): argument for storing 3d points has to have 3 rows and "+N+" columns.");
			System.exit(-1);
		}
		
		double h_k, theta_k, phi_k, phi_k_1 = 0;
		
		for (int k = 0; k < N; k++) {
			
			h_k = -1 + 2 * (double)k/(N-1); // -1 : 1
			
			theta_k = Math.acos(h_k);
			
			if(k==0 || k==(N-1)){
				
				phi_k   = 0;
				phi_k_1 = 0;
			
			}
			else{
				
				phi_k = phi_k_1 + 3.6 / ( Math.sqrt(N) * Math.sqrt(1-h_k*h_k));
				phi_k_1 = phi_k;
				
			}
			
			// cartesian coordinates
			points[0][k] = Math.sin(theta_k) * Math.cos(phi_k);
			points[1][k] = Math.sin(theta_k) * Math.sin(phi_k);
			points[2][k] = Math.cos(theta_k);
			
			// orientation vecs
			point_orients[0][k] = points[0][k];
			point_orients[1][k] = points[1][k];
			point_orients[2][k] = points[2][k];
			
			// scale it 
			points[0][k] = this.r * points[0][k];
			points[1][k] = this.r * points[1][k];
			points[2][k] = this.r * points[2][k];
			
			// set sphere center as origin - translation
			points[0][k] += this.x;
			points[1][k] += this.y;
			points[2][k] += this.z;
			
		}
		
	}

	public double[][] 	generate3DFullSpherePts(int N){ // 3xN
		
		double h_k, theta_k, phi_k, phi_k_1 = 0;
		
		double[][] points = new double[3][N];
		
		for (int k = 0; k < N; k++) {
			
			h_k = -1 + 2 * (double)k/(N-1); // -1 : 1
			
			theta_k = Math.acos(h_k);
			
			if(k==0 || k==(N-1)){
				
				phi_k   = 0;
				phi_k_1 = 0;
			
			}
			else{
				
				phi_k = phi_k_1 + 3.6 / ( Math.sqrt(N) * Math.sqrt(1-h_k*h_k));
				phi_k_1 = phi_k;
				
			}
			
			// cartesian coordinates
			points[0][k] = Math.sin(theta_k) * Math.cos(phi_k);
			points[1][k] = Math.sin(theta_k) * Math.sin(phi_k);
			points[2][k] = Math.cos(theta_k);
			
			// scale it & translate
			points[0][k] = r * points[0][k] + x;
			points[1][k] = r * points[1][k] + y;
			points[2][k] = r * points[2][k] + z;
			
		}
		
		return points;
		
	}
	
	public double[][] 	generate3DFullSpherePts_Nx3(int N){ 
		
		double h_k, theta_k, phi_k, phi_k_1 = 0;
		
		double[][] points = new double[N][3];
		
		for (int k = 0; k < N; k++) {
			
			h_k = -1 + 2 * (double)k/(N-1); // -1 : 1
			
			theta_k = Math.acos(h_k);
			
			if(k==0 || k==(N-1)){
				
				phi_k   = 0;
				phi_k_1 = 0;
			
			}
			else{
				
				phi_k = phi_k_1 + 3.6 / ( Math.sqrt(N) * Math.sqrt(1-h_k*h_k));
				phi_k_1 = phi_k;
				
			}
			
			// cartesian coordinates
			points[k][0] = Math.sin(theta_k) * Math.cos(phi_k);
			points[k][1] = Math.sin(theta_k) * Math.sin(phi_k);
			points[k][2] = Math.cos(theta_k);
			
			// scale it & translate
			points[k][0] = r * points[k][0] + x;
			points[k][1] = r * points[k][1] + y;
			points[k][2] = r * points[k][2] + z;
			
		}
		
		return points;
		
	}
	
	private int 		numberOfSolidAnglePoints(double angleRange, int K){
		
		int cnt = 1; // for the central one
		for (int k = 1; k <= K; k++) {
			double phi	= k*(angleRange/K);
			cnt += (int)Math.ceil(((2*Math.PI)*(Math.sin(phi)/(angleRange/K))  ));
		}
		return cnt;
		
	}
	
	public double[][] 	coordinatesSolidAngleSpherePts(double angleRange, int K, double phi, double theta) {
		// gives global coords
		double[][]  x_y_z = new double[numberOfSolidAnglePoints(angleRange, K)][3];
		
		Transf.sph2cart(r, 0.0, 0.0, x_y_z[0]);
		
		int cnt = 1;
		for (int k = 1; k <= K; k++) {
			double incl = k*(angleRange/K);
			int	N 		= (int)Math.ceil((  (2*Math.PI)*(Math.sin(incl)/(angleRange/K))  ));
			
			for (int n = 1; n <= N; n++) {
				double azim = n*((2*Math.PI)/N);
				
				Transf.sph2cart(r, incl, azim, x_y_z[cnt]);
				cnt++;
				
			}
		}
		
		Transf.rotateNx3(phi, theta, x_y_z);
		
		for (int i = 0; i < x_y_z.length; i++) {
			x_y_z[i][0] += this.x;
			x_y_z[i][1] += this.y;
			x_y_z[i][2] += this.z;
		}
		
		return x_y_z;
	}

	public double[][] 	cartesian2spherical(double[][] cart){// Nx3
		double[][] sph = new double[cart.length][3];
		for (int i = 0; i < cart.length; i++) {
			Transf.cart2sph(cart[i][0]-this.x, cart[i][1]-this.y, cart[i][2]-this.z, sph[i]); // does wrapping
		}
		return sph;
	}
	
	public double[] 	cartesian2spherical(double[] cart){
		double[] sph = new double[3];
		//for (int i = 0; i < cart.length; i++) {
			Transf.cart2sph(cart[0]-this.x, cart[1]-this.y, cart[2]-this.z, sph); // does wrapping
		//}
		return sph;
	}

	public double[][] 	spherical2cartesian(double[][] sph){
		double[][] cart = new double[sph.length][3];
		for (int i = 0; i < sph.length; i++) {
			Transf.sph2cart(sph[i][0], sph[i][1], sph[i][2], cart[i]);
			cart[i][0] += this.x; 
			cart[i][1] += this.y; 
			cart[i][2] += this.z; 
		}
		return cart;
	}
	
	public double[] 	spherical2cartesian(double[] sph){ // r, phi, theta
		double[] cart = new double[3];
		Transf.sph2cart(sph[0], sph[1], sph[2], cart);
		cart[0] += this.x; 
		cart[1] += this.y; 
		cart[2] += this.z;
		return cart;
	}
	// TODO: this method shouldn't be here but method in MeanShift3D
	public float 		avgValuePerDirection(double phi, double theta, IntensityCalc extracted_sphere_calc){
		
		float value 	= 0;
		//int count 		= 0;
		
		//for (double r = 0.95*getR(); r < getR(); r+=0.02*getR()) {
		double r = 0.95*getR();
			double x = this.x+Transf.sph2cart_x(r, phi, theta);
			double y = this.y+Transf.sph2cart_y(r, phi, theta);
			double z = this.z+Transf.sph2cart_z(r, phi, theta);
			value += extracted_sphere_calc.interpolateAt_new((float)x, (float)y, (float)z);
//			count++;
		//}
		
		return value;///(float)count;
	}
	
	public float 		sphereSurfaceValuePerDirection(double phi, double theta, IntensityCalc calc){ 
		
		float value 	= 0;
		int count 		= 0;
		
		for (double r = 0.9*getR(); r <= getR(); r+=0.05*getR()) {
			
			float x = (float)(this.x+Transf.sph2cart_x(r, phi, theta));
			float y = (float)(this.y+Transf.sph2cart_y(r, phi, theta));
			float z = (float)(this.z+Transf.sph2cart_z(r, phi, theta));
			value += calc.interpolateAt_new(x, y, z);
			count++;
			
		}
		
		return value/count;
	}

	public ImagePlus	avgValuesPerDirection(ImagePlus input_img, int nr_points, int resolution){
		
		ImagePlus output 		= NewImage.createByteImage(
				"average_values_per_direction", 2*resolution, resolution, 1, NewImage.FILL_BLACK);
		IntensityCalc img_calc = new IntensityCalc(input_img.getStack());
		
		// generate cartesian points on the sphere
		double[][] pts = generate3DFullSpherePts_Nx3(nr_points);
		// make them spherical
		double[][] pts_sph = cartesian2spherical(pts);
		for (int i = 0; i < pts_sph.length; i++) {
			// get indexes to be plotted on the image
			double current_phi = pts_sph[i][1];
			double current_theta = pts_sph[i][2];
			int row = ArrayHandling.value2index(current_phi, 	ArrayHandling.IdxMode.LAST_INCLUDED, 0.0, Math.PI, 	resolution);
			int col = ArrayHandling.value2index(current_theta, 	ArrayHandling.IdxMode.LAST_EXCLUDED, 0.0, 2*Math.PI, 2*resolution);
			 
			double g = avgValuePerDirection(current_phi, current_theta, img_calc);
			output.getStack().setVoxel(col, row, 0, g);
		}
		
		return output;
		
	}

	/*	##################################
	 *  GRAPHICS
	 * 	##################################
	 */	
	// if it is gray8  - valueToDraw 0-255
	// if it is gray16 - valueToDraw 0-65,535
	// use ImageConversions.setJet256ColorValue to draw over color images
	public ImagePlus 	drawOverGrayImage(ImagePlus template_image, int valueToDraw){
		
		int wt = template_image.getStack().getWidth();
		
		// check if it's grey image
		if(template_image.getType()==ImagePlus.GRAY8){
			// expects 0-255 8-bit scale
			valueToDraw = (valueToDraw>255) ? 255 : valueToDraw ;
			valueToDraw = (valueToDraw<  0) ?   0 : valueToDraw ;
		}
		else if(template_image.getType()==ImagePlus.GRAY16){
			// expects 0-65,535 16-bit scale
			valueToDraw = (valueToDraw>65535) ? 65535 : valueToDraw ;
			valueToDraw = (valueToDraw<  0  ) ?   0   : valueToDraw ;
		}
		else{
			System.err.println("Sphere:drawOverImage(): this image format is not suported - possible to use gray8, gray16 or rgb!");
			System.exit(-1);
		}
		
		int[][] template_image_array = ImageConversions.GraytoIntArray(template_image);
		
		// extract the points
		int[][] sphere_coordinates 	= 
				new int[3][Sphere.numberOfVoxInSphere((int) r)];
		
		int cnt_coords = extractCoords(template_image, sphere_coordinates);
		
		if(cnt_coords>0){
			
			// change the spots, assign them with valueToDraw
			for (int i = 0; i < cnt_coords; i++) {
				int row = sphere_coordinates[0][i];
				int col = sphere_coordinates[1][i];
				int layer = sphere_coordinates[2][i];
				template_image_array[layer][col+row*wt] = valueToDraw;
				
			}
			
		}else{
			
			System.out.println("Sphere:drawOverImage()\nThere was not enough points to draw anything!");
		
		}

		// to ImagePlus
		if(template_image.getType()==ImagePlus.GRAY8){
			return ImageConversions.toGray8(template_image_array, wt);
		}
		else if(template_image.getType()==ImagePlus.GRAY16){
			return ImageConversions.toGray16(template_image_array, wt);
			
		}
		else{
			
			System.err.println("Sphere:drawOverImage(): this image format is not suported - possible to use gray8, gray16 or rgb!");
			System.exit(-1);
			return ImageConversions.toGray8(template_image_array, wt); // just a dummy line to avoid error 
			
		}
		
	}
	
	public void 		drawOverColorImage(ImagePlus template_image, int valueToDraw){
		
		int w = template_image.getWidth();
		
		if((template_image.getType()!=ImagePlus.COLOR_RGB)){
			// set it to rgb
			template_image = ImageConversions.ImagePlusToRGB(template_image);
			System.out.println("converting ImagePlus to rgb...");
		}

		byte[] rgb_values = ColourTransf.Jet256(valueToDraw);
		
		byte[][][] img_array = ImageConversions.RgbToByteArray(template_image);
		
		// extract the points
		int[][] cyl_coordinates = new int[3][Sphere.numberOfVoxInSphere((int)Math.ceil(r))]; 
		
		int cnt = this.extractCoords(template_image, cyl_coordinates);
		
		if(cnt>0){
			
			// change the spots, assign them with valueToDraw
			for (int i = 0; i < cnt; i++) {
				
				int row 	= cyl_coordinates[0][i];
				int col 	= cyl_coordinates[1][i];
				int layer 	= cyl_coordinates[2][i];
				
				img_array[0][layer][ArrayHandling.sub2index_2d(row, col, w)] = rgb_values[0];
				img_array[1][layer][ArrayHandling.sub2index_2d(row, col, w)] = rgb_values[1];
				img_array[2][layer][ArrayHandling.sub2index_2d(row, col, w)] = rgb_values[2];
			}
			
		}else{
			
			System.out.println("Sphere:drawOverImage()\nThere was not enough points to draw anything!");
		
		}
		
		for (int i = 1; i <= template_image.getStack().getSize(); i++) {
			
			((ColorProcessor)template_image.getStack().getProcessor(i)).setRGB(
					img_array[0][i-1], 
					img_array[1][i-1], 
					img_array[2][i-1]
							);
		}

	}
	
	public void 		drawOverColorImage(ImagePlus template_image, int r, int g, int b){
		
		int w = template_image.getWidth();
		
		if((template_image.getType()!=ImagePlus.COLOR_RGB)){
			template_image = ImageConversions.ImagePlusToRGB(template_image);
			System.out.println("converting ImagePlus to rgb...");
		}

		byte[] rgb_values = new byte[]{(byte)r, (byte)g, (byte)b}; 
		
		byte[][][] img_array = ImageConversions.RgbToByteArray(template_image);
		
		// extract the points
		int[][] cyl_coordinates = new int[3][Sphere.numberOfVoxInSphere((int)Math.ceil(r))]; 
		
		int cnt = extractCoords(template_image, cyl_coordinates);
		
		if(cnt>0){
			
			// change the spots, assign them with valueToDraw
			for (int i = 0; i < cnt; i++) {
				
				int row 	= cyl_coordinates[0][i];
				int col 	= cyl_coordinates[1][i];
				int layer 	= cyl_coordinates[2][i];
				
				img_array[0][layer][ArrayHandling.sub2index_2d(row, col, w)] = rgb_values[0];
				img_array[1][layer][ArrayHandling.sub2index_2d(row, col, w)] = rgb_values[1];
				img_array[2][layer][ArrayHandling.sub2index_2d(row, col, w)] = rgb_values[2];
			}
			
		}else{
			
			System.out.println("Sphere:drawOverImage()\nThere was not enough points to draw anything!");
		
		}
		
		for (int i = 1; i <= template_image.getStack().getSize(); i++) {
			
			((ColorProcessor)template_image.getStack().getProcessor(i)).setRGB(
					img_array[0][i-1], 
					img_array[1][i-1], 
					img_array[2][i-1]
							);
		}

	}
	
	public static void 	drawOverColorImage(Sphere[] sph, ImagePlus template_image, int valueToDraw){
		for (int i = 0; i < sph.length; i++) {
			sph[i].drawOverColorImage(template_image, valueToDraw);
		}
	}
	
	/*	##################################
	 *  AUX
	 * 	##################################
	 */	
	
	public int 			numberOfVoxInSphere(){
		
		return Sphere.numberOfVoxInSphere((int)Math.ceil(r));

	}
	
	public static int 	numberOfVoxInSphere(int radius){
		int voxelNr = 0;
		switch (radius) {
		case 1: voxelNr = 7      ; break; 
		case 2: voxelNr = 33     ; break; 
		case 3: voxelNr = 123    ; break; 
		case 4: voxelNr = 257    ; break; 
		case 5: voxelNr = 515    ; break; 
		case 6: voxelNr = 925    ; break; 
		case 7: voxelNr = 1419   ; break; 
		case 8: voxelNr = 2109   ; break; 
		case 9: voxelNr = 3071   ; break; 
		case 10: voxelNr = 4169  ; break; 
		case 11: voxelNr = 5575  ; break; 
		case 12: voxelNr = 7153  ; break; 
		case 13: voxelNr = 9171  ; break; 
		case 14: voxelNr = 11513 ; break; 
		case 15: voxelNr = 14147 ; break; 
		case 16: voxelNr = 17077 ; break; 
		case 17: voxelNr = 20479 ; break; 
		case 18: voxelNr = 24405 ; break; 
		case 19: voxelNr = 28671 ; break; 
		case 20: voxelNr = 33401 ; break; 
        default: voxelNr = 33401 ; break;
		}
		
		return voxelNr;
		
	}

}