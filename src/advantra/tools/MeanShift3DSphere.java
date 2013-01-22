package advantra.tools;

import ij.ImagePlus;
import ij.gui.NewImage;
import advantra.general.ArrayHandling;
import advantra.general.DebugExport;
import advantra.general.Transf;
import advantra.processing.IntensityCalc;
import advantra.shapes.Sphere;

public class MeanShift3DSphere {

	/*
	 *  reference: 
	 *  Mean shift, mode seeking, and clustering, Yizong Cheng 
	 *  doi: 10.1109/34.400568
	 */
	
	ImagePlus	sphere_img;				// will be used as source for weights
	
	double		angle_range;
	Sphere		sphere_roi;				// sphere that captures the image
	
	double[][]	S;						
										/* 
										 * Nx2 
										 * 1. col inclination (phi)
										 * 2. col azimuth (theta)
										 */
	
	double[][]	T;						
										/*
										 * 2xN
										 * 1. col phi
										 * 2. col theta
										 */
	
	public MeanShift3DSphere(ImagePlus sphere_img, double angle_range, int N){
		
		this.sphere_img				= sphere_img;
		this.angle_range 			= angle_range;
		
		double mid_x				= (double)sphere_img.getHeight()/2 - 0.5;
		double mid_y				= (double)sphere_img.getWidth()/2  - 0.5;
		double mid_z				= (double)sphere_img.getStack().getSize()/2 - 0.5;
		double mid_r 				= mid_x; 
		
		// form the sphere to extract the starting directions - matrix S
		sphere_roi = new Sphere(mid_x, mid_y, mid_z, mid_r);
				
		//double[][] spherical_coords = sphere_roi.extractSphericalCoords(this.sphere_img); 
		double[][] cartesian_coords = sphere_roi.generate3DFullSpherePts(N);
		double[][] spherical_coords = new double[N][3];
		for (int i = 0; i < N; i++) {
			Transf.cart2sph(
					cartesian_coords[0][i]-sphere_roi.getCenterX(), 
					cartesian_coords[1][i]-sphere_roi.getCenterY(), 
					cartesian_coords[2][i]-sphere_roi.getCenterZ(), 
					spherical_coords[i]);
		}
		
		// initialize S, T will stay zero
		S = new double[N][2];
		T = new double[N][2];
		for (int i = 0; i < N; i++) {
			S[i][0] = spherical_coords[i][1];	S[i][1] = spherical_coords[i][2];
		}
		
	}
	
	public int 			getNumberOfPoints(){
		return S.length;
	}
	
	public double[][] 	getT(){
		return T;
	}

	public double[][] 	getS(){
		return S;
	}
	
	private double[] 	runOneIterAtPos(double[] curr_pos){ 
		// curr_pos[0] = phi, curr_pos[1] = theta
		double[] 	new_pos			= new double[2];
		double 		sum 			= 0;
		
		// extract surrounding directions in both cartesian & spherical - this is mean-shift 'sphere'
//		Sphere temp_sph = new Sphere(
//				this.sphere_roi.getCenterX(), 
//				this.sphere_roi.getCenterY(), 
//				this.sphere_roi.getCenterZ(), 
//				r);
		
		double[][] solid_angle_cart = 
				sphere_roi.coordinatesSolidAngleSpherePts(angle_range, 10, curr_pos[0], curr_pos[1]); // 10 steps in angle_range
		double[][] solid_angle_sph = sphere_roi.cartesian2spherical(solid_angle_cart);
		
//		System.out.println("curr_pos (cart):");
//		double[] x_y_z_around = new double[3];
//		Transf.sph2cart(r, curr_pos[0], curr_pos[1], x_y_z_around);
//		x_y_z_around[0]+=this.sphere_roi.getCenterX();
//		x_y_z_around[1]+=this.sphere_roi.getCenterY();
//		x_y_z_around[2]+=this.sphere_roi.getCenterZ();
//		ArrayHandling.print1DArray(x_y_z_around);
//		System.out.println("curr_pos (sph):");
//		System.out.println(r+" , "+curr_pos[0]+" , "+curr_pos[1]);
//		System.out.println("sph:");
//		ArrayHandling.print2DArray(solid_angle_sph);
//		System.out.println("cart:");
//		ArrayHandling.print2DArray(solid_angle_cart);
		
		// here it comes down to a problem - weighting angles when they're not wrapped is wrong
		// to avoid it, I am doing weighted mean of wrapped differences in angles phi & theta so that 
		// the value of the mean-shift vector is added on current value without risk of getting 
		// corrupted weighted means
		
		// take the interpolated values from the image to use as weights
		//IntensityCalc calculator = new IntensityCalc(sphere_img.getStack());
		
//		System.out.format("extracted\n");
		IntensityCalc calc = new IntensityCalc(sphere_img.getStack());
		for (int i = 0; i < solid_angle_sph.length; i++) {
			
			//double weight 	=  sphere_roi.avgValuePerDirection(solid_angle_sph[i][1], solid_angle_sph[i][2], sphere_img.getStack());
			float weight 	= 	sphere_roi.sphereSurfaceValuePerDirection(solid_angle_sph[i][1], solid_angle_sph[i][2], calc);		
//				calculator.interpolateAt(
//				solid_angle_cart[i][0], 
//				solid_angle_cart[i][1], 
//				solid_angle_cart[i][2]);
//			System.out.format("%f, %f, %f , %f \n", 
//					solid_angle_cart[i][0], 
//					solid_angle_cart[i][1], 
//					solid_angle_cart[i][2],
//					weight);			
			
			new_pos[0] 		+= weight * Transf.angle_wrap_pi(solid_angle_sph[i][1]-curr_pos[0]);	
			new_pos[1] 		+= weight * Transf.angle_wrap_pi(solid_angle_sph[i][2]-curr_pos[1]);
			sum				+= weight;
		}
		
		if(sum>0){
			new_pos[0] = new_pos[0]/sum + curr_pos[0]; 	// phi
			new_pos[1] = new_pos[1]/sum + curr_pos[1];	// theta
			// wrap azimuth only
			new_pos[1] = Transf.angle_wrap_2_pi(new_pos[1]);
		}
		else {// keep the old
			new_pos[0] = curr_pos[0];
			new_pos[1] = curr_pos[1];
		}
		
//		System.out.println("new_pos (sph):");
//		System.out.println(r+" , "+new_pos[0]+" , "+new_pos[1]);
		
		return new_pos;
		
	}
	
	public double[][] 	run(int max_iter, double epsilon){
		
		T				= new double[S.length][2];
		
		for (int i = 0; i < T.length; i++) {
			T[i][0] = S[i][0];	T[i][1] = S[i][1];
		}
		
		for (int i = 0; i < T.length; i++) { // 
			
			//long t11 = System.currentTimeMillis();
			int iter = 0;
			double d = Double.MAX_VALUE;
			
			do{
				
				double[] new_pos = runOneIterAtPos(T[i]);
				
				d = d2(new_pos, T[i]);
				T[i][0] = new_pos[0];	T[i][1] = new_pos[1];
				
				iter++;
				
			}
			while(iter < max_iter && d > epsilon*epsilon);
			//long t22 = System.currentTimeMillis();
			//System.out.format("point %d: %d iterations, epsilon=%f, elapsed %f sec\n", i, iter, epsilon, ((t22-t11)/1000f));
			
		}
		
		//saveConvergencePts(String.format("T_iter%04d.csv", max_iter));
		//saveConvergencePtsAsImage(String.format("T_iter%04d.tif", max_iter));
		return T;
		
	}
	
	public double[][] extractDirections(double range, int M){ 
		
		// 'range' describes neighborhood range size, 
		// 'M' number of samples within the cluster
		
		// T will be checked and extracted clusters stored in it, starting from the
		// element with first index

		boolean[] 	clustered 		= new boolean[T.length]; // initialized with false
		int[] 		cluster_size 	= new int[T.length];
		int nr_clusters = 0;
		int nr_clusters_nigher_than_M = 0;
		
		for (int i = 0; i < T.length; i++) {
			
			if(!clustered[i]){
				
				clustered[i] = true;
				cluster_size[nr_clusters]++;
				T[nr_clusters][0] = T[i][0];
				T[nr_clusters][1] = T[i][1];
				
				for (int j = i+1; j < T.length; j++) {
					
					if(!clustered[j] && d(T[i], T[j])<=range){
						
						clustered[j] = true;
						cluster_size[nr_clusters]++;
						
					}
					
				}
				
				if(cluster_size[nr_clusters]>M){
					nr_clusters_nigher_than_M++;
				}
				
				nr_clusters++;
				
			}
			
		}
		
		// take them out
		double[][] out_directions = null;
		int cnt = 0;
		
		out_directions = new double[nr_clusters_nigher_than_M][3];
			
		for (int i = 0; i < nr_clusters; i++) {
			if(cluster_size[i]>M){
				Transf.sph2cart(1.0, T[i][0], T[i][1], out_directions[cnt]); 
				cnt++;
			}
		}
			
		return out_directions;
		
	}
	
	private double d2(double[] a, double[] b){
		return Math.pow(a[0]-b[0], 2)+Math.pow(a[1]-b[1], 2);
	}
	
	private double d(double[] a, double[] b){
		// phi distance
		double phi_dist = a[0]-b[0];
		// wrap azimuth
		double az_dist = a[1]-b[1];
		double wrapped_az = Transf.angle_wrap_pi(az_dist);
		return Math.sqrt(Math.pow(phi_dist, 2)+Math.pow(wrapped_az, 2));
	}
	
	public void saveS(String csv_file_name){
		DebugExport f = new DebugExport(csv_file_name);
		
		for (int i = 0; i < S.length; i++) {
			double[] x_y_z = new double[3];
			Transf.sph2cart(sphere_roi.getR(), S[i][0], S[i][1], x_y_z);
			
			x_y_z[0] += this.sphere_roi.getCenterX();
			x_y_z[1] += this.sphere_roi.getCenterY();
			x_y_z[2] += this.sphere_roi.getCenterZ();
			
			f.writeLine(String.format("%f, %f, %f", x_y_z[0], x_y_z[1], x_y_z[2]));
		}
		f.closeDebug();
		System.out.println(csv_file_name+" saved...");
	}
	
	public void saveT(String csv_file_name){
		DebugExport f = new DebugExport(csv_file_name);
		
		for (int i = 0; i < T.length; i++) {
			double[] x_y_z = new double[3];
			Transf.sph2cart(sphere_roi.getR(), T[i][0], T[i][1], x_y_z);
			
			x_y_z[0] += this.sphere_roi.getCenterX();
			x_y_z[1] += this.sphere_roi.getCenterY();
			x_y_z[2] += this.sphere_roi.getCenterZ();
			
			f.writeLine(String.format("%f, %f, %f", x_y_z[0], x_y_z[1], x_y_z[2]));
		}
		f.closeDebug();
		System.out.println(csv_file_name+" saved...");
	}	
	
	public ImagePlus extractConvergence(int resolution){
		ImagePlus output 		= NewImage.createByteImage(
				"average_values_per_direction", 2*resolution, resolution, 1, NewImage.FILL_BLACK);
		// generate cartesian points on the sphere
		//double[][] pts = sphere_roi.generate3DFullSpherePts_Nx3(T.length);
		// make them spherical
		//double[][] pts_sph = sphere_roi.cartesian2spherical(pts);
		for (int i = 0; i < T.length; i++) {
			// get indexes to be plotted on the image
			double current_phi 		= T[i][0];
			double current_theta 	= T[i][1];
			int row = ArrayHandling.value2index(current_phi, 	ArrayHandling.IdxMode.LAST_INCLUDED, 0.0, Math.PI, 	resolution);
			int col = ArrayHandling.value2index(current_theta, 	ArrayHandling.IdxMode.LAST_EXCLUDED, 0.0, 2*Math.PI, 2*resolution);
			double g = sphere_roi.avgValuePerDirection(current_phi, current_theta, sphere_img.getStack());
			output.getStack().setVoxel(col, row, 0, g);
		}
		
		return output;

	}
	
}

//private boolean isNeighbour(double[] a, double[] b, double range){
//return d(a, b)<=range;
//}
//public void saveStartPtsAsImage(String img_file_name){
//ImagePlus img = this.sphere_img_rgb.duplicate();
//
//for (int i = 0; i < S.length; i++) {
//	double[] x_y_z = new double[3];
//	Transf.sph2cart(r[i], S[i][0], S[i][1], x_y_z);
//	x_y_z[0] += this.sphere_roi.getCenterX();
//	x_y_z[1] += this.sphere_roi.getCenterY();
//	x_y_z[2] += this.sphere_roi.getCenterZ();
//	(new Point(x_y_z)).drawOverColorImageStack(img, 255, 0, 0);
//}
//(new FileSaver(img)).saveAsTiffStack(img_file_name);
//System.out.println(img_file_name+" saved...");
//}
//public void saveConvergencePtsAsImage(String img_file_name){
//
//long timeInMillis = System.currentTimeMillis();
//
//
//ImagePlus img = this.sphere_img_rgb.duplicate();//new ImagePlus("conv", this.sphere_img_rgb.getStack());
//
//for (int i = 0; i < T.length; i++) {
//	double[] x_y_z = new double[3];
//	Transf.sph2cart(r[i], T[i][0], T[i][1], x_y_z);
//	
//	x_y_z[0] += this.sphere_roi.getCenterX();
//	x_y_z[1] += this.sphere_roi.getCenterY();
//	x_y_z[2] += this.sphere_roi.getCenterZ();
//	
//	(new Point(x_y_z)).drawOverColorImageStack(img, 0, 255, 0);
//}
//(new FileSaver(img)).saveAsTiffStack(img_file_name);
//
//System.out.println(img_file_name+" saved...  "+(System.currentTimeMillis()-timeInMillis)+" ms");
//}

//private void saveOutDirectionsAsImage(String img_file_name, double[][] directions){
//
//long timeInMillis = System.currentTimeMillis();
//
//
//ImagePlus img = this.sphere_img_rgb.duplicate();
//
//for (int i = 0; i < directions.length; i++) {
//	double X = this.sphere_roi.getCenterX()+this.sphere_roi.getR()*directions[i][0];
//	double Y = this.sphere_roi.getCenterY()+this.sphere_roi.getR()*directions[i][1];
//	double Z = this.sphere_roi.getCenterZ()+this.sphere_roi.getR()*directions[i][2];
//	(new Sphere(X, Y, Z, 3)).drawOverColorImage(img, 255, 0, 0);
//}
//(new FileSaver(img)).saveAsTiffStack(img_file_name);
//
//System.out.println(img_file_name+" saved...  "+(System.currentTimeMillis()-timeInMillis)+" ms");
//}