package advantra.tools;

import java.util.Vector;

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
	
	IntensityCalc	img_calc;			// will be used as tool to calculate weights
	
	double		angle_range;			// main mean-shift parameter
	Sphere		sphere_space;			// original sphere when taken from the image
	Sphere		sphere_roi;				// sphere that captures the image
	
	double[][]	S;						
										/* 
										 * Nx2 
										 * 1. col inclination (phi)
										 * 2. col azimuth (theta)
										 */
	
	double[][]	T;						
										/*
										 * Nx2
										 * 1. col phi
										 * 2. col theta
										 */
	
	double[][]	T_clust;
										/*
										 * extract clusters of T
										 * Nx2
										 * 1. col phi
										 * 2. col theta
										 */
	
	double[][]	cluster_dirs;
	double[][]	cluster_seed;			// global image stack coordinates
	double[][]	cluster_seed_test;
	int[][]		cluster_local_seeds;	// coordinates in local sphere image coordinates (need it for connectivity tests)
	
										/*
										 * NrClustx3
										 * 1. col: row coordinate (x)
										 * 2. col: col coordinate (y)
										 * 3. col: lay coordinate (z)
										 */
	
	public MeanShift3DSphere(ImagePlus sphere_img, Sphere sphere_space, double angle_range, int N){
		
		this.img_calc				= new IntensityCalc(sphere_img.getStack());
		this.angle_range 			= angle_range;
		this.sphere_space			= new Sphere(sphere_space);
		
		double mid_x				= (double)sphere_img.getHeight()/2 - 0.5;
		double mid_y				= (double)sphere_img.getWidth() /2 - 0.5;
		double mid_z				= (double)sphere_img.getStack().getSize()/2 - 0.5;
		double mid_r 				= mid_x; 
		
		// form the sphere to extract the starting directions - matrix S
		sphere_roi = new Sphere(mid_x, mid_y, mid_z, mid_r);
		double[][] cartesian_coords = sphere_roi.generate3DFullSpherePts(N);
		double[][] spherical_coords = new double[N][3];
		for (int i = 0; i < N; i++) {
			Transf.cart2sph(
					cartesian_coords[0][i]-sphere_roi.getCenterX(), 
					cartesian_coords[1][i]-sphere_roi.getCenterY(), 
					cartesian_coords[2][i]-sphere_roi.getCenterZ(), 
					spherical_coords[i]);
		}
		
		// initialize S, T will stay zero, T_clust and cluster_* are null before anything is calculated
		S = new double[N][2];
		T = new double[N][2];
		
		T_clust 				= null;
		cluster_dirs 			= null;
		cluster_seed 			= null;
		cluster_seed_test		= null;
		cluster_local_seeds 	= null;
		
		for (int i = 0; i < N; i++) {
			S[i][0] = spherical_coords[i][1];	S[i][1] = spherical_coords[i][2];
		}
		
	}
	
	public int 			getNumberOfPoints(){
		return S.length;
	}
	
	public double[][] 	getT_cartesian(){
		double[][] T_cartesian = new double[T.length][3];
		for (int i = 0; i < T.length; i++) {
			Transf.sph2cart(sphere_roi.getR(), T[i][0], T[i][1], T_cartesian[i]);
			T_cartesian[i][0] += sphere_roi.getCenterX();
			T_cartesian[i][1] += sphere_roi.getCenterY();
			T_cartesian[i][2] += sphere_roi.getCenterZ();
		}
		return T_cartesian;
	}

	public double[][] 	getClusterSeed(){
		return cluster_seed;
	}
	
	public double[][] 	getClusterSeedTest(){
		return cluster_seed_test;
	}
	
	public int[][] 		getClusterSeedLocal(){
		return cluster_local_seeds;
	}
	
	public double[][] 	getClusterDirs(){
		return cluster_dirs;
	}
	
	public double[][]	getS_cartesian(){
		double[][] S_cartesian = new double[S.length][3];
		for (int i = 0; i < S.length; i++) {
			Transf.sph2cart(sphere_roi.getR(), S[i][0], S[i][1], S_cartesian[i]);
			S_cartesian[i][0] += sphere_roi.getCenterX();
			S_cartesian[i][1] += sphere_roi.getCenterY();
			S_cartesian[i][2] += sphere_roi.getCenterZ();
		}
		return S_cartesian;
	}
	
	public double[][]	getT_Clust_cartesian(){ // row, col, lay
		double[][] T_clust_cartesian = new double[T_clust.length][3];
		for (int i = 0; i < T_clust.length; i++) {
			Transf.sph2cart(sphere_roi.getR(), S[i][0], S[i][1], T_clust_cartesian[i]);
			T_clust_cartesian[i][0] += sphere_roi.getCenterX();
			T_clust_cartesian[i][1] += sphere_roi.getCenterY();
			T_clust_cartesian[i][2] += sphere_roi.getCenterZ();
		}
		return T_clust_cartesian;
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
		
		int nr_d_angle = 10;
		double[][] solid_angle_cart = 
				sphere_roi.coordinatesSolidAngleSpherePts(angle_range, nr_d_angle, curr_pos[0], curr_pos[1]); 
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
		
		for (int i = 0; i < solid_angle_sph.length; i++) {
			
			// avgValuePerDirection or sphereSurfaceValuePerDirection
			float weight 	= 	
					sphere_roi.avgValuePerDirection(solid_angle_sph[i][1], solid_angle_sph[i][2], img_calc);		
			
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
		
		//T				= new double[S.length][2];
		
		for (int i = 0; i < T.length; i++) {
			T[i][0] = S[i][0];	T[i][1] = S[i][1];
		}
		
		//long t11 = System.currentTimeMillis();
		
		for (int i = 0; i < T.length; i++) { // 
			
			
			int iter = 0;
			double d = Double.MAX_VALUE;
			
			do{
				
				double[] new_pos = runOneIterAtPos(T[i]);
				
				d = d2(new_pos, T[i]);
				T[i][0] = new_pos[0];	T[i][1] = new_pos[1];
				
				iter++;
				
			}
			while(iter < max_iter && d > epsilon*epsilon);
			
		}
		//long t22 = System.currentTimeMillis();
		//System.out.format("run():%f sec.\n", ((t22-t11)/1000f));
		
		return T;
		
	}
	
	public ImagePlus 	extractConvergence(int resolution){
		ImagePlus output 		= NewImage.createByteImage(
				"average_values_per_direction", 2*resolution, resolution, 1, NewImage.FILL_BLACK);
		for (int i = 0; i < T.length; i++) {
			// get indexes to be plotted on the image
			double current_phi 		= T[i][0];
			double current_theta 	= T[i][1];
			int row = ArrayHandling.value2index(current_phi, 	ArrayHandling.IdxMode.LAST_INCLUDED, 0.0, Math.PI, 	resolution);
			int col = ArrayHandling.value2index(current_theta, 	ArrayHandling.IdxMode.LAST_EXCLUDED, 0.0, 2*Math.PI, 2*resolution);
			double g = sphere_roi.avgValuePerDirection(current_phi, current_theta, img_calc); 
			output.getStack().setVoxel(col, row, 0, g);
		}
		
		return output;

	}
	
	public void 		extractClusters_new(double angle_th, int M){
		
		//System.out.println("single linkage");
		
		//long t11 = System.currentTimeMillis();
		
		double[][] D 				= new double[T.length][T.length];
		double[] cartesian_i 		= new double[3];
		double[] cartesian_j 		= new double[3];
		
		// fill the distances - thresholded distance map
		for (int i = 0; i < T.length; i++) {
			
			Transf.sph2cart(1.0, T[i][0], T[i][1], cartesian_i);
			
			for (int j = i; j < T.length; j++) {
				if(i==j){
					D[i][j] = 0;
				}
				else{
					// ang between i and j
					Transf.sph2cart(1.0, T[j][0], T[j][1], cartesian_j);
					D[i][j] = angle3(cartesian_i, cartesian_j);// in radians
					if(D[i][j]<=angle_th){
						D[i][j] = 0;
					}
					D[j][i] = D[i][j];
					
				}
				
			}
		}
		
		//System.out.println("table calculated.");
		
		// group clusters
		
		Vector<Vector<Integer>> clusters = new Vector<Vector<Integer>>();
		
		boolean[] clustered = new boolean[T.length];
		for (int i = 0; i < clustered.length; i++) {
			clustered[i] = false;
		}
		
		//int nr_classes = 0;
		
		for (int i = 0; i < T.length; i++) { // T.length
			
			//System.out.print(".");
			
			if(!clustered[i]){
				
				clustered[i] = true;
				
				Vector<Integer> v = new Vector<Integer>();
				v.add(i);
				clusters.add(v);
				
				Vector<Integer> neighbors = new Vector<Integer>();
				
				// check neighbors for i -th
				for (int j = 0; j < T.length; j++) {
					if(j!=i && !clustered[j]){
						if(D[i][j]==0){
							neighbors.add(j);
							clustered[j] = true;
							//System.out.print("nb:"+j+", ");
							// each time it was zero - ad it to the list
							clusters.get(clusters.size()-1).add(j);
						}
					}
				}
				
				//System.out.println("nb list: "+neighbors.size());
				
				while(neighbors.size()>0){
					
					int take_neighbor = neighbors.get(0);
					neighbors.removeElementAt(0);
					
					// check neighbors for take_neighbor -th
					
					//if(clustered[]){
					for (int j = 0; j < T.length; j++) {
						if (j!=take_neighbor && !clustered[j]) {
							if(D[take_neighbor][j]==0){
								neighbors.add(j);
								clustered[j] = true;
								// each time it was zero - ad it to the list
								clusters.get(clusters.size()-1).add(j);
							}
						}
					}
					
					//System.out.println("nb list: "+neighbors.size());
					
				}
				
				//nr_classes++;
				
			}
		}
		//long t22 = System.currentTimeMillis();
		// return clusters as 'cluster_seed' (used further)
		//System.out.println("elapsed "+((t22-t11)/1000f)+" s.");

		// TODO: this is a bit stupid solution (two loops just to see the size), cluster_seed should be ArrayList!
		
		int nr_clusters_higher_than_M = 0;
		
		for (int i = 0; i < clusters.size(); i++) {
			if(clusters.get(i).size()>M){
				nr_clusters_higher_than_M++;
			}
		}
		
		cluster_seed_test 	= new double[nr_clusters_higher_than_M][3];
		double[] current_v 	= new double[3]; // cartesian
		double[] avg_v 		= new double[3]; // cartesian
		int idx = 0;
		
		for (int i = 0; i < clusters.size(); i++) {
			if(clusters.get(i).size()>M){
				
				
				
				
				
				
//				int sz = clusters.get(i).size();
//				// avg
//				avg_v[0] = avg_v[1] = avg_v[2] = 0;
//				for (int j = 0; j < sz; j++) {
//					int index_v = clusters.get(i).get(j);
//					Transf.sph2cart(1.0, T[index_v][0], T[index_v][1], current_v);
//					avg_v[0] += current_v[0];
//					avg_v[1] += current_v[1];
//					avg_v[2] += current_v[2];
//				}
//				
//				double avg_v_norm = Math.sqrt(avg_v[0]*avg_v[0]+avg_v[1]*avg_v[1]+avg_v[2]*avg_v[2]);
//				
//				avg_v[0] /= avg_v_norm;
//				avg_v[1] /= avg_v_norm;
//				avg_v[2] /= avg_v_norm;
//				
//				avg_v[0] *= sphere_roi.getR();
//				avg_v[1] *= sphere_roi.getR();
//				avg_v[2] *= sphere_roi.getR();
				
				int index_v = clusters.get(i).get(0);
				Transf.sph2cart(sphere_roi.getR(), T[index_v][0], T[index_v][1], avg_v);
				
				
				cluster_seed_test[idx][0] = avg_v[0] + sphere_roi.getCenterX();//sphere_space.getCenterX();
				cluster_seed_test[idx][1] = avg_v[1] + sphere_roi.getCenterY();
				cluster_seed_test[idx][2] = avg_v[2] + sphere_roi.getCenterZ();
				
				idx++;
			}
		}
		
	}
	
	public void 		extractClusters(double range, int M){ 
		
		// 'range' describes neighborhood range size, 
		// 'M' number of samples within the cluster
		
		// T will be checked and extracted clusters stored in it, starting from the
		// element with first index, results stored in T_clust

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
				
				for (int j = 0; j < T.length; j++) {
					
					if(!clustered[j] &&  d(T[i], T[j])<=range){  // // angle2(T[i], T[j])
						
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
		
		// take them out and form T_clust
		int cnt = 0;
		if(nr_clusters_nigher_than_M<=0){
			T_clust 			= null;
			cluster_dirs 		= null;
			cluster_seed 		= null;
			cluster_local_seeds = null;
			return;
		}
		
		T_clust 		= new double[nr_clusters_nigher_than_M][2];
		cluster_dirs 	= new double[nr_clusters_nigher_than_M][3];
		cluster_seed 	= new double[nr_clusters_nigher_than_M][3];
		cluster_local_seeds = new int[nr_clusters_nigher_than_M][3];
		
//		System.out.println("total: "+nr_clusters);
//		System.out.println(">M: "+nr_clusters_nigher_than_M);
		
		double[] cartesian_aux = new double[3];
		for (int i = 0; i < nr_clusters; i++) {
			if(cluster_size[i]>M){
				/*
				 */
				T_clust[cnt][0] = T[i][0];
				/*
				 */
				// TODO: expell cluster_dirs
				Transf.sph2cart(1.0, 				T[i][0], T[i][1], cluster_dirs[cnt]); 
				/*
				 */
				Transf.sph2cart(1.1*sphere_space.getR(), 	T[i][0], T[i][1], cluster_seed[cnt]);
				cluster_seed[cnt][0] += sphere_space.getCenterX();
				cluster_seed[cnt][1] += sphere_space.getCenterY();
				cluster_seed[cnt][2] += sphere_space.getCenterZ();
				/*
				 */
				// TODO: expell cluster_local_seeds
				Transf.sph2cart(sphere_roi.getR(), 	T[i][0], T[i][1], cartesian_aux);
				cluster_local_seeds[cnt][0] = (int)Math.round(cartesian_aux[0] + sphere_roi.getCenterX());
				cluster_local_seeds[cnt][1] = (int)Math.round(cartesian_aux[1] + sphere_roi.getCenterY());	
				cluster_local_seeds[cnt][2] = (int)Math.round(cartesian_aux[2] + sphere_roi.getCenterZ());
				
				cnt++;
			}
		}
		
//		System.out.println("cluster_dirs by mean-shift:");
//		ArrayHandling.print2DArray(cluster_dirs);
//		System.out.println("cluster_seed by mean-shift:");
//		ArrayHandling.print2DArray(cluster_seed);
//		System.out.println("cluster_local_seeds by mean-shift:");
//		ArrayHandling.print2DArray(cluster_local_seeds);
		
		
	}
	
	private double 		d2(double[] a, double[] b){
		return Math.pow(a[0]-b[0], 2)+Math.pow(a[1]-b[1], 2);
	}
	
	private double 		d(double[] a, double[] b){
		// phi distance
		double phi_dist = a[0]-b[0];
		// wrap azimuth
		double az_dist = a[1]-b[1];
		double wrapped_az = Transf.angle_wrap_pi(az_dist);
		return (wrapped_az>phi_dist)? wrapped_az : phi_dist ;
		//return Math.sqrt(Math.pow(phi_dist, 2)+Math.pow(wrapped_az, 2));
	}
	
	private double 		angle2(double[] T1, double[] T2){
		double[] a = new double[3];
		double[] b = new double[3];
		Transf.sph2cart(1.0, T1[0], T1[1], a);
		Transf.sph2cart(1.0, T2[0], T2[1], b);
		return angle3(a, b);
		
	}
	
	private double 		angle3(double[] a, double[] b){
		
		double[] a_n  = new double[3];
		double a_norm = Math.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]); 
		if(a_norm>0){
			a_n[0] = a[0]/a_norm;	a_n[1] = a[1]/a_norm;	a_n[2] = a[2]/a_norm;
		}
		else{
			a_n[0] = a[0];	a_n[1] = a[1];	a_n[2] = a[2];
		}
		
		double[] b_n  = new double[3];
		double b_norm = Math.sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]); 
		if(b_norm>0){
			b_n[0] = b[0]/b_norm;	b_n[1] = b[1]/b_norm;	b_n[2] = b[2]/b_norm;
		}
		else{
			b_n[0] = b[0];	b_n[1] = b[1];	b_n[2] = b[2];
		}
		
		double ang = Math.acos(dotProd3(a_n, b_n));
		return ang;
		
	}
	
	private double 		dotProd3(double[] a, double[] b){
		return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
	}
	
	public void 		saveS(String csv_file_name){
		DebugExport f = new DebugExport(csv_file_name);
		
		for (int i = 0; i < S.length; i++) {
			double[] x_y_z = new double[3];
			Transf.sph2cart(sphere_roi.getR(), S[i][0], S[i][1], x_y_z);
			
			x_y_z[0] += sphere_roi.getCenterX();
			x_y_z[1] += sphere_roi.getCenterY();
			x_y_z[2] += sphere_roi.getCenterZ();
			
			f.writeLine(String.format("%f, %f, %f", x_y_z[0], x_y_z[1], x_y_z[2]));
		}
		f.closeDebug();
		System.out.println(csv_file_name+" saved...");
	}
	
	public void 		saveT(String csv_file_name){
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
	
}