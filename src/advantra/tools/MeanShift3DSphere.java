package advantra.tools;

import java.awt.Color;
import java.util.Vector;

import ij.ImagePlus;
import ij.gui.Plot;
import ij.io.FileSaver;
import advantra.general.Transf;
import advantra.processing.IntensityCalc;
import advantra.shapes.Sphere;

public class MeanShift3DSphere extends Thread {

	/*
	 *  reference: 
	 *  Mean shift, mode seeking, and clustering, Yizong Cheng 
	 *  doi: 10.1109/34.400568
	 */
	
	private static double		angle_range;			// main mean-shift parameter
	private static Sphere		sphere_space;			// original sphere when taken from the image
	private static Sphere		sphere_roi;				// sphere that captures the image
	
	private static double[][]	S;						
										/* 
										 * Nx2 
										 * 1. col inclination (phi)
										 * 2. col azimuth (theta)
										 */
	
	private static double[][]	T;						
										/*
										 * Nx2
										 * 1. col phi
										 * 2. col theta
										 */
	
	public static double[][]	T_clust;
										/*
										 * extract clusters of T
										 * Nx2
										 * 1. col phi
										 * 2. col theta
										 */
	
	public static double[][]	cluster_seed;			// global image stack coordinates
	//int[][]		cluster_local_seeds;				// coordinates in local sphere image coordinates (need it for connectivity tests)
	
										/*
										 * NrClustx3
										 * 1. col: row coordinate (x)
										 * 2. col: col coordinate (y)
										 * 3. col: lay coordinate (z)
										 */
	
	private static int max_iter; 
	private static double epsilon;
	private static IntensityCalc img_calc;
	
	/*
	 * non-static elements define parallel execution
	 */
	
	private int n0, n1; // currently these apply to run()
	
	public MeanShift3DSphere(int n0, int n1){
		this.n0 = n0;
		this.n1 = n1;
	}
	
	public static void load(ImagePlus sphere_img, Sphere sphere_space1, double angle_range1, int N, int max_itr, double eps){
		
//		new FileSaver(sphere_img).saveAsTiffStack("ms_sphere.tif");
		
		img_calc 				= new IntensityCalc(sphere_img.getStack());
		angle_range 			= angle_range1;
		sphere_space			= sphere_space1;//new Sphere(sphere_space);
		
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
		
		epsilon = eps;
		max_iter = max_itr;
		
		// initialize S, T will stay zero, T_clust and cluster_* are null before anything is calculated
		S = new double[N][2];
		T = new double[N][2];
		
		for (int i = 0; i < N; i++) {
			S[i][0] = spherical_coords[i][1];	S[i][1] = spherical_coords[i][2];
		}
		
	}
	
	public static int 			getNumberOfPoints(){
		return S.length;
	}

	private float 				avgValuePerDirection(double phi, double theta){
		
		float value 	= 0;
		
		int count 		= 0;
		
		for (double r = 0.90*sphere_roi.getR(); r <= sphere_roi.getR(); r+=0.03*sphere_roi.getR()) {
			//double r = 0.95*sphere_roi.getR();
			double x = sphere_roi.getCenterX()+Transf.sph2cart_x(r, phi, theta);
			double y = sphere_roi.getCenterY()+Transf.sph2cart_y(r, phi, theta);
			double z = sphere_roi.getCenterZ()+Transf.sph2cart_z(r, phi, theta);
			value += img_calc.interpolateAt((float)x, (float)y, (float)z);
			count++;
		}
		
		return value/(float)count;
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
			
			float weight 	= 	
					avgValuePerDirection(solid_angle_sph[i][1], solid_angle_sph[i][2]);		
			
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
		
		return new_pos;
		
	}

	// this one runs in parallel
	public void 	run(){
		
		
		
		for (int i = n0; i < n1; i++) {
			
			T[i][0] = S[i][0];	T[i][1] = S[i][1];
			
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
		
	}
	
	public static void 			extractClusters(double range, int M){ 
		
		// 'range' describes neighborhood range size, 
		// 'M' number of samples within the cluster
		// T will be checked and extracted clusters stored in it, starting from the
		// element with first index, results stored in T_clust

		boolean[] 	seed_clustered 		= new boolean[T.length]; 
		int[] 		seed_cluster_idx 	= new int[T.length];
		Vector<Integer> cluster_sizes	= new Vector<Integer>();
		int nr_clusters 				= 0;
		
		for (int i = 0; i < T.length; i++) {
			for (int j = 0; j < T.length; j++) {
				
				if(!seed_clustered[j] && j!=i){ 
					
					
					if(d(T[i], T[j])<range){
						
						if(seed_clustered[i]){
							// i was clustered
							seed_clustered[j] = true;
							seed_cluster_idx[j] = seed_cluster_idx[i];
							
							int tmp = cluster_sizes.get(seed_cluster_idx[i]);
							tmp++;
							cluster_sizes.setElementAt(tmp, seed_cluster_idx[i]);
							
						}
						else{
							// i was not clustered
							seed_clustered[i] = seed_clustered[j] = true;
							seed_cluster_idx[i] = seed_cluster_idx[j] = nr_clusters;
							cluster_sizes.add(2);
							nr_clusters++;
						}
					}
				}
			}
			
			if(!seed_clustered[i]){
				seed_clustered[i] = true;
				seed_cluster_idx[i] = nr_clusters;
				cluster_sizes.add(1);//set(nr_clusters, 1);
				nr_clusters++;
			}
			
		}
		
//		System.out.println("\ncluster sizes for "+T.length+" conv. pts, M thr. "+M+" and range "+range);
//		for (int k = 0; k < cluster_sizes.size(); k++) {
//			System.out.print("c["+k+"] = "+cluster_sizes.get(k)+" el. ");
//		}
		
		int nr_clusters_higher_than_M = 0;
		// so that it is known how much to allocate
		for (int i = 0; i < cluster_sizes.size(); i++) {
			if(cluster_sizes.get(i)>M){
				nr_clusters_higher_than_M++;
			}
		}
		
		if(nr_clusters_higher_than_M<=0){
			T_clust 			= null;
			cluster_seed 		= null;
			return;
		}
		
		T_clust 		= new double[nr_clusters_higher_than_M][2];
		cluster_seed 	= new double[nr_clusters_higher_than_M][3];
		
		int cnt = 0;
		for (int i = 0; i < cluster_sizes.size(); i++) {
			if(cluster_sizes.get(i)>M){
				
				// i marks the cluster index
				// take the first one with this cluster index because they are fairly close
				int take_one = -1;
				for (int k = 0; k < seed_cluster_idx.length; k++) {
					if(seed_cluster_idx[k]==i){
						take_one = k;
					}
				}
				
				T_clust[cnt][0] = T[take_one][0]; // not necessary to keep T_clust
				T_clust[cnt][1] = T[take_one][1];
				
				Transf.sph2cart(1.0*sphere_space.getR(), T_clust[cnt][0], T_clust[cnt][1], cluster_seed[cnt]);
				
				cluster_seed[cnt][0] += sphere_space.getCenterX();
				cluster_seed[cnt][1] += sphere_space.getCenterY();
				cluster_seed[cnt][2] += sphere_space.getCenterZ();
				cnt++;
			}
		}
		
		// plot for debug
		
//		if(nr_clusters_higher_than_M!=2){
//		
//		float[] x0 = new float[T.length];
//		float[] y0 = new float[T.length];
//		
//		float[] x1 = new float[T.length];
//		float[] y1 = new float[T.length];
//		
//		for (int i = 0; i < T.length; i++) {
//			
//			x0[i] = (float)S[i][1];
//			y0[i] = (float)S[i][0];
//			
//			x1[i] = (float)T[i][1];
//			y1[i] = (float)T[i][0];
//			
//		}
//		
//		float[] x2 = new float[nr_clusters_higher_than_M];
//		float[] y2 = new float[nr_clusters_higher_than_M];
//		
//		for (int i = 0; i < nr_clusters_higher_than_M; i++) {
//			x2[i] = (float)T_clust[i][1];
//			y2[i] = (float)T_clust[i][0];
//		}
//		
//		Plot plot = new Plot(String.format("Mean-shift,range %f",range), "phi", "theta", new double[0], new double[0]);
//	    plot.setFrameSize(600,300);
//	    plot.setLimits(0, 2*Math.PI, 0, 1*Math.PI);
//	    plot.setLineWidth(2);
//	    plot.setColor(Color.gray);
//	    plot.addPoints(x0, y0, Plot.DOT);
//	    plot.setColor(Color.blue);
//	    plot.addPoints(x1, y1, Plot.X);
//	    plot.setColor(Color.red);
//	    plot.addPoints(x2, y2, Plot.BOX);
//	    plot.show(); 
//
//	    printS();
//	    printT();
//	    printConv();
//	    
//	    
//		}

		
	}
	
	private static double 		d2(double[] a, double[] b){
		return Math.pow(a[0]-b[0], 2)+Math.pow(a[1]-b[1], 2);
	}
	
	public static double 		d(double[] a, double[] b){
		
		double phi1 	= a[0];
		double theta1 	= a[1];
		
		double phi2 	= b[0];
		double theta2 	= b[1];
		
		// phi inclination
		double dphi = Math.abs(phi2-phi1);
		// theta azimuth
		if(theta1>theta2){
			double dummy = theta1; theta1 = theta2; theta2 = dummy;
		}
		double dtheta1 = theta2 - theta1;
		double dtheta2 = theta1 + 2*Math.PI - theta2;
		
		double dtheta = (dtheta1<dtheta2)? dtheta1 : dtheta2 ;

		return (dtheta>dphi)? dtheta : dphi ;
	}

	/*
	 * 
	 * for plotting (no other functionality)
	 * 
	 */
	
	public static void 	printT(){
		System.out.println("T[][] in cartesian");
		double[][] T_cartesian = new double[T.length][3];
		for (int i = 0; i < T.length; i++) {
			Transf.sph2cart(sphere_roi.getR(), T[i][0], T[i][1], T_cartesian[i]);
			T_cartesian[i][0] += sphere_roi.getCenterX()+1;
			T_cartesian[i][1] += sphere_roi.getCenterY()+1;
			T_cartesian[i][2] += sphere_roi.getCenterZ()+1;
			System.out.format("%.2f, %.2f, %.2f \n", T_cartesian[i][1], T_cartesian[i][0], T_cartesian[i][2]);
		}
		
	}

	public static void	printS(){
		System.out.println("S[][] in cartesian");
		double[][] S_cartesian = new double[S.length][3];
		for (int i = 0; i < S.length; i++) {
			Transf.sph2cart(sphere_roi.getR(), S[i][0], S[i][1], S_cartesian[i]);
			S_cartesian[i][0] += sphere_roi.getCenterX()+1;
			S_cartesian[i][1] += sphere_roi.getCenterY()+1;
			S_cartesian[i][2] += sphere_roi.getCenterZ()+1;
			System.out.format("%.2f, %.2f, %.2f \n", S_cartesian[i][1], S_cartesian[i][0], S_cartesian[i][2]);
		}
		
	}

	public static void	printConv(){
		System.out.println("Convergence[][] in cartesian");
		double[][] T_clust_cartesian = new double[T_clust.length][3];
		for (int i = 0; i < T_clust.length; i++) {
			Transf.sph2cart(sphere_roi.getR(), T_clust[i][0], T_clust[i][1], T_clust_cartesian[i]);
			T_clust_cartesian[i][0] += sphere_roi.getCenterX()+1;
			T_clust_cartesian[i][1] += sphere_roi.getCenterY()+1;
			T_clust_cartesian[i][2] += sphere_roi.getCenterZ()+1;
			System.out.format("%.2f, %.2f, %.2f \n", T_clust_cartesian[i][1], T_clust_cartesian[i][0], T_clust_cartesian[i][2]);
		}
		
	}
	
}


/*
private static double 		angle2(double[] T1, double[] T2){
	double[] a = new double[3];
	double[] b = new double[3];
	Transf.sph2cart(1.0, T1[0], T1[1], a);
	Transf.sph2cart(1.0, T2[0], T2[1], b);
	return angle3(a, b);
	
}
*/
/*
private static double 		angle3(double[] a, double[] b){
	
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
*/
/*
private static double 		dotProd3(double[] a, double[] b){
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
*/
//public void 		extractClusters_new(double angle_th, int M){
//
////System.out.println("single linkage");
//
////long t11 = System.currentTimeMillis();
//
//double[][] D 				= new double[T.length][T.length];
//double[] cartesian_i 		= new double[3];
//double[] cartesian_j 		= new double[3];
//
//// fill the distances - thresholded distance map
//fextractClusters_newor (int i = 0; i < T.length; i++) {
//	
//	Transf.sph2cart(1.0, T[i][0], T[i][1], cartesian_i);
//	
//	for (int j = i; j < T.length; j++) {
//		if(i==j){
//			D[i][j] = 0;
//		}
//		else{
//			// ang between i and j
//			Transf.sph2cart(1.0, T[j][0], T[j][1], cartesian_j);
//			D[i][j] = angle3(cartesian_i, cartesian_j);// in radians
//			if(D[i][j]<=angle_th){
//				D[i][j] = 0;
//			}
//			D[j][i] = D[i][j];
//			
//		}
//		
//	}
//}
//
////System.out.println("table calculated.");
//
//// group clusters
//
//Vector<Vector<Integer>> clusters = new Vector<Vector<Integer>>();
//
//boolean[] clustered = new boolean[T.length];
//for (int i = 0; i < clustered.length; i++) {
//	clustered[i] = false;
//}
//
////int nr_classes = 0;
//
//for (int i = 0; i < T.length; i++) { // T.length
//	
//	//System.out.print(".");
//	
//	if(!clustered[i]){
//		
//		clustered[i] = true;
//		
//		Vector<Integer> v = new Vector<Integer>();
//		v.add(i);
//		clusters.add(v);
//		
//		Vector<Integer> neighbors = new Vector<Integer>();
//		
//		// check neighbors for i -th
//		for (int j = 0; j < T.length; j++) {
//			if(j!=i && !clustered[j]){
//				if(D[i][j]==0){
//					neighbors.add(j);
//					clustered[j] = true;
//					//System.out.print("nb:"+j+", ");
//					// each time it was zero - ad it to the list
//					clusters.get(clusters.size()-1).add(j);
//				}
//			}
//		}
//		
//		//System.out.println("nb list: "+neighbors.size());
//		
//		while(neighbors.size()>0){
//			
//			int take_neighbor = neighbors.get(0);
//			neighbors.removeElementAt(0);
//			
//			// check neighbors for take_neighbor -th
//			
//			//if(clustered[]){
//			for (int j = 0; j < T.length; j++) {
//				if (j!=take_neighbor && !clustered[j]) {
//					if(D[take_neighbor][j]==0){
//						neighbors.add(j);
//						clustered[j] = true;
//						// each time it was zero - ad it to the list
//						clusters.get(clusters.size()-1).add(j);
//					}
//				}
//			}
//			
//			//System.out.println("nb list: "+neighbors.size());
//			
//		}
//		
//		//nr_classes++;
//		
//	}
//}
////long t22 = System.currentTimeMillis();
//// return clusters as 'cluster_seed' (used further)
////System.out.println("elapsed "+((t22-t11)/1000f)+" s.");
//
//// : this is a bit stupid solution (two loops just to see the size), cluster_seed should be ArrayList!
//
//int nr_clusters_higher_than_M = 0;
//
//for (int i = 0; i < clusters.size(); i++) {
//	if(clusters.get(i).size()>M){
//		nr_clusters_higher_than_M++;
//	}
//}
//
//cluster_seed_test 	= new double[nr_clusters_higher_than_M][3];
//double[] current_v 	= new double[3]; // cartesian
//double[] avg_v 		= new double[3]; // cartesian
//int idx = 0;
//
//for (int i = 0; i < clusters.size(); i++) {
//	if(clusters.get(i).size()>M){
//		
//		
//		
//		
//		
//		
////		int sz = clusters.get(i).size();
////		// avg
////		avg_v[0] = avg_v[1] = avg_v[2] = 0;
////		for (int j = 0; j < sz; j++) {
////			int index_v = clusters.get(i).get(j);
////			Transf.sph2cart(1.0, T[index_v][0], T[index_v][1], current_v);
////			avg_v[0] += current_v[0];
////			avg_v[1] += current_v[1];
////			avg_v[2] += current_v[2];
////		}
////		
////		double avg_v_norm = Math.sqrt(avg_v[0]*avg_v[0]+avg_v[1]*avg_v[1]+avg_v[2]*avg_v[2]);
////		
////		avg_v[0] /= avg_v_norm;
////		avg_v[1] /= avg_v_norm;
////		avg_v[2] /= avg_v_norm;
////		
////		avg_v[0] *= sphere_roi.getR();
////		avg_v[1] *= sphere_roi.getR();
////		avg_v[2] *= sphere_roi.getR();
//		
//		int index_v = clusters.get(i).get(0);
//		Transf.sph2cart(sphere_roi.getR(), T[index_v][0], T[index_v][1], avg_v);
//		
//		
//		cluster_seed_test[idx][0] = avg_v[0] + sphere_roi.getCenterX();//sphere_space.getCenterX();
//		cluster_seed_test[idx][1] = avg_v[1] + sphere_roi.getCenterY();
//		cluster_seed_test[idx][2] = avg_v[2] + sphere_roi.getCenterZ();
//		
//		idx++;
//	}
//}
//
//}


//public ImagePlus 	extractConvergence(int resolution){
//ImagePlus output 		= NewImage.createByteImage(
//		"average_values_per_direction", 2*resolution, resolution, 1, NewImage.FILL_BLACK);
//for (int i = 0; i < T.length; i++) {
//	// get indexes to be plotted on the image
//	double current_phi 		= T[i][0];
//	double current_theta 	= T[i][1];
//	int row = ArrayHandling.value2index(current_phi, 	ArrayHandling.IdxMode.LAST_INCLUDED, 0.0, Math.PI, 	resolution);
//	int col = ArrayHandling.value2index(current_theta, 	ArrayHandling.IdxMode.LAST_EXCLUDED, 0.0, 2*Math.PI, 2*resolution);
//	double g = avgValuePerDirection(current_phi, current_theta); 
//	output.getStack().setVoxel(col, row, 0, g);
//}
//
//return output;
//
//}


//public void 		saveS(String csv_file_name){
//DebugExport f = new DebugExport(csv_file_name);
//
//for (int i = 0; i < S.length; i++) {
//	double[] x_y_z = new double[3];
//	Transf.sph2cart(sphere_roi.getR(), S[i][0], S[i][1], x_y_z);
//	
//	x_y_z[0] += sphere_roi.getCenterX();
//	x_y_z[1] += sphere_roi.getCenterY();
//	x_y_z[2] += sphere_roi.getCenterZ();
//	
//	f.writeLine(String.format("%f, %f, %f", x_y_z[0], x_y_z[1], x_y_z[2]));
//}
//f.closeDebug();
//System.out.println(csv_file_name+" saved...");
//}

//public void 		saveT(String csv_file_name){
//DebugExport f = new DebugExport(csv_file_name);
//
//for (int i = 0; i < T.length; i++) {
//	double[] x_y_z = new double[3];
//	Transf.sph2cart(sphere_roi.getR(), T[i][0], T[i][1], x_y_z);
//	
//	x_y_z[0] += this.sphere_roi.getCenterX();
//	x_y_z[1] += this.sphere_roi.getCenterY();
//	x_y_z[2] += this.sphere_roi.getCenterZ();
//	
//	f.writeLine(String.format("%f, %f, %f", x_y_z[0], x_y_z[1], x_y_z[2]));
//}
//f.closeDebug();
//System.out.println(csv_file_name+" saved...");
//}	