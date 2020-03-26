package advantra.tools;

import java.util.ArrayList;

import ij.ImagePlus;

public class MeanShift {
	
	/*
	 *  reference: 
	 *  Mean shift, mode seeking, and clustering, Yizong Cheng 
	 *  doi: 10.1109/34.400568
	 */
	
	int 		image_width;
	int 		image_height;
	
	double[][]	S;					// finite set, data, sample
	double[][]	T;					// cluster centers
	
	double[] 	w; 					// weights (image intensities) 
	
	int 		h_spatial;
	int 		h_spatial2; 		// squared values

	static double Tdist 	= 0.001; 		
	KernelType	knType 	= KernelType.FLAT;
	
	static enum KernelType {
		FLAT, TRUNCATED_GAUSSIAN_1
	}
	
	static enum MeanShiftType {
		NO_BLURRING, BLURRING, NO_BLURRING_WRAPPED_COLUMNS
	}
	
	/*
	 * constructor
	 */
	
	public MeanShift(ImagePlus img, int h_spatial){
		
		/*
		 *  check if it is gray8
		 */
		boolean isGray8 = img!=null && img.getType()==ImagePlus.GRAY8;
		if(!isGray8){
			System.err.println("Image was not gray8. MeanShift class was not constructed...");
			System.exit(1);
		}
		
		image_height 	= img.getHeight();
		image_width 	= img.getWidth();
		
		byte[] image_values = (byte[])img.getProcessor().getPixels();
		
		w				= new double[image_height*image_width];
		for (int i = 0; i < image_height*image_width; i++) {
			w[i] 		= (double)(image_values[i] & 0xff);
		}
		
		int count = 0;
		S				= new double[image_height*image_width][2];
		T				= new double[image_height*image_width][2];
		for (int i = 0; i < image_height; i++) {
			for (int j = 0; j < image_width; j++) {
				S[count][0] = (double)i;
				S[count][1] = (double)j;
				count++;
			}
		}
		
		this.h_spatial 		= h_spatial;
		h_spatial2 			= h_spatial*h_spatial;
		knType 				= KernelType.FLAT;
		
	}
	
	public void set(ImagePlus img, int h_spatial){
		
		/*
		 *  check if it is gray8
		 */
		boolean isGray8 = img!=null && img.getType()==ImagePlus.GRAY8;
		if(!isGray8){
			System.err.println("Image was not gray8. MeanShift class was not constructed...");
			System.exit(1);
		}
		
		image_height 	= img.getHeight();
		image_width 	= img.getWidth();
		
		byte[] image_values = (byte[])img.getProcessor().getPixels();
		
		w				= new double[image_height*image_width];
		for (int i = 0; i < image_height*image_width; i++) {
			w[i] 		= (double)(image_values[i] & 0xff);
		}
		
		int count = 0;
		S				= new double[image_height*image_width][2];
		T				= new double[image_height*image_width][2];
		for (int i = 0; i < image_height; i++) {
			for (int j = 0; j < image_width; j++) {
				S[count][0] = (double)i;
				S[count][1] = (double)j;
				count++;
			}
		}
		
		this.h_spatial 		= h_spatial;
		h_spatial2 			= h_spatial*h_spatial;
		knType 				= KernelType.FLAT;
		
	}
	
	public void reset(){
		S				= new double[image_height*image_width][2];
		int count = 0;
		for (int i = 0; i < image_height; i++) {
			for (int j = 0; j < image_width; j++) {
				S[count][0] = (double)i;
				S[count][1] = (double)j;
				count++;
			}
		}
	}
	
	/*
	 * non-blurring
	 */
	
	private double[] 	runOneIterAtPos(double[] curr_pos){ 
		double[] 	new_pos		= new double[2];
		double 		sum 		= 0;
		for (int x = -h_spatial; x <= h_spatial; x++) {
			if ((  curr_pos[0] + x >=0) && ( curr_pos[0] + x <= image_height-1) ){
				for (int y = -h_spatial  ; y <= h_spatial; y++) {
					if(( curr_pos[1] + y >=0) && (curr_pos[1] + y <=image_width-1)) {
						if (x * x + y * y <= h_spatial2) {
							//double value_read 		= w[x0 * image_width + y0];
							double value_read = interpolateIntensity(curr_pos[0]+x, curr_pos[1]+y);
							// possible kernel
							// / (1 + x * x) / (1 + y * y);
							// / (1 + x * x) / (1 + y * y);
							// / (1 + x * x) / (1 + y * y);
							switch(knType){
							case FLAT:
								new_pos[0] 	+= value_read * x;
								new_pos[1] 	+= value_read * y;
								sum 		+= value_read 	 ;
								
								break;
							case TRUNCATED_GAUSSIAN_1:
								new_pos[0] 	+= value_read * x * Math.exp(-(x*x+y*y));
								new_pos[1] 	+= value_read * y * Math.exp(-(x*x+y*y));
								sum 		+= value_read 	 ;
								break;
							default:
								new_pos[0] 	+= value_read * x;
								new_pos[1] 	+= value_read * y;
								sum 		+= value_read 	 ;
								break;
							}
						}
					}
				}
			}
		}
		if(sum>0){
			new_pos[0] = new_pos[0]/sum + curr_pos[0];
			new_pos[1] = new_pos[1]/sum + curr_pos[1];
		}
		else {
			new_pos[0] = curr_pos[0];
			new_pos[1] = curr_pos[1];
		}
		return new_pos;
	}

	/*
	 * wrapped columns - so that it can be applied to planisphere image where columns are circular
	 */
	
	private double[] 	runOneIterAtPos_WrapY(double[] curr_pos){ 
		double[] 	new_pos		= new double[2];
		double 		sum 		= 0;
		
		for (int x = -h_spatial; x <= h_spatial; x++) {
			if ((  curr_pos[0] + x >=0) && ( curr_pos[0] + x <= image_height-1) ){
				for (int y = -h_spatial  ; y <= h_spatial; y++) {
						if (x * x + y * y <= h_spatial2) {
							double value_read = interpolateIntensity_WrapY(curr_pos[0]+x, curr_pos[1]+y); 
							// possible kernel
							// / (1 + x * x) / (1 + y * y);
							// / (1 + x * x) / (1 + y * y);
							// / (1 + x * x) / (1 + y * y);
							switch(knType){
							case FLAT:
								new_pos[0] 	+= value_read * x;
								new_pos[1] 	+= value_read * y;
								sum 		+= value_read 	 ;
								
								break;
							case TRUNCATED_GAUSSIAN_1:
								new_pos[0] 	+= value_read * x * Math.exp(-(x*x+y*y));
								new_pos[1] 	+= value_read * y * Math.exp(-(x*x+y*y));
								sum 		+= value_read 	 ;
								break;
							default:
								new_pos[0] 	+= value_read * x;
								new_pos[1] 	+= value_read * y;
								sum 		+= value_read 	 ;
								break;
							}
						}
				}
			}
		}
		if(sum>0){
			new_pos[0] = new_pos[0]/sum + curr_pos[0];
			new_pos[1] = new_pos[1]/sum + curr_pos[1];
			// wrap it again, this time as a position
			new_pos[1] = (new_pos[1]<-0.5)?(new_pos[1]+image_width):((new_pos[1]>=image_width-0.5)?(new_pos[1]-image_width):new_pos[1]);
		}
		else {
			new_pos[0] = curr_pos[0];
			new_pos[1] = curr_pos[1];
		}
		
		return new_pos;
	}
	
	public 	double[][]	runAtPoss(double[][] poss, int max_iter, double epsilon){
		
		T				= new double[poss.length][2];
		for (int i = 0; i < poss.length; i++) {
			T[i][0] = poss[i][0];	T[i][1] = poss[i][1];
		}
		
		for (int i = 0; i < poss.length; i++) {
			int iter = 0;
			double d = Double.MAX_VALUE;
			do {
				double[] new_pos = runOneIterAtPos(T[i]);
				d = d2(new_pos, poss[i]); 
				T[i][0] = new_pos[0];
				T[i][1] = new_pos[1];
				iter++;
			}
			while (iter < max_iter&& d > epsilon*epsilon);
			
		}
		return T;
	}
	
	public double[][] 	run(int max_iter, double epsilon){
		
		T				= new double[image_height*image_width][2];
		for (int i = 0; i < T.length; i++) {
			T[i][0] = S[i][0];	T[i][1] = S[i][1];
		}
		
		for (int i = 0; i < image_height*image_width; i++) {
			int iter = 0;
			double d = Double.MAX_VALUE;
			do{
				double[] new_pos = runOneIterAtPos(T[i]);
				d = d2(new_pos, T[i]);
				T[i][0] = new_pos[0];
				T[i][1] = new_pos[1];
				iter++;
			}
			while(iter < max_iter && d > epsilon*epsilon);
		}
		
		return T;
	}

	public double[][] 	run_WrapY(int max_iter, double epsilon){
		
		T				= new double[image_height*image_width][2];
		for (int i = 0; i < T.length; i++) {
			T[i][0] = S[i][0];	T[i][1] = S[i][1];
		}
		
		for (int i = 0; i < image_height*image_width; i++) {
			int iter = 0;
			double d = Double.MAX_VALUE;
			
			do{
				
				double[] new_pos = runOneIterAtPos_WrapY(T[i]);
				
				d = d2(new_pos, T[i]);
				T[i][0] = new_pos[0];
				T[i][1] = new_pos[1];
				iter++;
			}
			while(iter < max_iter && d > epsilon*epsilon);
			
			
		}
		
		return T;
	}
	
	/*
	 * blurring
	 */
	
	public 	double[][] 	runBlurring(int max_iter, double epsilon){
		
		T				= new double[image_height*image_width][2];
		double d2_max	= Double.MIN_VALUE;
		int iter 		= 0;
		
		do{
			
			d2_max = Double.MIN_VALUE;
			
			for (int i = 0; i < S.length; i++) {
				double sum = 0;
				for (int j = 0; j < S.length; j++) {
					double dist2 = Math.pow(S[j][0]-S[i][0], 2) + Math.pow(S[j][1]-S[i][1], 2);
					if(dist2<=h_spatial2){
						
						
						T[i][0] 	+= w[j] * S[j][0];
						T[i][1] 	+= w[j] * S[j][1];
						sum 		+= w[j];
					
					}
				}
				if(sum>0){
					T[i][0] = T[i][0]/sum;
					T[i][1] = T[i][1]/sum;
				}
				else{
					T[i][0] = S[i][0];
					T[i][1] = S[i][1];
				}
				
				if(d2(S[i], T[i])>d2_max){
					d2_max = d2(S[i], T[i]);
				}
				
			}
			
			for (int i = 0; i < S.length; i++) {
				S[i][0] = T[i][0];
				S[i][1] = T[i][1];
			}
			
			iter++;
		
		}
		while(iter<=max_iter && d2_max>epsilon*epsilon); //
		
		return T;
		
	}

	public double[][] extractConvPoints(double range, int M){
		
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
					
					if(!clustered[j] && dist(T[i], T[j])<=range){
						
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
		double[][] conv_pts = null; 
		int cnt = 0;
		conv_pts = new double[nr_clusters_nigher_than_M][2];
		
		for (int i = 0; i < nr_clusters; i++) {
			if(cluster_size[i]>M){
				
				conv_pts[cnt][0] = T[i][0]; conv_pts[cnt][1] = T[i][1];
				cnt++;
			}
		}
		
		return conv_pts; // nr_pointsx2
		
	}
	
	public int[][] extractClust(int M){
		
		// create list
		ArrayList<int[]>	T_;
		T_			= new ArrayList<int[]>();
		
		// initialize the list with T values
		for (int i = 0; i < T.length; i++) {
			T_.add(new int[]{ (int)Math.round(T[i][0]) , (int)Math.round(T[i][1]) });
		}
		
		
		// final list - after cropping those with small number of elements
		ArrayList<int[]>	out_list;
		out_list	= new ArrayList<int[]>();
		
		while(!T_.isEmpty()){
			
			// take the first element
			int[] 	check_element 		= T_.get(0);
			int		check_element_count	= 1;
			T_.remove(0);
			
			int i = 0;
			while(i<T_.size()){
				int[] take_element = T_.get(i);
				
				if(equal(take_element, check_element)){
					T_.remove(i);
					check_element_count++;
				}
				else{
					i++;
				}
				
			}
			
			if(check_element_count>M){
				out_list.add(new int[]{ check_element[0], check_element[1], check_element_count });
			}
			
		}
		// out_list contains the selection of the coordinates, each corersponding to a direction
		
		int[][] out = new int[out_list.size()][3];
		for (int i = 0; i < out_list.size(); i++) {
			int[] take_dir_spherical = out_list.get(i);
			
			out[i][0] = take_dir_spherical[0];
			out[i][1] = take_dir_spherical[1];
			out[i][2] = take_dir_spherical[2];
			
		}
		
		return out;
		
	}

	
	public void printRoundedT(){
		for (int i = 0; i < T.length; i++) {
			System.out.format("%d, %d \n", (int)Math.round(T[i][0]), (int)Math.round(T[i][1]));
		}
	}
	
	public void printS(){
		for (int i = 0; i < S.length; i++) {
			System.out.format("%f, %f \n", S[i][0], S[i][1]);
		}
	}

	public void printT(){
		for (int i = 0; i < T.length; i++) {
			System.out.format("%f, %f \n", T[i][0], T[i][1]);
		}
	}
	
	public double d2(double[] a, double[] b){
		return Math.pow(a[0]-b[0], 2)+Math.pow(a[1]-b[1], 2);
	}
	
	/*
	 * private
	 */
	
	private boolean equal(int[] a, int[] b){
		return a[0]==b[0] && a[1]==b[1];
	}
	
	private double interpolateIntensity(double p1, double p2){// row,col
		// p1, p2 is a position of a pixel in an image 
		// p1 = [-0.5, H-0.5)
		// p2 = [-0.5, W-0.5) divide into 4 regions, each with different way of calculating bilinear interpolation
		// this was necessary due to changed way of calculating at the borders
		
		boolean isIn = (p1>=0 && p1<=(image_height-1) && p2>=0 && p2<=(image_width-1));
		boolean isCorner = 
				(p1<0 & p2<0) || 
				(p1>(image_height-1) & p2<0) || 
				(p1<0 & p2>image_width-1) || 
				(p1>(image_height-1) & p2>(image_width-1));
		boolean p1_border = (p1<0 || p1>(image_height-1));
		//boolean p2_border = (p2<0 || p2>(image_width-1));
		
		if(isIn){
			int[] p11 = {(int)Math.floor(p1),	(int)Math.floor(p2)}; 	double I11 = w[p11[0]*image_width+p11[1]];
			int[] p12 = {(int)Math.floor(p1), 	(int)Math.ceil(p2)}; 	double I12 = w[p12[0]*image_width+p12[1]];	
			int[] p21 = {(int)Math.ceil(p1), 	(int)Math.floor(p2)}; 	double I21 = w[p21[0]*image_width+p21[1]];
			int[] p22 = {(int)Math.ceil(p1), 	(int)Math.ceil(p2)}; 	double I22 = w[p22[0]*image_width+p22[1]];
			
			double a = ((p12[1]-p11[1])>0)?(p2-p11[1])/(p12[1]-p11[1]) : 0.5;
			double b = ((p21[0]-p11[0])>0)?(p1-p11[0])/(p21[0]-p11[0]) : 0.5;
			
			return (1-a)*(1-b)*I11 + (1-a)*b*I21 + a*(1-b)*I12 + a*b*I22;
		}
		else if(isCorner){
			return w[(int)(Math.round(p1)*image_width+Math.round(p2))];
		}
		else if(p1_border){
			int[] p11 = {(int)Math.round(p1),	(int)Math.floor(p2)}; 	double I11 = w[p11[0]*image_width+p11[1]];
			int[] p12 = {(int)Math.round(p1), 	(int)Math.ceil(p2)}; 	double I12 = w[p12[0]*image_width+p12[1]];	
			double a = ((p12[1]-p11[1])>0)?(p2-p11[1])/(p12[1]-p11[1]) : 0.5;
			return (1-a)*I11 + a*I12;
		}
		else{
			int[] p11 = {(int)Math.floor(p1),	(int)Math.round(p2)}; 	double I11 = w[p11[0]*image_width+p11[1]];
			int[] p21 = {(int)Math.ceil(p1), 	(int)Math.round(p2)}; 	double I21 = w[p21[0]*image_width+p21[1]];	
			double b = ((p21[0]-p11[0])>0)?(p1-p11[0])/(p21[0]-p11[0]) : 0.5;
			return (1-b)*I11 + b*I21;
		}
		
		
	}
	
	private double interpolateIntensity_WrapY(double p1, double p2){// row,col
		
		// p1, p2 is a real position of a pixel in an image 
		// p1 = [0.0, 	H-1.0)
		// p2 = [-0.5, 	W-0.5) 
		// width is image width - columns will be wrapped - first corresponds to the last one
		
		// find four surrounding points
		int[] p11 = {(int)Math.floor(p1),	(int)Math.floor(p2)}; 	
		int[] p12 = {(int)Math.floor(p1), 	(int)Math.ceil(p2)}; 		
		int[] p21 = {(int)Math.ceil(p1), 	(int)Math.floor(p2)}; 	
		int[] p22 = {(int)Math.ceil(p1), 	(int)Math.ceil(p2)};
		// bilinear coeffs a,b
		double a = ((p12[1]-p11[1])>0)?(p2-p11[1])/(p12[1]-p11[1]) : 0.5;
		double b = ((p21[0]-p11[0])>0)?(p1-p11[0])/(p21[0]-p11[0]) : 0.5;
		// wrap cols of surrounding pts
		p11[1] = (p11[1]<-0.5)? p11[1]+image_width : ((p11[1]>=image_width-0.5)?(p11[1]-image_width) : p11[1]);
		p12[1] = (p12[1]<-0.5)? p12[1]+image_width : ((p12[1]>=image_width-0.5)?(p12[1]-image_width) : p12[1]);
		p21[1] = (p21[1]<-0.5)? p21[1]+image_width : ((p21[1]>=image_width-0.5)?(p21[1]-image_width) : p21[1]);
		p22[1] = (p22[1]<-0.5)? p22[1]+image_width : ((p22[1]>=image_width-0.5)?(p22[1]-image_width) : p22[1]);
		
		double I11 = w[p11[0]*image_width+p11[1]];
		double I12 = w[p12[0]*image_width+p12[1]];
		double I21 = w[p21[0]*image_width+p21[1]];
		double I22 = w[p22[0]*image_width+p22[1]];
		
		return (1-a)*(1-b)*I11 + (1-a)*b*I21 + a*(1-b)*I12 + a*b*I22;
	}

	private double dist(double[] a, double[] b){
		return Math.sqrt(Math.pow((a[0]-b[0]), 2)+Math.pow((a[1]-b[1]), 2));
	}

}