package advantra.processing;

import ij.ImageStack;

public class IntensityCalc {
	
	ImageStack 	img;
	// store values in an array for faster access
	int H;
	int W;
	int L;
	float[][] 	img_array;
	
	public  IntensityCalc(ImageStack img){
		this.img = img;
		H = img.getHeight();
		W = img.getWidth();
		L = img.getSize();
		this.img_array = new float[L][H*W];
		for (int i = 0; i < L; i++) {
			img_array[i] = (float [])img.getProcessor(i+1).convertToFloat().getPixels();
		}
	}
	
	public 	double interpolateAt(double p1, double p2, double p3){// p1, p2, p3 ---> row, col, lay
		double value = 0;
		boolean isIn = p1>=0 && p1<=(img.getHeight()-1) && p2>=0 && p2<=(img.getWidth()-1) && p3>=0 && p3<=(img.getSize()-1);
		if(isIn){
			int[] p11 = {(int)Math.floor(p1),	(int)Math.floor(p2)}; 	
			int[] p12 = {(int)Math.floor(p1), 	(int)Math.ceil(p2)}; 		
			int[] p21 = {(int)Math.ceil(p1), 	(int)Math.floor(p2)}; 	
			int[] p22 = {(int)Math.ceil(p1), 	(int)Math.ceil(p2)};
			
			int		l1 = (int)Math.floor(	p3);
			int		l2 = (int)Math.ceil(	p3);
			// TODO: change the way values are taken to speed it up
			double I11_1 = img.getVoxel(p11[1], p11[0], l1) ;
			double I12_1 = img.getVoxel(p12[1], p12[0], l1) ;
			double I21_1 = img.getVoxel(p21[1], p21[0], l1) ;
			double I22_1 = img.getVoxel(p22[1], p22[0], l1) ;

			double I11_2 = img.getVoxel(p11[1], p11[0], l2) ;
			double I12_2 = img.getVoxel(p12[1], p12[0], l2) ;
			double I21_2 = img.getVoxel(p21[1], p21[0], l2) ;
			double I22_2 = img.getVoxel(p22[1], p22[0], l2) ;
			
			double a = (p12[1]!=p11[1])?(p2-p11[1])/(p12[1]-p11[1]) : 0.5;	//col
			double b = (p21[0]!=p11[0])?(p1-p11[0])/(p21[0]-p11[0]) : 0.5;	//row
			double c = (l1!=l2)?(p3-l1)/(l2-l1) : 0.5;						//lay
			
			double I_1 = (1-a)*(1-b)*I11_1 + (1-a)*b*I21_1 + a*(1-b)*I12_1 + a*b*I22_1;
			double I_2 = (1-a)*(1-b)*I11_2 + (1-a)*b*I21_2 + a*(1-b)*I12_2 + a*b*I22_2;
			
			value = (1-c)*I_1+c*I_2;
			
		}
		return value;
	}

	public 	float interpolateAt_new(float p1, float p2, float p3){// p1, p2, p3 ---> row, col, lay
		float value = 0;
		boolean isIn = p1>=0 && p1<=(img.getHeight()-1) && p2>=0 && p2<=(img.getWidth()-1) && p3>=0 && p3<=(img.getSize()-1);
		if(isIn){
			int[] p11 = {(int)Math.floor(p1),	(int)Math.floor(p2)}; 	
			int[] p12 = {(int)Math.floor(p1), 	(int)Math.ceil(p2)}; 		
			int[] p21 = {(int)Math.ceil(p1), 	(int)Math.floor(p2)}; 	
			int[] p22 = {(int)Math.ceil(p1), 	(int)Math.ceil(p2)};
			
			int		l1 = (int)Math.floor(	p3);
			int		l2 = (int)Math.ceil(	p3);
			
			float I11_1 = img_array[l1][p11[0]*W+p11[1]]; //.getVoxel(p11[1], p11[0], l1) ;
			float I12_1 = img_array[l1][p12[0]*W+p12[1]]; //.getVoxel(p12[1], p12[0], l1) ;
			float I21_1 = img_array[l1][p21[0]*W+p21[1]]; //.getVoxel(p21[1], p21[0], l1) ;
			float I22_1 = img_array[l1][p22[0]*W+p22[1]]; //.getVoxel(p22[1], p22[0], l1) ;

			float I11_2 = img_array[l2][p11[0]*W+p11[1]]; //.getVoxel(p11[1], p11[0], l2) ;
			float I12_2 = img_array[l2][p12[0]*W+p12[1]]; //.getVoxel(p12[1], p12[0], l2) ;
			float I21_2 = img_array[l2][p21[0]*W+p21[1]]; //.getVoxel(p21[1], p21[0], l2) ;
			float I22_2 = img_array[l2][p22[0]*W+p22[1]]; //.getVoxel(p22[1], p22[0], l2) ;
			
			float a = (p12[1]!=p11[1])?(p2-p11[1])/(p12[1]-p11[1]) : 0.5f;	//col
			float b = (p21[0]!=p11[0])?(p1-p11[0])/(p21[0]-p11[0]) : 0.5f;	//row
			float c = (l1!=l2)?(p3-l1)/(l2-l1) : 0.5f;						//lay
			
			float I_1 = (1-a)*(1-b)*I11_1 + (1-a)*b*I21_1 + a*(1-b)*I12_1 + a*b*I22_1;
			float I_2 = (1-a)*(1-b)*I11_2 + (1-a)*b*I21_2 + a*(1-b)*I12_2 + a*b*I22_2;
			
			value = (1-c)*I_1+c*I_2;
			
		}
		return value;
	}
	
	public static double[] averageIn_Out(// TODO: check whether this one is used at all
			final double[][] extracted_cross_projection_coords,
			final double[]   extracted_cross_projection_radial_distances,
			final int        extracted_cross_projection_coords_length,
			final double     cross_projection_blob_center_x,
			final double     cross_projection_blob_center_y,
			final double     cross_projection_blob_radius,
			final int   []   extracted_cross_projection_intensities
			){
		
		double[] outputAverages = new double[2];
		
		double sumIn 	= 0;
		double sumOut 	= 0;
		int cntOut 		= 0;
		int cntIn 		= 0;

//		extracted_cross_projection_radial_distances[i] <= cross_projection_blob_radius
//		){
		// averageIn calculate
		for (int i = 0; i < extracted_cross_projection_coords_length; i++) {
			if(
					Math.sqrt(
					Math.pow(extracted_cross_projection_coords[0][i]-cross_projection_blob_center_x, 2) + 
					Math.pow(extracted_cross_projection_coords[1][i]-cross_projection_blob_center_y, 2)
					) <= cross_projection_blob_radius
							){
				// it is inside the blob (inside neurite body)
				sumIn += extracted_cross_projection_intensities[i];
				cntIn ++;
			}
			else{
				// it is outside the blob
				sumOut += extracted_cross_projection_intensities[i];
				cntOut++;
			}
			 
		}
		
		if(cntIn>0 && cntOut>0){
			
			outputAverages[0] 	= sumIn 	/ cntIn; 	// in
			outputAverages[1] 	= sumOut 	/ cntOut;	// out
			
		}
		else{
			
			outputAverages[0] 	= 0;
			outputAverages[1] 	= 0;
			System.out.println("IntensityCalculations:averageIn_Out():\n" +
					"ERROR: counted "+cntIn+" voxels in and "+cntOut+" voxels out.\n"+
					"Setting averages to zero...");
		}
		
		return outputAverages;
		
	}
	
	
}
