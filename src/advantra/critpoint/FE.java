package advantra.critpoint;

import java.util.ArrayList;

import ij.ImagePlus;

public class FE {
//
//	FEparam 		extraction_params;
//	ImagePlus 		input_img; 
//	double[][]		locs;
	
	public static double[][] extract(FEparam params, ImagePlus[] imgs, ArrayList<double[][]> locs2D){
		
		int nr_locations = 0;
		for (int i = 0; i < locs2D.size(); i++) {
			nr_locations += locs2D.get(i).length;
		}
		
		
		
		double[][] feat = new double[nr_locations][];
		extraction_params = params;
		input_img = img;
		locs = new double[locs2D.length][2];
		for (int i = 0; i < locs.length; i++) {
			locs[i][0] = locs2D[i][0];
			locs[i][1] = locs2D[i][1];
		}
	}
	
	public FE(FEparam params, ImagePlus img){
		extraction_params = params;
		input_img = img;
		locs = new double[input_img.getHeight()*input_img.getWidth()][2];
		int i = 0;
		for (int row = 0; row < input_img.getHeight(); row++) {
			for (int col = 0; col < input_img.getWidth(); col++) {
				locs[i][0] = row;
				locs[i][1] = col;
				i++;
			}
		}
	}
	
	
}
