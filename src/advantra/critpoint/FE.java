package advantra.critpoint;

import java.util.ArrayList;

import weka.core.Instances;
import advantra.feature.DifferentialFeatures;

import ij.ImagePlus;

public class FE {

	FEparam 		extraction_params;
	
	public FE(FEparam parameters){
		extraction_params = parameters;
	}
	
	public double[][] extractFeats(ArrayList<ImagePlus> imgs, ArrayList<Instances> locs2D){

		double[][] feat = null;
		
		// extraction
		for (int i = 0; i < imgs.size(); i++) {
			
			// diff
			DifferentialFeatures df = new DifferentialFeatures(
					imgs.get(i), 
					extraction_params.sigma_start, 
					extraction_params.sigma_end, 
					extraction_params.nr_sigmas);
			System.out.println("extracting features..."+imgs.get(i).getTitle());
			
			df.calculateFeatures();
			
			double[][] extraction_locations = new double[locs2D.get(i).numInstances()][2];
			for (int i1 = 0; i1 < extraction_locations.length; i1++) {
				extraction_locations[i1] = locs2D.get(i).instance(i1).toDoubleArray();
			}
			// extraction_locations contains [col,row]
			double[][] diff_feats = df.exportFeatures(extraction_locations, extraction_params.getDiffChoice());
			
			// some other feats
			// output is [diff_feats other_feats]
			
			// concatenate cols and then concatenate rows
			// [feat; [[diff_feats other_feats]]]
			//System.out.println("feat: "+diff_feats.length+" x "+diff_feats[0].length);
			feat = concatenateRows(feat, diff_feats);
			
		}
		
		System.out.println();
		
		return feat;
	}
	
	public String[] extractFeatLabels(boolean with_scales){
		DifferentialFeatures df = new DifferentialFeatures(
				null, 
				extraction_params.sigma_start, 
				extraction_params.sigma_end, 
				extraction_params.nr_sigmas);
		
		String[] all_labels = df.exportFeatureLabels(with_scales, extraction_params.getDiffChoice());
		// there's only differential now
		return all_labels;
	}
	
	public static String[] 	exportAllLabels(){
		// differential
		String[] all_labels = DifferentialFeatures.exportAllLabels();
		
		return all_labels;
	}

 	public double[][] extractFeats(ImagePlus img){
		
		// locations are now all the pixel locations within the image
		int curr_img_h = img.getHeight();
		int curr_img_w = img.getWidth();
		double[][] locs2D = new double[curr_img_h*curr_img_w][2];
		
		for (int loc_idx = 0; loc_idx < curr_img_h*curr_img_w; loc_idx++) {
			int[] row_col = index2location(loc_idx, curr_img_w);
			locs2D[loc_idx][0] = (double)row_col[1];//col
			locs2D[loc_idx][1] = (double)row_col[0];//row
		}
		
		// extraction
		DifferentialFeatures df = new DifferentialFeatures(
				img, 
				extraction_params.sigma_start, 
				extraction_params.sigma_end, 
				extraction_params.nr_sigmas);
		
		df.calculateFeatures();
		
		//System.out.println("locations 2d : "+locs2D.length+" , "+locs2D[0].length);
		//Utils.arrayToString(locs2D);
		double[][] diff_feats = df.exportFeatures(locs2D, extraction_params.getDiffChoice());
		
		return diff_feats; // return only diff. feats
		
	}

	private double[][] concatenateRows(double[][] in11, double[][] in21){
		
		if(in11==null){
			double[][] out = new double[in21.length][in21[0].length];
			for (int i = 0; i < in21.length; i++) {
				 for (int j = 0; j < in21[0].length; j++) {
					 out[i][j] = in21[i][j];
				 }
			}
			return out;
		}
		
		if(in11[0].length!=in21[0].length){
			System.out.println("concatenateRows: htey have to have same amount of columns");
			return null;
		}
		
		double[][] out = new double[in11.length+in21.length][in11[0].length];
		
		 for (int i = 0; i < in11.length; i++) {
			for (int j = 0; j < in11[0].length; j++) {
				out[i][j] = in11[i][j];
			}
		 }
		 
		 for (int i = 0; i < in21.length; i++) {
			 for (int j = 0; j < in21[0].length; j++) {
				 out[i+in11.length][j] = in21[i][j];
			 }
		 }
		
		return out;
		
	}
	
	public static int location2index(int row, int col, int width){
		return row*width+col;
	}
	
	public static int[] index2location(int index, int width){
		return new int[]{(index/width), (index%width)};
	}
	
}
