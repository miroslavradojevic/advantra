package advantra.critpoint;

import java.util.ArrayList;

import advantra.feature.DifferentialFeatures;

import ij.ImagePlus;

public class FE {

	FEparam 		extraction_params;
	
	public FE(FEparam parameters){
		extraction_params = parameters;
	}
	
	public double[][] extractFeats(ArrayList<ImagePlus> imgs, ArrayList<double[][]> locs2D){

		double[][] feat = null;
		
		// extraction
		for (int i = 0; i < imgs.size(); i++) {
			
			// diff
			DifferentialFeatures df = new DifferentialFeatures(
					imgs.get(i), 
					extraction_params.sigma_start, 
					extraction_params.sigma_end, 
					extraction_params.nr_sigmas);
			double[][] diff_feats = df.exportFeatures(locs2D.get(i), extraction_params.getDiffChoice());
			
			// some other feats
			// output is [diff_feats other_feats]
			
			// concatenate cols and then concatenate rows
			// [feat; [[diff_feats other_feats]]]
			System.out.println("feat: "+diff_feats.length+" x "+diff_feats[0].length);
			feat = concatenateRows(feat, diff_feats);
			
		}
		
		return feat;
	}
	
	public String[] extractFeatLabels(){
		String[] all_labels = DifferentialFeatures.exportFeatureLabels(extraction_params.nr_sigmas, extraction_params.getDiffChoice());// there's only differential now
		return all_labels;
	}

	public double[][] extractFeats(ArrayList<ImagePlus> imgs){

		// locations are now all the pixel locations within the image
		ArrayList<double[][]> locs2D = new ArrayList<double[][]>(imgs.size());
		for (int i = 0; i < imgs.size(); i++) {
			int curr_img_h = imgs.get(i).getHeight();
			int curr_img_w = imgs.get(i).getWidth();
			double[][] locs = new double[curr_img_h*curr_img_w][2];
			for (int loc_idx = 0; loc_idx < curr_img_h*curr_img_w; loc_idx++) {
				int[] curr_loc = index2location(loc_idx, curr_img_w);
				locs[loc_idx][0] = (double)curr_loc[0];
				locs[loc_idx][1] = (double)curr_loc[1];
			}
			locs2D.add(locs);
		}
		
		double[][] feat = null;
		
		// extraction
		for (int i = 0; i < imgs.size(); i++) {
			// diff
			DifferentialFeatures df = new DifferentialFeatures(
					imgs.get(i), 
					extraction_params.sigma_start, 
					extraction_params.sigma_end, 
					extraction_params.nr_sigmas);
			
			double[][] diff_feats = df.exportFeatures(locs2D.get(i), extraction_params.getDiffChoice());
			
			// some other feats
			// output is [diff_feats other_feats]
			
			// concatenate cols and then concatenate rows
			// [feat; [[diff_feats other_feats]]]
				
			feat = concatenateRows(feat, diff_feats);
				
		}
		
		return feat;
		
	}
	
	public double[][] extractFeats(ImagePlus img){

		// locations are now all the pixel locations within the image
		int curr_img_h = img.getHeight();
		int curr_img_w = img.getWidth();
		double[][] locs2D = new double[curr_img_h*curr_img_w][2];
		for (int loc_idx = 0; loc_idx < curr_img_h*curr_img_w; loc_idx++) {
			int[] curr_loc = index2location(loc_idx, curr_img_w);
			locs2D[loc_idx][0] = (double)curr_loc[0];
			locs2D[loc_idx][1] = (double)curr_loc[1];
		}
		
		// extraction
		DifferentialFeatures df = new DifferentialFeatures(
				img, 
				extraction_params.sigma_start, 
				extraction_params.sigma_end, 
				extraction_params.nr_sigmas);
		
		/// DEBUG
		
//		System.out.println("rows:");
//		for (int i = 0; i < locs2D.length; i++) {
//			System.out.print(locs2D[i][0]+" , ");
//		}
//		System.out.println("cols:");
//		for (int i = 0; i < locs2D.length; i++) {
//			System.out.print(locs2D[i][1]+" , ");
//		}
		
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
