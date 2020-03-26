package advantra.critpoint;

import java.util.ArrayList;

public class FEparam {

	/*
	 * class that defines which features are extracted
	 * general parameters for feature extraction (scales)
	 * parameters for features that will be added in future
	 */
	
	//multi-scale
	public double 	sigma_start;
	public double 	sigma_end;
	public int 		nr_sigmas;
	
	// diff feats
	static int nr_diff_feats = 14;
	static int some_other_feats = 0;
	
	ArrayList<boolean[]> enable_feat; // each element of the list is one type of features
	
	// other feats (to be added)
	
	public FEparam(){
		sigma_start = 2;
		sigma_end = 3;
		nr_sigmas =2;
		enable_feat = new ArrayList<boolean[]>();
	}
	
	public FEparam(double s1, double s2, int s_nr, boolean[] diff_feat_enable){ // other feature types can be added here as arguments
		sigma_start 	= s1;
		sigma_end 		= s2;
		nr_sigmas 		= s_nr;
		enable_feat		= new ArrayList<boolean[]>();
		
		// diff feats
		if(diff_feat_enable!=null) {
			boolean[] diff_feat = new boolean[nr_diff_feats];
			for (int i = 0; i < nr_diff_feats; i++) {
				diff_feat[i] = diff_feat_enable[i];
			}
			enable_feat.add(diff_feat);
		}
		// add some other feats
		
	}
	
	public int getNrFeatures(){
		int nr_feat = 0;
		for (int i = 0; i < enable_feat.size(); i++) {
			if(enable_feat.get(i)!=null) {
				for (int j = 0; j < enable_feat.get(i).length; j++) {
					if(enable_feat.get(i)[j]){
						nr_feat++;
					}
				}
			}
		}
		return nr_feat;
	}
	
	public boolean[] getDiffChoice(){
		return enable_feat.get(0);
	}
	
}
