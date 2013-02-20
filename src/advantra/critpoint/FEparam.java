package advantra.critpoint;

public class FEparam {

	//multi-scale
	public double sigma_start;
	public double sigma_end;
	public int nr_sigmas;
	
	// diff feats
	static int nr_diff_feats = 14;
	public boolean[] enable_diff_feat;
	
	// other feats (to be added)
	// ...
	
	public FEparam(){
		sigma_start = 2;
		sigma_end = 3;
		nr_sigmas =2;
		enable_diff_feat = new boolean[nr_diff_feats];
	}
	
	public FEparam(double s1, double s2, int s_nr, boolean[] which_diff){
		sigma_start = s1;
		sigma_end = s2;
		nr_sigmas = s_nr;
		for (int i = 0; i < enable_diff_feat.length; i++) {
			enable_diff_feat[i] = which_diff[i];
		}
	}
	
	public int getNrFeatures(){
		
		return 1;
	}
	
	
	
}
