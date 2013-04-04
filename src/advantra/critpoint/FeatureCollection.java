package advantra.critpoint;

import java.util.ArrayList;

public class FeatureCollection {

	ArrayList<String>					sources 	= 	new ArrayList<String>();
	ArrayList<Double>				 	scales 		= 	new ArrayList<Double>();
	ArrayList<String[]>				 	labels 		= 	new ArrayList<String[]>();
	ArrayList<ArrayList<double[][]>> 	feats 		= 	new ArrayList<ArrayList<double[][]>>();
	
	int pool_size = 4;
	
	/*
	 * 1. differential invariants
	 * 2. gabor features
	 */
	
	
	public FeatureCollection(int nr_scales, int nr_sources){
		
		
		
	}
	
	public int getNrFeats(){
		
		int nr = 0;
		
		for (int i = 0; i < scales.size(); i++) {
//			nr += sca;
		}
		
		return 1;
		
	}
	
	
	
}
