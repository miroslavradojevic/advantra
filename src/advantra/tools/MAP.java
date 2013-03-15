package advantra.tools;

import advantra.general.ArrayHandling;
import advantra.processing.IntensityCalc;
import advantra.shapes.Sphere;
import advantra.trace.Hypothesis;
import advantra.trace.TinyBranch;

public class MAP extends Thread implements TinyBranch {

	private static 	IntensityCalc 	img_calc;
	public 	static 	Hypothesis[] 	trace_hypotheses;
	private static 	double[]		trace_rads;
	
	private static	Sphere			s;
	
	private static  double 			dr;
	private static  double 			dh;
	
	public 	static 	int				total_hyps;
	
	public  static 	Hypothesis   	in_hypothesis;
	public 	static 	Hypothesis		out_hypothesis; 
	
	private int n0, n1;
	
	public MAP(int n0, int n1){
		this.n0 = n0;
		this.n1 = n1;
	}
	
	public static void loadImage(IntensityCalc img_calc1){
		
		img_calc = img_calc1;
		
		
	}
	
	/*
	 * different modes
	 */
	
	public static void atPoint(double[] pnt){
		
		dr = dh = 1;

		in_hypothesis = null;//new Hypothesis(pnt, new double[]{1, 0, 0}, 1.0, 2);
		s = new Sphere(0, 0, 0, 1);
		// direction (1,0,0)
		double[][] orts = s.generate3DSemiSphereDirsNx3(N_orientations, 1, 0, 0);
		
		trace_rads = ArrayHandling.linspace(radius_init, radius_step, radius_limit);
		total_hyps = trace_rads.length*N_orientations;
		trace_hypotheses = new Hypothesis[total_hyps];
		
		// form the hypotheses
		int cnt = 0;
		for (int idx_ort = 0; idx_ort < N_orientations; idx_ort++) {
			for (int idx_rad = 0; idx_rad < trace_rads.length; idx_rad++) {
				
				trace_hypotheses[cnt] = new Hypothesis();//WATCH OUT! there has to be allocation here
				trace_hypotheses[cnt].setHypothesis(
						pnt[0],
						pnt[1],
						pnt[2],
						orts[idx_ort][0],
						orts[idx_ort][1],
						orts[idx_ort][2],
						trace_rads[idx_rad], 
						k);
				cnt++;

			}
		}
		
	}
	
	public static void atPoint(double[] pnt, double[] ort){
		
		dr = dh = 1;
		
		in_hypothesis = null;//new Hypothesis(pnt, new double[]{1, 0, 0}, 1.0, 2);
		s = new Sphere(0, 0, 0, 1);
		// direction defined by ort input
		double[][] orts = s.generate3DSemiSphereDirsNx3(N_orientations, ort[0], ort[1], ort[2]);
		
		trace_rads 			= ArrayHandling.linspace(radius_init, radius_step, radius_limit);
		total_hyps 			= trace_rads.length*N_orientations;
		trace_hypotheses 	= new Hypothesis[total_hyps];
		
		// form the hypotheses
		int cnt = 0;
		for (int idx_ort = 0; idx_ort < N_orientations; idx_ort++) {
			for (int idx_rad = 0; idx_rad < trace_rads.length; idx_rad++) {
				
				trace_hypotheses[cnt] = new Hypothesis();//WATCH OUT! there has to be allocation here
				trace_hypotheses[cnt].setHypothesis(
						pnt[0],
						pnt[1],
						pnt[2],
						orts[idx_ort][0],
						orts[idx_ort][1],
						orts[idx_ort][2],
						trace_rads[idx_rad], 
						k);
				cnt++;

			}
		}
		
	}
	
	public static void atHyp(Hypothesis hyp){
		
		dr = dh = 1;
		
		in_hypothesis = hyp;//new Hypothesis(pnt, new double[]{1, 0, 0}, 1.0, 2);
		s = new Sphere(
				in_hypothesis.getPositionX(), 
				in_hypothesis.getPositionY(), 
				in_hypothesis.getPositionZ(), 
				jump_ahead); //  *current_hyp_estimate.getNeuriteRadius()
		
		//3xN_orientations
		double[][] points = s.generate3DSemiSpherePts(
				N_orientations, 
				in_hypothesis.getOrientationX(), 
				in_hypothesis.getOrientationY(), 
				in_hypothesis.getOrientationZ()); 

		trace_rads 			= ArrayHandling.linspace(radius_init, radius_step, radius_limit);
		total_hyps 			= trace_rads.length*N_orientations;
		trace_hypotheses 	= new Hypothesis[total_hyps];
		
		// form the hypotheses
		int cnt = 0;
		for (int idx_ort = 0; idx_ort < N_orientations; idx_ort++) {
			for (int idx_rad = 0; idx_rad < trace_rads.length; idx_rad++) {
				
				trace_hypotheses[cnt] = new Hypothesis();//WATCH OUT! there has to be allocation here
				
				trace_hypotheses[cnt].setHypothesis(
						points[0][idx_ort], 
						points[1][idx_ort], 
						points[2][idx_ort], 
						(points[0][idx_ort]-s.getCenterX()), 
						(points[1][idx_ort]-s.getCenterY()), 
						(points[2][idx_ort]-s.getCenterZ()), 
						trace_rads[idx_rad], 
						k); 
				
				cnt++;

			}
		}
		
	}
	
	public void run(){ // will calculate posteriors
			
		for (int i = n0; i < n1; i++) {	
				
				if(in_hypothesis==null){
					
					// all have uniform weight
					trace_hypotheses[i].setPrior(1/(double)trace_hypotheses.length);
					
				}
				else{
					
					trace_hypotheses[i].calculatePrior(in_hypothesis, radius_std, direction_std);
					
				}
				
				trace_hypotheses[i].calculateLikelihood(img_calc, dr, dh);
				
				trace_hypotheses[i].calculatePosterior();
				
			}

	}
	
	/*
	 * select the one with max posterior
	 */
	
	public static Hypothesis takeMap(){
		
		double max_posterior = Double.MIN_VALUE;
		double sum_posterior = 0;
		int idx_max = 0;
		
		
		for (int i = 0; i < total_hyps; i++) {
			if(trace_hypotheses[i].isPosteriorCalculated()){
				
				sum_posterior += trace_hypotheses[i].getPosterior();
				
				if(trace_hypotheses[i].getPosterior()>max_posterior){
					max_posterior = trace_hypotheses[i].getPosterior();
					idx_max = i;
				}
				
			}
		}
		
		return (sum_posterior>0.0001)? trace_hypotheses[idx_max] : null;
		
	}
	
}
