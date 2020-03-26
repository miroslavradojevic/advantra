package advantra.trace;

import advantra.general.ArrayHandling;
import advantra.processing.IntensityCalc;

public class BranchTrace implements Tracing {
	
												// dimensionality
	double[] 			seed_point; 			// 3
	double[]			seed_direction;			// 3
	double				seed_radius;			// 1
	
	double[]			trace_rads;				// radiuses used when tracing (linspace using radius_init, radius_step, radius_limit)
	Hypothesis[]		trace_hyps;				// hypotheses for the trace
	
	double[][]			centerlines; 			// N x3
	double[]			radiuses;				// N
	
	Hypothesis 			current_hyp_estimate;	// actual hypothesis 
	
	int 				count;					// counts points on the branch
	
	public 		BranchTrace(){ 
		
		this.seed_point 	= new double[3];
		this.seed_direction = new double[3];
		this.seed_radius = 0;
		
		trace_rads = ArrayHandling.linspace(radius_init, radius_step, radius_limit);
		trace_hyps = new Hypothesis[trace_rads.length*N_orientations];
		for (int i = 0; i < trace_rads.length*N_orientations; i++) {
			trace_hyps[i] = new Hypothesis();
		}
		
		centerlines 	= new double[N][3];
		radiuses		= new double[N];
		
		current_hyp_estimate = new Hypothesis();
		
//		new_seeds		= null;

//		sphere_img = null;
		
//		before_conv = null;
		
		count = 0;
	}

	public 		BranchTrace(final Hypothesis seed_hyp){ 
		
		seed_point 	= new double[3];
		seed_point[0]	= seed_hyp.getPositionX();
		seed_point[1]	= seed_hyp.getPositionY();
		seed_point[2]	= seed_hyp.getPositionZ();
		
		seed_direction = new double[3];
		seed_direction[0]	= seed_hyp.getOrientationX();
		seed_direction[1]	= seed_hyp.getOrientationY();
		seed_direction[2]	= seed_hyp.getOrientationZ();
		
		seed_radius 		= seed_hyp.getNeuriteRadius();
		
		trace_rads = ArrayHandling.linspace(radius_init, radius_step, radius_limit);
		trace_hyps = new Hypothesis[trace_rads.length*N_orientations];
		for (int i = 0; i < trace_rads.length*N_orientations; i++) {
			trace_hyps[i] = new Hypothesis();
		}
		
		centerlines 	= new double[N][3];
		radiuses		= new double[N];
		
		current_hyp_estimate = new Hypothesis();

//		new_seeds 		= null;
		
//		sphere_img      = null;
		
//		before_conv = null;
		
		count = 0;
	}
	
	public void 				setStart(Hypothesis		seed_hyp){ 
		
		seed_point[0]	= seed_hyp.getPositionX();
		seed_point[1]	= seed_hyp.getPositionY();
		seed_point[2]	= seed_hyp.getPositionZ();
		
		seed_direction[0]	= seed_hyp.getOrientationX();
		seed_direction[1]	= seed_hyp.getOrientationY();
		seed_direction[2]	= seed_hyp.getOrientationZ();
		
		seed_radius 		= seed_hyp.getNeuriteRadius();
		
		// add particular hypothesis to the estimation thread
		current_hyp_estimate.setHypothesis(seed_hyp);
		
		centerlines[count][0] 	= current_hyp_estimate.getPositionX();
		centerlines[count][1] 	= current_hyp_estimate.getPositionY();
		centerlines[count][2] 	= current_hyp_estimate.getPositionZ();
			
		radiuses[count] 		= current_hyp_estimate.getNeuriteRadius();
						
		count++;
		
	}
	
	public void 				storeHypothesis(Hypothesis h){
		current_hyp_estimate.setHypothesis(h);
		
		centerlines[count][0] 	= current_hyp_estimate.getPositionX();
		centerlines[count][1] 	= current_hyp_estimate.getPositionY();
		centerlines[count][2] 	= current_hyp_estimate.getPositionZ();
			
		radiuses[count] 		= current_hyp_estimate.getNeuriteRadius();
						
		count++;
		
//		System.out.format("currently:\n");
//		for (int i = 0; i < count; i++) {
//			System.out.format("%f, %f, %f ---- %f \n", 
//					centerlines[i][0], 
//					centerlines[i][1], 
//					centerlines[i][2], 
//					radiuses[i]);
//		}
	}
	
	public void 				setPriors(Hypothesis[] hyps, double prior_value){
		// the prior value will be distributed to all hypotheses from the list
		// problem.. might not sum up to 1 this way and it should...
		for (int i = 0; i < hyps.length; i++) {
			hyps[i].setPrior(prior_value);
		}
	}
	
	public void 				setLikelihoods(Hypothesis[] hyps, double likelihood_value){
		// all hypotheses are equally likely this way
		for (int i = 0; i < hyps.length; i++) {
			hyps[i].setLikelihood(likelihood_value);
		}
	}
	
	public void 				calculateLikelihoods(IntensityCalc img_calc){ 
		for (int i = 0; i < trace_hyps.length; i++) {
			trace_hyps[i].calculateLikelihood(img_calc, 1.0, 1.0);
		}
	}
	
	public double 				calculatePosteriors(){ 
		
		double sum_posteriors = 0;
		
		for (int i = 0; i < trace_hyps.length; i++) {
			trace_hyps[i].calculatePosterior();					// calculate the posterior
			sum_posteriors += trace_hyps[i].getPosterior();		// add it to the sum
		}
		// normalize posteriors (so that they sum to 1)
		if(sum_posteriors>0){
			for (int i = 0; i < trace_hyps.length; i++) {
				trace_hyps[i].scalePosterior(1/sum_posteriors);
			}
		}
		else{
			System.out.print("POSTERIOR_SUM==0");
		}
		
		return sum_posteriors;
	}
	
	public double[] 			getLastCenterpoint(){
		return (count>=1)?centerlines[count-1]:null;
	}
	
	public double 				getLastRadius(){
		return (count>=1)?radiuses[count-1]:null;
	}
	
	public double[] 			getSeedPoints(){
		return seed_point;
	}
	
	public double 				getSeedPosX(){
		return seed_point[0];
	}
	
	public double 				getSeedPosY(){
		return seed_point[1];
	}
	
	public double 				getSeedPosZ(){
		return seed_point[2];
	}
	
	public double[] 			getSeedDirection(){
		return seed_direction;
	}
	
	public double[] 			getCurrentPos(){
		return current_hyp_estimate.getPosition();
	}

	public double 				getCurrentPosX(){
		return current_hyp_estimate.getPositionX();
	}

	public double 				getCurrentPosY(){
		return current_hyp_estimate.getPositionY();
	}
	
	public double 				getCurrentPosZ(){
		return current_hyp_estimate.getPositionZ();
	}

	public double[] 			getCurrentOrient(){
		return current_hyp_estimate.getOrientation();
	}

	public double 				getCurrentOrientX(){
		return current_hyp_estimate.getOrientationX();
	}

	public double 				getCurrentOrientY(){
		return current_hyp_estimate.getOrientationY();
	}
	
	public double 				getCurrentOrientZ(){
		return current_hyp_estimate.getOrientationZ();
	}	
	
	public Hypothesis 			getCurrentTraceHypothesis(){
		return current_hyp_estimate;
	}
	
	public double[] 			getCurrentTraceHypothesisCenterpoint(){
		return current_hyp_estimate.getPosition();
	}
	
	public double 				getCurrentTraceHypothesisNeuriteRadius(){
		return current_hyp_estimate.getNeuriteRadius();
	}
	
	public double 				getCurrentTraceHypothesisRadius(){
		return current_hyp_estimate.getHypothesisRadius();
	}
	
	public Hypothesis 			getMaximumPosteriorHypothesis(){ 
		
		if(!trace_hyps[0].isPosteriorCalculated()){
			System.err.println("BranchTrace:getMaximumPosteriorHypothesis(): \n" +
					"posterior of the first hypothesis is not calculated!  Stopping...");
			return null;
		}
		
		Hypothesis out = new Hypothesis();
		
		double 	max_posterior 	= trace_hyps[0].getPosterior();
		int 	index_max 		= 0;
		
		for (int i = 1; i < trace_hyps.length; i++) {
			
			if(!trace_hyps[i].isPosteriorCalculated()){
				System.err.println("BranchTrace:getMaximumPosteriorHypothesis(): \n" +
						"posterior of hypothesis "+i+" not calculated! Return null...");
				return null;
			}
			
			if(trace_hyps[i].getPosterior()>max_posterior){
				max_posterior 	= trace_hyps[i].getPosterior();
				index_max		= i;
			}
		}
		out.setHypothesis(trace_hyps[index_max]);
		return out; 
		
	}
	
	public	static int 			getIterationLimit(){
		return N;
	}
	
}