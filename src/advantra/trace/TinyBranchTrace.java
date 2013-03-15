package advantra.trace;

import java.util.Vector;

import advantra.general.ArrayHandling;
import advantra.processing.IntensityCalc;
import advantra.shapes.Sphere;
import advantra.tools.MeanShift3DSphere;
import ij.ImagePlus;


public class TinyBranchTrace implements TinyBranch {
	
												// dimensionality
	double[] 			seed_point; 			// 3
	double[]			seed_direction;			// 3
	double				seed_radius;			// 1
	
	double[]			trace_rads;				// radiuses used when tracing (linspace using radius_init, radius_step, radius_limit)
	Hypothesis[]		trace_hyps;				// hypotheses for the trace
	
	double[][]			centerlines; 			// N x3
	double[]			radiuses;				// N
	
	Hypothesis 			current_hyp_estimate;	// actual hypothesis 
	
	public double[][] 	new_seeds; 				// br_dirs x3 (because they're 3d space coordinates)
//	double[][]			new_seeds_test;
	
	ImagePlus 			sphere_img; 
	//double[][] 			after_conv;
	double[][]			before_conv;
	//double[][]			new_directions;			// br_dirs x3 (3d space directions)
	//int[][]				new_seeds_coords;		// br_dirs x3 (3d sphere image local coordinates)
	
	int 				count;					// counts points on the branch
	
	public 		TinyBranchTrace(){ 
		
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
		
		new_seeds		= null;
//		new_seeds_test 	= null;

		sphere_img = null;
		
		//after_conv = null;
		before_conv = null;
		
		count = 0;
	}

	public 		TinyBranchTrace(final Hypothesis seed_hyp){ 
		
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

		new_seeds 		= null;
//		new_seeds_test  = null;
		
		sphere_img      = null;
		
		//after_conv = null;
		before_conv = null;
		
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
	
	public void storeHypothesis(Hypothesis h){
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
	
	public double[][] 			getNewSeeds(){
		return new_seeds;
	}
	
	public Vector<Hypothesis> 	getNewSeedHypotheses(IntensityCalc img_calc, boolean manual_start){
		
		Vector<Hypothesis> new_hyps = new Vector<Hypothesis>();
		
		if(new_seeds==null){
			
			return new_hyps;
			
		}
		else{
			
			boolean expelled = false; 
			
			for (int i = 0; i < new_seeds.length; i++) {
				
				//if(new_seeds[i]!=null){
					
					if(!manual_start) {
						
						/* check if it belongs to the current branch already
						 * chk checkings at max
						 */

						int chk = 0;
						
						for (int k = count-1; k >=0; k--) {
							
							double dst2 = 
									Math.pow((new_seeds[i][0]-centerlines[k][0]), 2)+
									Math.pow((new_seeds[i][1]-centerlines[k][1]), 2)+
									Math.pow((new_seeds[i][2]-centerlines[k][2]), 2);
							
							chk++;
							
							if(chk>3) break;
							
							if(dst2<1.3*radiuses[k]*1.3*radiuses[k]) {
								
								new_seeds[i] 	= null;
								expelled = true;
								break;
								
							}
							
						}
						
					}
					
					if(new_seeds[i]!=null) { // actually makes sense in case it is was a manual start
					
						double[] S_vec = new double[3]; // seed orientation
						S_vec[0] = new_seeds[i][0]-current_hyp_estimate.getPositionX();
						S_vec[1] = new_seeds[i][1]-current_hyp_estimate.getPositionY();
						S_vec[2] = new_seeds[i][2]-current_hyp_estimate.getPositionZ();
						
						// not necessary to normalize
						double S_vec_norm = Math.sqrt(S_vec[0]*S_vec[0]+S_vec[1]*S_vec[1]+S_vec[2]*S_vec[2]);
						S_vec[0] /= S_vec_norm;
						S_vec[1] /= S_vec_norm;
						S_vec[2] /= S_vec_norm;
						
						double r = current_hyp_estimate.getNeuriteRadius();
						
						Hypothesis hyp_to_add = new Hypothesis(new_seeds[i], S_vec, r, k);
						new_hyps.add(hyp_to_add);

					}

				//}
				 
			}
			
			if(!expelled && !manual_start){
				new_hyps.clear();
				return new_hyps;
			}
			else{
				return new_hyps;
			}
			
			 
		}
		
	}
	
	public int 					getCount(){
		return count;
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
	
	public void					predict(){
		
		// form the points on sphere... they will be placed on the semi-sphere in the direction of tangent 
		//double[][] points 			= new double[3][N_orientations];
		//double[][] point_orients	= new double[3][N_orientations];
		
		// form them at semi-sphere with radius of that semi-sphere depending on the neurite radius linearly with scale_prediction 
		Sphere prediction_sphere = new Sphere(
				current_hyp_estimate.getPositionX(), 
				current_hyp_estimate.getPositionY(), 
				current_hyp_estimate.getPositionZ(), 
				jump_ahead);  //  *current_hyp_estimate.getNeuriteRadius()
		
		double[][] points = prediction_sphere.generate3DSemiSpherePts(
						N_orientations, 
						current_hyp_estimate.getOrientationX(), 
						current_hyp_estimate.getOrientationY(), 
						current_hyp_estimate.getOrientationZ()); 
		
		// set the hypotheses... 
		int count_hypotheses = 0;
		
		for (int hypo_idx = 0; hypo_idx < N_orientations; hypo_idx++) {
			for (int radius_idx = 0; radius_idx < trace_rads.length; radius_idx++) {
				
				// (re)set the Hypothesis
				trace_hyps[count_hypotheses].setHypothesis(
						points[0][hypo_idx], 
						points[1][hypo_idx], 
						points[2][hypo_idx], 
						(points[0][hypo_idx]-prediction_sphere.getCenterX())/prediction_sphere.getR(), 
						(points[1][hypo_idx]-prediction_sphere.getCenterY())/prediction_sphere.getR(), 
						(points[2][hypo_idx]-prediction_sphere.getCenterZ())/prediction_sphere.getR(), 
						trace_rads[radius_idx], 
						k); 
				
				count_hypotheses++;
				
			}
		}
		
	}
	
	public void 				calculatePriors(){ 
		// it will set priors in hyps w.r.t. ref_hyp
		for (int i = 0; i < trace_hyps.length; i++) {
			trace_hyps[i].calculatePrior(current_hyp_estimate, radius_std, direction_std);
		}
	}

	public double				relativeJumpDist(){ 					// relative means it's in radiuses
		return TinyBranch.jump_ahead;
	}
	
	public double 				getK(){
		return TinyBranch.k;
	}
	
	/*
	 * MEAN-SHIFT to calculate 
	 */
	
	public void					calculateNewSeeds(IntensityCalc img_calc){ 
		
		int 	MS_PTS 				= 8*20;
		int 	CPU_NR				= 8;
		int 	MS_PTS_TH 			= 10; // depends on sensitivity
		int 	MS_MAX_ITER 		= 100;
		double 	MS_EPS 				= 0.001; 
		double 	MS_NEIGHBOUR 		= 2;//deg
		double  MS_NEIGHBOUR_RAD	= (MS_NEIGHBOUR/180)*Math.PI;
		double 	MS_ANGLE_RANGE_DEG 	= 30;
		double 	MS_ANGLE_RANGE_RAD 	= (MS_ANGLE_RANGE_DEG/180)*Math.PI;
		
		/*
		 *  around getCurrentTraceHypothesisCenterpoint() 
		 */
		double rr = current_hyp_estimate.getNeuriteRadius();
		Sphere sphere_from_image 		= new Sphere(getCurrentTraceHypothesisCenterpoint(), ((rr<=2)? (3*rr) : (2*rr)));
		sphere_img     					= sphere_from_image.extract(img_calc, TinyBranch.extract_sphere_resolution); 
		
		MeanShift3DSphere.load(sphere_img, sphere_from_image, MS_ANGLE_RANGE_RAD, MS_PTS, MS_MAX_ITER, MS_EPS);
		MeanShift3DSphere ms_jobs[] = new MeanShift3DSphere[CPU_NR];
		for (int i = 0; i < ms_jobs.length; i++) {
			ms_jobs[i] = new MeanShift3DSphere(i * MS_PTS / CPU_NR,  (i+1) * MS_PTS / CPU_NR);
			ms_jobs[i].start();
		}		
		for (int i = 0; i < ms_jobs.length; i++) {
			try {
				ms_jobs[i].join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}		

		
		MeanShift3DSphere.extractClusters(MS_NEIGHBOUR_RAD, MS_PTS_TH); 
		
		new_seeds 			= MeanShift3DSphere.cluster_seed;
		
//		System.out.println("new_seeds from ms:");
//		ArrayHandling.print2DArray(new_seeds);
		
	}
	
	public ImagePlus 			getSphereImage(){
		return sphere_img;
	}
	
	public double[][] 			getBeforeConv(){
		return before_conv;
	}

}

//public double[][]			getNewSeedsTest(){
//return new_seeds_test;
//}

//public double[][] 			getNewSeedsCopy(){
//double[][] out = new double[new_seeds.length][new_seeds[0].length];
//for (int i = 0; i < new_seeds.length; i++) {
//	for (int j = 0; j < new_seeds[0].length; j++) {
//		out[i][j] = new_seeds[i][j];
//	}
//}
//return out;
//}
//private void 				setInitialPriors(){ 
//for (int i = 0; i < trace_hyps.length; i++) {
//	trace_hyps[i].setPrior(1/(double)trace_hyps.length);
//}
//}
//private void 				createInitialHypothesesAndCalculateLikelihoods(double[] seed_point,	IntensityCalc img_calc){
//
//double[][] points 			= new double[3][N_orientations];
//double[][] point_orients	= new double[3][N_orientations];
//
////(new Sphere(seed_point[0], seed_point[1], seed_point[2], 1)).generate3DFullSpherePts(
////		N_orientations, 
////		points, 
////		point_orients);
//
//(new Sphere(seed_point[0], seed_point[1], seed_point[2], 1)).generate3DSemiSpherePts(
//		N_orientations, 
//		1, 0, 0, 
//		points, 
//		point_orients);
//// create initial hypotheses
//int count_hypotheses = 0;
//for (int hypo_idx = 0; hypo_idx < N_orientations; hypo_idx++) {
//	for (int radius_idx = 0; radius_idx < trace_rads.length; radius_idx++) {
//		
//		// set the Hypothesis
//		trace_hyps[count_hypotheses].setHypothesis(
//				seed_point[0], seed_point[1], seed_point[2], 
//				point_orients[0][hypo_idx], 
//				point_orients[1][hypo_idx], 
//				point_orients[2][hypo_idx], 
//				trace_rads[radius_idx], 
//				k,
//				2
//				); 
//
//		trace_hyps[count_hypotheses].calculateLikelihood(img_calc, 1, 1);
//		
////		System.out.println(
////				"l'hood["+count_hypotheses+"] = "+trace_hyps[count_hypotheses].getLikelihood());
//		count_hypotheses++;
//		
//	}
//}
//
//}
//public Hypothesis 			hyp_at(double[] 		point3d, IntensityCalc 		img_calc){ 
//	createInitialHypothesesAndCalculateLikelihoods(point3d, img_calc);
//	setInitialPriors(); 
//	calculatePosteriors(); 
//	Hypothesis estimated_hyp 	= getMaximumPosteriorHypothesis();
//	return estimated_hyp;
//}