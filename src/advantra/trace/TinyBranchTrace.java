package advantra.trace;

import java.util.ArrayList;
import advantra.general.ArrayHandling;
import advantra.general.ImageConversions;
import advantra.processing.IntensityCalc;
import advantra.shapes.Cylinder;
import advantra.shapes.OrientedProjectivePlane;
import advantra.shapes.RegionOfInterest.RoiType;
import advantra.shapes.Sphere;
import advantra.tools.MeanShift;

import ij.ImagePlus;
import ij.io.FileSaver;


public class TinyBranchTrace implements TinyBranch {
	
	ImagePlus 			img_traced;				// image being traced
	IntensityCalc		img_calc;				// TODO: both ImagePlus and IntensityCalc are perhaps redundant
	
	ImagePlus			img_traced_output;		// for debug, visualizations...
	
												// dimensionality
	double[] 			seed_point; 			// 3
	double[]			seed_direction;			// 3
	double				seed_radius;			// 1
	
	double[]			trace_rads;				// radiuses used when tracing (linspace using radius_init, radius_step, radius_limit)
	Hypothesis[]		trace_hyps;				// hypotheses for the trace
	
	double[][]			centerlines; 			// N x3
	double[]			radiuses;				// N
	
	Hypothesis 			current_hyp_estimate;	// actual hypothesis 
	
	double[][] 			new_seeds; 				// br_dirs x3 (because they're 3d coordinates)
	double[][]			new_directions;			// br_dirs x3
	
	int 				count;					// counts points on the branch
	
	public static enum CharacteristicPoint {
	    END, BODY, BRANCH, TRACE_CANCELLED 
	}
	
	public 		TinyBranchTrace(){ 
		
		this.img_traced 		= null;
		this.img_calc			= new IntensityCalc();
		this.img_traced_output 	= null;
		
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
		new_directions 	= null;
		
		count = 0;
	}

	public 		TinyBranchTrace(
			final ImagePlus 		img_traced 
			){ 
		
		this.img_traced 		= img_traced;
		this.img_calc			= new IntensityCalc(img_traced.getStack());
		this.img_traced_output 	= ImageConversions.ImagePlusToRGB(img_traced);
		
		this.seed_point 	= new double[3];
		this.seed_direction = new double[3];
		
		this.seed_radius = 0;
		
		trace_rads = ArrayHandling.linspace(radius_init, radius_step, radius_limit);
		trace_hyps = new Hypothesis[trace_rads.length*N_orientations];
		for (int i = 0; i < trace_rads.length*N_orientations; i++) {
			trace_hyps[i] = new Hypothesis();
		}
		
		centerlines 	= new double[N][3];
		radiuses	= new double[N];
		
		current_hyp_estimate = new Hypothesis();

		new_directions 	= null;
		new_seeds		= null;
		
		count = 0;
	}
	
	public 		TinyBranchTrace(
			final ImagePlus 		img_traced, 
			final Hypothesis		seed_hyp
			){ 
		
		this.img_traced 		= img_traced;
		this.img_calc			= new IntensityCalc(img_traced.getStack());
		this.img_traced_output 	= ImageConversions.ImagePlusToRGB(img_traced);
		
		this.seed_point 	= new double[3];
		this.seed_point[0]	= seed_hyp.getPositionX();
		this.seed_point[1]	= seed_hyp.getPositionY();
		this.seed_point[2]	= seed_hyp.getPositionZ();
		
		this.seed_direction = new double[3];
		this.seed_direction[0]	= seed_hyp.getOrientationX();
		this.seed_direction[1]	= seed_hyp.getOrientationY();
		this.seed_direction[2]	= seed_hyp.getOrientationZ();
		
		this.seed_radius 		= seed_hyp.getNeuriteRadius();
		
		trace_rads = ArrayHandling.linspace(radius_init, radius_step, radius_limit);
		trace_hyps = new Hypothesis[trace_rads.length*N_orientations];
		for (int i = 0; i < trace_rads.length*N_orientations; i++) {
			trace_hyps[i] = new Hypothesis();
		}
		
		centerlines 	= new double[N][3];
		radiuses		= new double[N];
		
		current_hyp_estimate = new Hypothesis();

		new_seeds 		= null;
		new_directions 	= null;
		
		count = 0;
	}
	
	public void setStart(Hypothesis		seed_hyp){ 
		
		seed_point[0]	= seed_hyp.getPositionX();
		seed_point[1]	= seed_hyp.getPositionY();
		seed_point[2]	= seed_hyp.getPositionZ();
		
		seed_direction[0]	= seed_hyp.getOrientationX();
		seed_direction[1]	= seed_hyp.getOrientationY();
		seed_direction[2]	= seed_hyp.getOrientationZ();
		
		seed_radius 		= seed_hyp.getNeuriteRadius();
		
		storeHypothesis(seed_hyp, false);
		
	}
	
	public Hypothesis hyp_at(
			double[] 			point3d
			){ 
		
		//long t11 = System.currentTimeMillis();
		createInitialHypothesesAndCalculateLikelihoods(point3d);
		//long t12 = System.currentTimeMillis();
		//System.out.format("createInitialHypothesesAndCalculateLikelihoods -> %f sec.   \n", ((t12-t11)/1000f) );
		
		//t11 = System.currentTimeMillis();
		setInitialPriors(); 
		//t12 = System.currentTimeMillis();
		//System.out.format("setInitialPriors -> %f sec.   \n", ((t12-t11)/1000f) );
		
		//t11 = System.currentTimeMillis();
		calculatePosteriors(); 
		//t12 = System.currentTimeMillis();
		//System.out.format("calculatePosteriors -> %f sec.   \n", ((t12-t11)/1000f) );
		
		//t11 = System.currentTimeMillis();
		Hypothesis estimated_hyp 	= getMaximumPosteriorHypothesis();
		//t12 = System.currentTimeMillis();
		//System.out.format("getMaximumPosteriorHypothesis -> %f sec.   \n", ((t12-t11)/1000f) );
		
		//System.out.println("est.l'hood = "+estimated_hyp.getLikelihood());
		
		return estimated_hyp;
		
	}
	
	public CharacteristicPoint trace(){ 

		Hypothesis 		estimated_hyp	 	= new Hypothesis();
		
		estimated_hyp.setHypothesis(current_hyp_estimate);
		
		int loop_count = 0;
		
		CharacteristicPoint p = CharacteristicPoint.BODY;
		
		while(loop_count<N-1){  
			
			predict(); 									// make predictions & calculate priors (these two can go together as same method)
			
			calculatePriors(); 							// set priors according to previous iteration's estimation
			
			/*
			 * //TODO: speed this stage up
			 */
			calculateLikelihoods(); 					// calculate likelihoods
			
			double sum_posts = calculatePosteriors(); 	// calculate posteriors 
			
			System.out.print(".");
			
			if(sum_posts!=0){
				
				estimated_hyp 	= getMaximumPosteriorHypothesis(); // map
				storeHypothesis(estimated_hyp, false);
				
//				////////////////
//				OrientedProjectivePlane plane = new OrientedProjectivePlane(
//						current_hyp_estimate.getPositionX(), 
//						current_hyp_estimate.getPositionY(), 
//						current_hyp_estimate.getPositionZ(), 
//						current_hyp_estimate.getNeuriteRadius(), 
//						3,
//						current_hyp_estimate.getOrientationX(), 
//						current_hyp_estimate.getOrientationY(), 
//						current_hyp_estimate.getOrientationZ()
//						);
//				
//				ImagePlus plane_img = plane.extract(img_traced, 32);
//				String image_plane_name = String.format("plane_at_%f_%f_%f_count%d.tif", 
//						current_hyp_estimate.getPositionX(),
//						current_hyp_estimate.getPositionY(),
//						current_hyp_estimate.getPositionZ(),
//						0
//						);
//				(new FileSaver(plane_img)).saveAsTiff(image_plane_name);
//				
//				System.out.format("saved plane image: %s   \n", image_plane_name);
//				////////////////
				
				/*
				 * extract critical points
				 */
				
				if(count>3){ // check for bifurcations after count is ...
				
					// 3D mean-shift used to extract directions
					//extractDirections_MeanShift3D();
					
					// 2D
					extractDirections_MeanShift2D_Plane();
					
					// check whether it is not the characteristic point
					if(new_seeds==null){
						p = CharacteristicPoint.END;
						break;
					}
					else if(new_seeds.length==1){
						p = CharacteristicPoint.BODY;
					}
					else if(new_seeds.length>1){
						p = CharacteristicPoint.BRANCH;
						break;
					}
				}
				loop_count++;
			}
			else{
				System.out.format("### Stop. All posteriors were zero for current point. Branch detection was not tried.\n");
				p = CharacteristicPoint.TRACE_CANCELLED;
				break;
			}
			
		}
		
		System.out.println("finished after "+loop_count+" iterations...");
		
		return p;
		
	}
	
	public void setPriors(Hypothesis[] hyps, double prior_value){
		// the prior value will be distributed to all hypotheses from the list
		// problem.. might not sum up to 1 this way and it should...
		for (int i = 0; i < hyps.length; i++) {
			hyps[i].setPrior(prior_value);
		}
	}

	private void extractDirections_MeanShift2D_Plane(){
		
		double 		detection_distance 	= TinyBranch.check_bifurcations*current_hyp_estimate.getNeuriteRadius();
		double[] 	detection_point		= current_hyp_estimate.getPosition();	
		double[] 	detection_orient	= current_hyp_estimate.getOrientation();
		
		OrientedProjectivePlane plane = new OrientedProjectivePlane(
				detection_point, 
				detection_distance,
				2,
				detection_orient);
		
		ImagePlus output 		= plane.extract(img_traced, extract_plane_resolution);	
		
		
		MeanShift ms2dplane = new MeanShift(output, 8);
		
		int max_iter = 200;	double epsilon = 0.00001; ms2dplane.run(max_iter, epsilon);
		
		double neighbourhood = 1.0;
		double[][] conv_pts = ms2dplane.extractConvPoints(neighbourhood, (int)(0.25*Math.pow(extract_plane_resolution, 2)));
		
		if(conv_pts==null){
			(new FileSaver(output)).saveAsTiff(String.format("plane_at_%.1f_%.1f_%.1f__%d_directions.tif", 
					detection_point[0], detection_point[1], detection_point[2], 0));
			System.out.println("exporting "+String.format("plane_at_%.1f_%.1f_%.1f__%d_directions.tif", 
					detection_point[0], detection_point[1], detection_point[2], 0));
		}
		else{
			(new FileSaver(output)).saveAsTiff(String.format("plane_at_%.1f_%.1f_%.1f__%d_directions.tif", 
					detection_point[0], detection_point[1], detection_point[2], conv_pts.length));
			System.out.println("exporting "+String.format("plane_at_%.1f_%.1f_%.1f__%d_directions.tif", 
					detection_point[0], detection_point[1], detection_point[2], conv_pts.length));
			System.out.println("mean-shift detection output (conv_points):");
			ArrayHandling.print2DArray(conv_pts);
		}
		
		if(conv_pts==null){
			System.out.format("\nno directions found  \n");
			new_seeds = null; 
			new_directions = null;
		}
		else if(conv_pts.length==1){
			new_seeds 		= new double[1][3];
			new_directions 	= new double[1][3];
			new_seeds[0][0] = (conv_pts[0][0]*(2f/(extract_plane_resolution-1)) - 1) * (2*detection_distance);
			new_seeds[0][1] = (conv_pts[0][1]*(2f/(extract_plane_resolution-1)) - 1) * (2*detection_distance);
			new_seeds[0][2] = detection_distance;
			double norm = Math.sqrt(Math.pow(new_seeds[0][0], 2)+Math.pow(new_seeds[0][1], 2)+Math.pow(new_seeds[0][2], 2));
			
			new_directions[0][0] = new_seeds[0][0]/norm;
			new_directions[0][1] = new_seeds[0][1]/norm;
			new_directions[0][2] = new_seeds[0][2]/norm;
			
			new_seeds[0][0] += detection_point[0];
			new_seeds[0][1] += detection_point[1];
			new_seeds[0][2] += detection_point[2];
			
		}
		else if(conv_pts.length>=2){
			
			System.out.format("\nbranch point detected by mean-shift: %d directions \n", conv_pts.length);
			
			new_seeds 		= new double[conv_pts.length][3];
			new_directions 	= new double[conv_pts.length][3];
			
			for (int i = 0; i < conv_pts.length; i++) {
				new_seeds[i][0] = (conv_pts[i][0]*(2f/(extract_plane_resolution-1)) - 1) * (2*detection_distance);
				new_seeds[i][1] = (conv_pts[i][1]*(2f/(extract_plane_resolution-1)) - 1) * (2*detection_distance);
				new_seeds[i][2] = detection_distance;
				double norm = Math.sqrt(Math.pow(new_seeds[i][0], 2)+Math.pow(new_seeds[i][1], 2)+Math.pow(new_seeds[i][2], 2));
				
				new_directions[i][0] = new_seeds[i][0]/norm;
				new_directions[i][1] = new_seeds[i][1]/norm;
				new_directions[i][2] = new_seeds[i][2]/norm;
				
				new_seeds[i][0] += detection_point[0];
				new_seeds[i][1] += detection_point[1];
				new_seeds[i][2] += detection_point[2];
			}
				
		}
		
	}
	
	private double[] subtract(double[] a3, double[] b3){
		double[] out3 = new double[3];
		out3[0] = a3[0] - b3[0];
		out3[1] = a3[1] - b3[1];
		out3[2] = a3[2] - b3[2];
		return out3;
	}
	
	private void 	normalize(double[] in3){
		double norm = Math.sqrt(Math.pow(in3[0], 2) + Math.pow(in3[1], 2) + Math.pow(in3[2], 2));
		if(norm>0){
			in3[0] /= norm;
			in3[1] /= norm;
			in3[2] /= norm;
		}
		
	}
	
	private double 	dotProd3(double[] a, double[] b){
		return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
	}
	
	public void setLikelihoods(Hypothesis[] hyps, double likelihood_value){
		// all hypotheses are equally likely this way
		for (int i = 0; i < hyps.length; i++) {
			hyps[i].setLikelihood(likelihood_value);
		}
	}
	
	private void calculateLikelihoods(){ 
		
		for (int i = 0; i < trace_hyps.length; i++) {
			//trace_hyps[i].calculateLikelihood(img_traced.getStack());
			//trace_hyps[i].calculateLikelihood(img_traced.getStack(), likelihood_cylinder_samples);
		}
		
	}

	private double calculatePosteriors(){ 
		
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
			
			double[] point = this.getCurrentCenterpoint();
			System.err.println("\nTraceBranch:calculatePosteriors(): " +
					"sum of the posteriors was "+sum_posteriors+" at ("+point[0]+" , "+point[1]+" , "+point[2]+")");
		}
		
		return sum_posteriors;
	}
	
	public void addNewBranchesToQueue(ArrayList<Hypothesis> branch_queue){ 
		
		for (int i = 0; i < new_directions.length; i++) {
			Hypothesis hyp_to_queue = new Hypothesis();
			hyp_to_queue.setHypothesis(
					centerlines[count-1][0]+check_bifurcations*radiuses[count-1]*new_directions[i][0], 
					centerlines[count-1][1]+check_bifurcations*radiuses[count-1]*new_directions[i][1], 
					centerlines[count-1][2]+check_bifurcations*radiuses[count-1]*new_directions[i][2],
					new_directions[i][0], 
					new_directions[i][1], 
					new_directions[i][2], 
					radiuses[count-1], 
					k);
			branch_queue.add(hyp_to_queue);
			
		}
		
	}
	
	public double[] getLastCenterpoint(){
		return centerlines[count-1];
	}
	
	public double getLastRadius(){
		return radiuses[count-1];
	}
	
	public double[] getSeedPoints(){
		return seed_point;
	}
	
	public double[] getSeedDirection(){
		return seed_direction;
	}
	
	public double[] getCurrentPos(){
		return current_hyp_estimate.getPosition();
	}

	public double getCurrentPosX(){
		return current_hyp_estimate.getPositionX();
	}

	public double getCurrentPosY(){
		return current_hyp_estimate.getPositionY();
	}
	
	public double getCurrentPosZ(){
		return current_hyp_estimate.getPositionZ();
	}

	public double[] getCurrentOrient(){
		return current_hyp_estimate.getOrientation();
	}

	public double getCurrentOrientX(){
		return current_hyp_estimate.getOrientationX();
	}

	public double getCurrentOrientY(){
		return current_hyp_estimate.getOrientationY();
	}
	
	public double getCurrentOrientZ(){
		return current_hyp_estimate.getOrientationZ();
	}	
	
	public double[][] getNewOrients(){
		return new_directions;
	}
	
	public int numberOfNewOrients() {
		return new_directions.length;
	}

	public double[][] getNewSeeds(){
		return new_seeds;
	}
	
	public int numberOfNewSeeds() {
		return new_seeds.length;
	}
	
	public double[] getNewDirection(int index){
		
		if(index<0 || index>=new_directions.length){
			return null;
		}
		else{
			return new_directions[index];
		}
		
	}
	
	public double[] getNewSeed(int index){
		
		if(index<0 || index>=new_seeds.length){
			return null;
		}
		else{
			return new_seeds[index];
		}
		
	}
	
	public int getCount(){
		return count;
	}
	
	public Hypothesis 	getCurrentTraceHypothesis(){
		return current_hyp_estimate;
	}
	
	public double[] 	getCurrentCenterpoint(){
		return current_hyp_estimate.getPosition();
	}
	
	public double 		getCurrentRadius(){
		return current_hyp_estimate.getNeuriteRadius();
	}
	
	private Hypothesis 			getMaximumPosteriorHypothesis(){ 
		
		if(!trace_hyps[0].isPosteriorCalculated()){
			System.err.println("BranchTrace:getMaximumPosteriorHypothesis(): \n" +
					"posterior of the first hypothesis is not calculated!  Stopping...");
			System.exit(1);
		}
		
		Hypothesis out = new Hypothesis();
		
		double 	max_posterior 	= trace_hyps[0].getPosterior();
		int 	index_max 		= 0;
		
		for (int i = 1; i < trace_hyps.length; i++) {
			
			if(!trace_hyps[i].isPosteriorCalculated()){
				System.err.println("BranchTrace:getMaximumPosteriorHypothesis(): \n" +
						"posterior of hypothesis "+i+" not calculated! Stopping...");
				System.exit(1);
			}
			
			if(trace_hyps[i].getPosterior()>max_posterior){
				
				max_posterior 	= trace_hyps[i].getPosterior();
				index_max		= i;
			}
		}
		out.setHypothesis(trace_hyps[index_max]);
		
		return out; 
	}
	
	private void 				storeHypothesis			(Hypothesis estimated_hyp, boolean printHyp){
		
		// add particular hypothesis to the estimation thread
		current_hyp_estimate.setHypothesis(estimated_hyp);
		
		centerlines[count][0] 	= current_hyp_estimate.getPositionX();
		centerlines[count][1] 	= current_hyp_estimate.getPositionY();
		centerlines[count][2] 	= current_hyp_estimate.getPositionZ();
			
		radiuses[count] 		= current_hyp_estimate.getNeuriteRadius();
			
		if(printHyp){
				
			estimated_hyp.print();
			
		}
						
		count++;
	}
	
	public void 				drawTrace(ImagePlus template_rgb_image, RoiType which_shape, int color_r, int color_g, int color_b, double shift_for_plotting){
		
		switch(which_shape){
		case SPHERE:
			for (int cnt = 0; cnt < count; cnt++) {
				(new Sphere(
						centerlines[cnt][0]+shift_for_plotting, 
						centerlines[cnt][1]+shift_for_plotting, 
						centerlines[cnt][2]+shift_for_plotting, 
						radiuses[cnt])).drawOverColorImage(
						template_rgb_image, color_r, color_g, color_b);
			}
			break;
		case CYLINDER:
			for (int cnt = 0; cnt < count; cnt++) {
				Cylinder c = trace_hyps[cnt].getNeuriteCylinder();
				c.translatePos(shift_for_plotting);
				c.drawOverColorImage(
						template_rgb_image, color_r, color_g, color_b);
			}
			break;
		default:
			System.err.println(which_shape+" shape is not supported");
			System.exit(1);
			break;
		}
		
	}

	public	static int 			getIterationLimit(){
		return N;
	}
	
	private void 				setInitialPriors(){ 
		for (int i = 0; i < trace_hyps.length; i++) {
			trace_hyps[i].setPrior(1/(double)trace_hyps.length);
		}
	}
	
	private void 				createInitialHypotheses(
			double[] 	seed_point
			){
		
		/* 
		 * creates array of initial Hypotheses
		 * all around the same point, different directions, radiuses
		 */
		
		double[][] points 			= new double[3][N_orientations];
		double[][] point_orients	= new double[3][N_orientations];
		
		(new Sphere(seed_point[0], seed_point[1], seed_point[2], 1)).generate3DFullSpherePts(
				N_orientations, 
				points, 
				point_orients);
		
		// create initial hypotheses
		int count_hypotheses = 0;
		for (int hypo_idx = 0; hypo_idx < N_orientations; hypo_idx++) {
			for (int radius_idx = 0; radius_idx < trace_rads.length; radius_idx++) {
				
				// set the Hypothesis
				trace_hyps[count_hypotheses].setHypothesis(
						seed_point[0], seed_point[1], seed_point[2], 
						point_orients[0][hypo_idx], 
						point_orients[1][hypo_idx], 
						point_orients[2][hypo_idx], 
						trace_rads[radius_idx], 
						k,
						4
						); //2*trace_rads[radius_idx]
				count_hypotheses++;
				
			}
		}
		
	}
	
	private void 				createInitialHypothesesAndCalculateLikelihoods(
			double[] 	seed_point
			){
		
		double[][] points 			= new double[3][N_orientations];
		double[][] point_orients	= new double[3][N_orientations];
		
		(new Sphere(seed_point[0], seed_point[1], seed_point[2], 1)).generate3DFullSpherePts(
				N_orientations, 
				points, 
				point_orients);
		
		// create initial hypotheses
		int count_hypotheses = 0;
		for (int hypo_idx = 0; hypo_idx < N_orientations; hypo_idx++) {
			for (int radius_idx = 0; radius_idx < trace_rads.length; radius_idx++) {
				
				// set the Hypothesis
				trace_hyps[count_hypotheses].setHypothesis(
						seed_point[0], seed_point[1], seed_point[2], 
						point_orients[0][hypo_idx], 
						point_orients[1][hypo_idx], 
						point_orients[2][hypo_idx], 
						trace_rads[radius_idx], 
						k,
						4
						); 

				trace_hyps[count_hypotheses].calculateLikelihood_new(img_calc, 1.0, 1.0);
				
//				System.out.println(
//						"l'hood["+count_hypotheses+"] = "+trace_hyps[count_hypotheses].getLikelihood());
				count_hypotheses++;
				
			}
		}
		
	}

	private void				predict(){
		
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
	
//	/*
//	 * make this one to speed up the estimations
//	 */
//	//TODO finish this to speed up 
//	private void 				likelihoods(int nr_samples){
//		
//		
//		for (int i = 0; i < nr_samples; i++) {
//			Random generator = new Random();
//			double r = generator.nextDouble()*trace_rads[trace_rads.length-1];
//			double a = generator.nextDouble()*2*Math.PI;
//			double h = (generator.nextDouble()-0.5)*4;
//			double cs_x = r*Math.cos(a);
//			double cs_y = r*Math.sin(a);
//		}
//		
//	}
	
	private void 				calculatePriors(){ 
		// it will set priors in hyps w.r.t. ref_hyp
		for (int i = 0; i < trace_hyps.length; i++) {
			
			trace_hyps[i].calculatePrior(current_hyp_estimate, radius_std, direction_std);
			
		}
	}

	public double				relativeJumpDist(){ 					// relative means it's in radiuses
		return TinyBranch.jump_ahead;
	}
	
	public double				relativeBifurcationSearchDistance(){ 	// relative means it's in radiuses
		return TinyBranch.check_bifurcations;
	}
	
	public double 				getK(){
		return TinyBranch.k;
	}
	
}

//private void extractDirections_MeanShift3D(){
//
//Sphere sp = new Sphere(
//		current_hyp_estimate.getPositionX(), 
//		current_hyp_estimate.getPositionY(), 
//		current_hyp_estimate.getPositionZ(), 
//		current_hyp_estimate.getNeuriteRadius()*check_bifurcations);
//
//ImagePlus output 		= sp.extract(img_traced, extract_sphere_resolution);	
//
//double angle_range_deg 	= 30;
//double angle_range_rad 	= (angle_range_deg/180)*Math.PI;
//MeanShift3DSphere ms3dSph = new MeanShift3DSphere(output, angle_range_rad, number_of_convergence_points);
//int max_iter = 200;	double epsilon = 0.00001; ms3dSph.run(max_iter, epsilon);
//double neighbourhood = 0.01;
//double[][] out_dirs = ms3dSph.extractDirections(neighbourhood, threshold_convergence_points);
//	
//if(out_dirs==null){
//	System.out.format("\n no directions found by mean-shift 3d  \n");
//	sp.drawOverColorImage(img_traced_output, 0, 0, 255); // blue
//	new_seeds = null; 
//}
//else if(out_dirs.length==1){
//	System.out.format("\nend point detected by mean-shift 3d\n");
//	
//	/////
//	(new FileSaver(output)).saveAsTiffStack(String.format("end_pt_at_%d_start_%f_%f_%f__dir_%f_%f_%f.tif", 
//			count, seed_point[0], seed_point[1], seed_point[2], seed_direction[0], seed_direction[1], seed_direction[2]));
//	sp.drawOverColorImage(img_traced_output, 0, 255, 0); // green
//	/////
//	
//	new_seeds = null;
//}
//else if(out_dirs.length>2){
//	
//	System.out.format("\n$$$ branch point detected by mean-shift 3d: %d directions \n", out_dirs.length);
//	
//	/////
//	(new FileSaver(output)).saveAsTiffStack(String.format("bifurcation_%d_thresholded_directs_at_%d_start_%f_%f_%f__dir_%f_%f_%f.tif", 
//			out_dirs.length, count, seed_point[0], seed_point[1], seed_point[2], seed_direction[0], seed_direction[1], seed_direction[2]));
//	sp.drawOverColorImage(img_traced_output, 255, 0, 0); // red
//	/////
//	
//	/*
//	 * compare 'out_dirs' directions with the direction towards previous
//	 * exclude the one that is the closest to previous direction
//	 * and store the final group of directions to current_planisph_directions
//	 * important that storeHypothesis() happened before so that count is set to proper value
//	 */
//	
//	// expel one from new_seeds - take out the one closest to the direction towards previous trace point
//		
//		double[] dir_towards_previous = new double[3];
//		dir_towards_previous = subtract(centerlines[count-2], centerlines[count-1]);
//		normalize(dir_towards_previous);
//		
//		int 	index_expell_direction 	= 0;
//		double	dot_prod_direction_max		= dotProd3(dir_towards_previous, out_dirs[0]);
//		
//		for (int i = 1; i < out_dirs.length; i++) {
//			double dot_prod_direction = dotProd3(dir_towards_previous, out_dirs[i]);
//			if(dot_prod_direction>dot_prod_direction_max){
//				dot_prod_direction_max = dot_prod_direction;
//				index_expell_direction = i;
//			}
//		}
//		new_seeds = new double[out_dirs.length-1][3];
//		int cnt = 0;
//		for (int i = 0; i < out_dirs.length; i++) {
//			if(i!=index_expell_direction){
//				new_seeds[cnt][0] = out_dirs[i][0];
//				new_seeds[cnt][1] = out_dirs[i][1];
//				new_seeds[cnt][2] = out_dirs[i][2];
//				cnt++;
//			}
//			
//		}
//		
//}
//
//}
