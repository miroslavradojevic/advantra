package advantra.trace;

import java.util.ArrayList;
import advantra.general.ArrayHandling;
import advantra.general.ImageConversions;
import advantra.processing.IntensityCalc;
import advantra.shapes.Cylinder;
import advantra.shapes.RegionOfInterest.RoiType;
import advantra.shapes.Sphere;
import advantra.tools.Find_Connected_Regions;
import advantra.tools.MeanShift3DSphere;
import advantra.tools.OtsuBinarisation;

import ij.ImagePlus;


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
	
	double[][] 			new_seeds; 				// br_dirs x3 (because they're 3d space coordinates)
	double[][]			new_directions;			// br_dirs x3 (3d space directions)
	int[][]				new_seeds_coords;		// br_dirs x3 (3d sphere image local coordinates)
	
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
		new_seeds_coords = null;
		
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
		new_seeds_coords = null;
		
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
		new_seeds_coords = null;
		
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
		
		storeHypothesis(seed_hyp, false);  // stores in current_hyp_estimate, centerlines, radiuses
		
	}
	
	public Hypothesis 			hyp_at(
			double[] 			point3d
			){ 
		createInitialHypothesesAndCalculateLikelihoods(point3d);
		setInitialPriors(); 
		calculatePosteriors(); 
		Hypothesis estimated_hyp 	= getMaximumPosteriorHypothesis();
		return estimated_hyp;
	}
	
	public void 				setPriors(Hypothesis[] hyps, double prior_value){
		// the prior value will be distributed to all hypotheses from the list
		// problem.. might not sum up to 1 this way and it should...
		for (int i = 0; i < hyps.length; i++) {
			hyps[i].setPrior(prior_value);
		}
	}
	
	private double[] 			subtract(double[] a3, double[] b3){
		double[] out3 = new double[3];
		out3[0] = a3[0] - b3[0];
		out3[1] = a3[1] - b3[1];
		out3[2] = a3[2] - b3[2];
		return out3;
	}
	
	private void 				normalize(double[] in3){
		double norm = Math.sqrt(Math.pow(in3[0], 2) + Math.pow(in3[1], 2) + Math.pow(in3[2], 2));
		if(norm>0){
			in3[0] /= norm;
			in3[1] /= norm;
			in3[2] /= norm;
		}
		
	}
	
	private double 				dotProd3(double[] a, double[] b){
		return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
	}
	
	public void 				setLikelihoods(Hypothesis[] hyps, double likelihood_value){
		// all hypotheses are equally likely this way
		for (int i = 0; i < hyps.length; i++) {
			hyps[i].setLikelihood(likelihood_value);
		}
	}
	
	public void 				calculateLikelihoods(){ 
		for (int i = 0; i < trace_hyps.length; i++) {
			trace_hyps[i].calculateLikelihood_new(img_calc, 1.0, 1.0);
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
	
	// TODO: check if it is used any more
	public void 				addNewBranchesToQueue(ArrayList<Hypothesis> branch_queue){ 
		
		for (int i = 0; i < new_directions.length; i++) {
			Hypothesis hyp_to_queue = new Hypothesis();
			hyp_to_queue.setHypothesis(
					centerlines[count-1][0]+check_bifurcations*current_hyp_estimate.getHypothesisRadius()*new_directions[i][0], 
					centerlines[count-1][1]+check_bifurcations*current_hyp_estimate.getHypothesisRadius()*new_directions[i][1], 
					centerlines[count-1][2]+check_bifurcations*current_hyp_estimate.getHypothesisRadius()*new_directions[i][2],
					new_directions[i][0], 
					new_directions[i][1], 
					new_directions[i][2], 
					radiuses[count-1], 
					k);
			branch_queue.add(hyp_to_queue);
			
		}
		
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
	
	public double[][] 			getNewDirections(){
		return new_directions;
	}

	public double[][] 			getNewSeeds(){
		return new_seeds;
	}
	
	public Hypothesis[] 		getNewSeedHypotheses(){
		
		if(new_directions==null){
			return null;
		}
		else{
			Hypothesis[] new_hyps = new Hypothesis[new_directions.length];
			for (int i = 0; i < new_directions.length; i++) {
				if(new_directions[i]!=null){
					new_hyps[i] = new Hypothesis(
							new_seeds[i], 
							new_directions[i], 
							current_hyp_estimate.getNeuriteRadius(), //ERROR!!! current_hyp_estimate.getHypothesisRadius(), 
							current_hyp_estimate.getK());
				}
				else{
					new_hyps[i] = null;
				}
				 
			}
			return new_hyps;
		}
		
	}
	
	public double[] 			getNewSeed(int index){
		
		if(index<0 || index>=new_seeds.length){
			return null;
		}
		else{
			return new_seeds[index];
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
	
	public void 				storeHypothesis			(Hypothesis estimated_hyp, boolean printHyp){
		
		// add particular hypothesis to the estimation thread
		current_hyp_estimate.setHypothesis(estimated_hyp);
		
		centerlines[count][0] 	= current_hyp_estimate.getPositionX();
		centerlines[count][1] 	= current_hyp_estimate.getPositionY();
		centerlines[count][2] 	= current_hyp_estimate.getPositionZ();
			
		radiuses[count] 		= current_hyp_estimate.getNeuriteRadius();
			
		if(printHyp) estimated_hyp.print();
						
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
	
	public double				relativeBifurcationSearchDistance(){ 	// relative means it's in radiuses
		return TinyBranch.check_bifurcations;
	}
	
	public double 				getK(){
		return TinyBranch.k;
	}
	
	public boolean				bifurcation_detection_MeanShift3D(ImagePlus sphere_img, Sphere sphere_from_image, boolean manual_start){ 
		
		// returns whether to stop the trace or not
		if(count<2 && !manual_start){  
			new_seeds 			= null;
			new_directions 		= null;
			new_seeds_coords 	= null;
			System.out.print("INITIAL-BLK");
			return false; // valid after 3 steps
		}
		
		int 	MS_PTS 				= 100;
		int 	MS_PTS_TH 			= 20;
		int 	MS_MAX_ITER 		= 100;
		double 	MS_EPS 				= 0.001; 
		double 	MS_NEIGHBOUR 		= 0.2; // TODO: this can be calculated wrt. circular angle & radius
		double 	MS_ANGLE_RANGE_DEG 	= 30;
		double 	MS_ANGLE_RANGE_RAD 	= (MS_ANGLE_RANGE_DEG/180)*Math.PI;
		
		MeanShift3DSphere ms3dSph = new MeanShift3DSphere(sphere_img, sphere_from_image, MS_ANGLE_RANGE_RAD, MS_PTS);
		ms3dSph.run(MS_MAX_ITER, MS_EPS);
		ms3dSph.extractClusters(MS_NEIGHBOUR, MS_PTS_TH); // result is stored in fields of the ms3dSph class
		
		double[][] 	cluster_dirs 			= ms3dSph.getClusterDirs();
		double[][] 	cluster_seed 			= ms3dSph.getClusterSeed();
		int[][] 	cluster_local_coords 	= ms3dSph.getClusterSeedLocal();
		
		if(cluster_dirs==null) 				System.out.print("MS:NO-DIRS");
		else if(cluster_dirs.length==1) 	System.out.print("MS:END");
		
		if(cluster_dirs==null || cluster_dirs.length==1){
			new_seeds 		= null;
			new_directions 	= null;
			new_seeds_coords = null;
			return true;
		}
		
		new_directions 		= new double[cluster_dirs.length][3];
		new_seeds 			= new double[cluster_dirs.length][3];
		new_seeds_coords 	= new int[cluster_dirs.length][3];

		for (int i = 0; i < cluster_dirs.length; i++) {
			// new_directions
			new_directions[i][0]= cluster_dirs[i][0];	
			new_directions[i][1]= cluster_dirs[i][1];	
			new_directions[i][2]= cluster_dirs[i][2];
			// new_seeds
			new_seeds[i][0] 	= cluster_seed[i][0];	
			new_seeds[i][1] 	= cluster_seed[i][1];	
			new_seeds[i][2] 	= cluster_seed[i][2];
			// sphere img coords
			new_seeds_coords[i][0]= cluster_local_coords[i][0];	
			new_seeds_coords[i][1]= cluster_local_coords[i][1];	
			new_seeds_coords[i][2]= cluster_local_coords[i][2];
		}
		
		if(manual_start){	
			System.out.print("MANUAL_SEED");
			return true;   // stop tracing & take all, it was just manual start
		}
		
//		System.out.println("count:"+count);
//		System.out.println("centerlines(count-1):");
//		ArrayHandling.print1DArray(centerlines[count-1]);
		
		if(cluster_dirs.length==2){
			new_seeds 		= null;
			new_directions 	= null;
			new_seeds_coords= null;
			return false;
		}
		else{ // cluster_dirs.length>2
			
			/*
			 * calculate direction towards previous
			 */
			
			double[] dir_towards_previous = new double[3];
			dir_towards_previous = subtract(centerlines[count-2], centerlines[count-1]); // centerlines[count-1] marks last centerpoint
			normalize(dir_towards_previous);
				
			/*
			 * get index to expel on the basis of previous direction
			 */
				
			int 	index_expell_direction 	= 0;
			double	dot_prod_direction_max		= dotProd3(dir_towards_previous, cluster_dirs[0]);
				
			for (int i = 1; i < cluster_dirs.length; i++) {
				double dot_prod_direction = dotProd3(dir_towards_previous, cluster_dirs[i]);
				if(dot_prod_direction>dot_prod_direction_max){
					dot_prod_direction_max = dot_prod_direction;
					index_expell_direction = i;
				}
			}
				
			new_directions[index_expell_direction] 		= null;
			new_seeds[index_expell_direction]			= null;
			new_seeds_coords[index_expell_direction]	= null;
				
				if(false){
					/*
					 * expel one that is not connected to the body
					 */
					OtsuBinarisation otsu = new OtsuBinarisation(sphere_img);
					ImagePlus sphere_img_binaized = otsu.run();
					Find_Connected_Regions conn = new Find_Connected_Regions(sphere_img_binaized, true);
					conn.run("");
					for (int i = 0; i < new_seeds_coords.length; i++) {
						if(new_seeds_coords[i]!=null && !conn.belongsToBiggestRegion(new_seeds_coords[i])){  
							new_directions[i] 		= null;
							new_seeds[i]			= null;
							new_seeds_coords[i]		= null;
						}
					}
				}
				
			int nr_new_branches = 0;
			for (int i = 0; i < new_seeds.length; i++) {
				if(new_seeds[i]!=null){
					nr_new_branches++;
				}
			}
				
			if(nr_new_branches>1){
				System.out.print("("+nr_new_branches+"new)");
			}
			else{
				new_seeds 			= null;
				new_directions 		= null;
				new_seeds_coords	= null;
				System.out.print("SUPP'SSED_BIF");
			}
				
			return (nr_new_branches>1)? true : false;  // stop in case there is more than one, otherwise 
				
		}
		
	}
	
}