package advantra.trace;

import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.NewImage;
import ij.io.FileSaver;
import advantra.general.Transf;
import advantra.processing.IntensityCalc;
import advantra.shapes.Cylinder;
import advantra.shapes.Sphere;

public class Hypothesis {

	// hypothesis of a branch segment at particular iteration - tube model
	Cylinder 	hypothesis_cylinder;
	double		k;				// radius of the big tube - how much background around neurite
								// k>1
	
	// hypothesis features
	double 		avgIn;			// average neurite intensity 
	double 		avgOut;			// average background intensity
	boolean		intensityCalculated;
	
	double 		prior;			// prior probability [0,1]
	boolean 	priorCalculated;
	
	double 		posterior;		// posterior probability [0,1]
	boolean     posteriorCalculated;
	
	double 		likelihood;		// likelihood of 'being' a vessel
	boolean 	likelihoodCalculated;
	
	public Hypothesis(){ 		// dummy initialization
		
		k						= 2.0;
		
		hypothesis_cylinder 	= new Cylinder();
		
		// features set as well
		avgIn 					= 0;
		avgOut 					= 0;
		intensityCalculated 	= false;
		
		prior 					= 0;	
		priorCalculated 		= false;
		
		posterior 				= 0;	
		posteriorCalculated 	= false;
		
		likelihood 				= 0;	
		likelihoodCalculated	= false;
		
	}

	public Hypothesis(
			double[] 	centerpoint,
			double[] 	orientation,
			double 		neurite_radius,
			double 		k
			){ 		
		
		// define fixed height 4
		hypothesis_cylinder 	= new Cylinder(centerpoint, k*neurite_radius, 4, orientation);
		
		this.k					= (k<1)?1.0:k;
		
		// features set as well
		avgIn 					= 0;
		avgOut 					= 0;
		intensityCalculated 	= false;
		
		prior 					= 0;	
		priorCalculated 		= false;
		
		posterior 				= 0;	
		posteriorCalculated 	= false;
		
		likelihood 				= 0;	
		likelihoodCalculated	= false;
		
	}
	
	/*
	 * setHypothesis sets the geometry parameters, 
	 * not the features with this one - features are extracted 
	 * using image data as arguments
	 */
	
	public void setHypothesis( 
			// sets the Hypothesis but not it's features
			double[] 	centerpoint, 
			double[] 	orientation, 
			double 		neurite_radius, 
			double 		k
			){
		//double model_height = (2*neurite_radius<3)?3:(2*neurite_radius);
		// define fixed height 4
		this.hypothesis_cylinder.setCylinder(centerpoint, k*neurite_radius, 4, orientation);
		this.k 				= (k<1)?1.0:k;
		
		// flags cancelled 
		setFlags(false);
	}

	public void setHypothesis( 
			// sets the Hypothesis but not it's features
			double 	centerpoint_x, 
			double 	centerpoint_y,
			double 	centerpoint_z,
			double	orientation_x,
			double 	orientation_y,
			double 	orientation_z,
			double 	neurite_radius, 
			double 	k
			){
		//double model_height = (2*neurite_radius<3)?3:(2*neurite_radius);
		//this.hypothesis_cylinder.setCylinder(centerpoint, k*neurite_radius, model_height, orientation);
		
		this.hypothesis_cylinder.setCylinder(
				centerpoint_x, centerpoint_y, centerpoint_z, 
				k*neurite_radius, 4, 
				orientation_x, orientation_y, orientation_z);
		
		this.k 				= (k<1)?1.0:k;
		
		// flags cancelled 
		setFlags(false);
	}	
	
	public void setHypothesis( 
			// fixed height
			double 	centerpoint_x, 
			double 	centerpoint_y,
			double 	centerpoint_z,
			
			double	orientation_x,
			double 	orientation_y,
			double 	orientation_z,
			double 	neurite_radius, 
			double 	k,
			
			double 	height
			){
		//double model_height = (height<3)?3:(height);
		// define fixed height 4
		this.hypothesis_cylinder.setCylinder(
				centerpoint_x, centerpoint_y, centerpoint_z, 
				k*neurite_radius, 4, 
				orientation_x, orientation_y, orientation_z);
		
		this.k 				= (k<1)?1.0:k;
		
		// flags cancelled 
		setFlags(false);
	}
	
	public void setHypothesis(Hypothesis other_hypothesis){
		this.hypothesis_cylinder.setCylinder(other_hypothesis.getCylinder());
		this.k  			= other_hypothesis.getK();
		this.likelihood		= other_hypothesis.likelihood;
		this.posterior		= other_hypothesis.posterior;
		// flags cancelled 
		setFlags(false);
	}

	public void setFlags(boolean value){
		intensityCalculated 	= value;
		priorCalculated 		= value;
		posteriorCalculated 	= value;
		likelihoodCalculated	= value;
	}
	
	/*
	 * GET COMPONENTS
	 */
	
	public Cylinder getCylinder(){
		return hypothesis_cylinder;
	}
	
	public Cylinder getNeuriteCylinder(){
		Cylinder neurite_cyl = new Cylinder();
		neurite_cyl.setCylinder(this.hypothesis_cylinder);
		neurite_cyl.setR(this.getNeuriteRadius());
		return neurite_cyl;
	}
	
	public double[] getPosition(){
		return this.getCylinder().getPos();
	}

	public double getPositionY(){
		return this.getCylinder().getCenterY();
	}

	public double getPositionZ(){
		return this.getCylinder().getCenterZ();
	}
	
	public double getPositionX(){
		return this.getCylinder().getCenterX();
	}
	
	public double[] getOrientation(){
		return this.getCylinder().getV();
	}

	public double getOrientationX(){
		return this.getCylinder().getOrientX();
	}

	public double getOrientationY(){
		return this.getCylinder().getOrientY();
	}

	public double getOrientationZ(){
		return this.getCylinder().getOrientZ();
	}

	public double getHypothesisRadius(){
		return this.getCylinder().getR();
	}

	public double getNeuriteRadius(){
		
		if(k>0){
			return getCylinder().getR()/k;
		}
		else{
			System.err.println("Hypothesis parameter k was not >0. Stopping...");
			System.exit(1);
			return 99;
		}

		
	}
	
	public double getK(){
		return k;
	}
	
	public double getAvgIn(){
		return avgIn;
	}

	public double getAvgOut(){
		return avgOut;
	}

	public double getPosterior(){
		if(!isPosteriorCalculated()){
			System.err.println("Posterior is not calculated... returning "+posterior);
		}
		return posterior;
	}
	
	public double getLikelihood(){
		return likelihood;
	}

	public double getPrior(){
		return prior;
	}
	
	public boolean isIntensityCalculated(){
		return this.intensityCalculated;
	}
	
	public boolean isPriorCalculated(){
		return this.priorCalculated;
	}
	
	public boolean isLikelihoodCalculated(){
		return this.likelihoodCalculated;
	}
	
	public boolean isPosteriorCalculated(){
		return this.posteriorCalculated;
	}
	
	/*
	 * SET COMPONENTS
	 */
	
	public void setPosition(double[] position){
		this.hypothesis_cylinder.setPos(position);
		setFlags(false);
	}
	
	public void setOrientation(double[] orientation){
		this.hypothesis_cylinder.setV(orientation);
		setFlags(false);
	}
	
	public void setModelRadius(double radius){
		this.getCylinder().setR(radius);
		//this.getCylinder().setH((2*radius<3)?3:2*radius);
		setFlags(false);
	}

	public void setK(double k){
		this.k = (k<1)?1.0:k;
		setFlags(false);
	}
	
	public void invertDirection(){
		double current_vx = this.getCylinder().getOrientX();
		double current_vy = this.getCylinder().getOrientY();
		double current_vz = this.getCylinder().getOrientZ();
		//System.out.format("before %f  %f   %f \n", current_vx, current_vy, current_vz);
		this.getCylinder().setV(-current_vx, -current_vy, -current_vz);
		
	}

	/*
	 * hypothesis features: needs an image or previous state estimation as an input
	 */
	
	/*
	 *  PRIOR how much it is likely to be a  vessel (before observation)
	 */
	
	/*
	 * TODO: add the flag for prior/posterior that says whether they've been normalized
	 * this way they are just values that need to be normalized to probabilities
	 */
	public void calculatePrior(Hypothesis hyp_prev_estimate, double std_radius, double  std_orientation_angle_diff){ 
		// cyl_prev_estimate describes previous neurite chunk geometry estimate, 
		// provides with radius&direction of the previous estimate
		double angle_diff = 
				Transf.angle_between_vectors_3d_radians(
						this.hypothesis_cylinder.getV(), 
						hyp_prev_estimate.getOrientation());
		double prior_value = 
				Math.exp(-Math.pow(angle_diff, 2)								/(2*Math.pow(std_orientation_angle_diff, 2)))
				*
				Math.exp(-Math.pow((this.getNeuriteRadius()-hyp_prev_estimate.getNeuriteRadius()), 2)/(2*Math.pow(std_radius, 2)))
				;
		
		if(prior_value==Double.NaN){ // TODO: check this... apparently such comparison is not working
			System.err.println("Hypothesis:calculatePrior(): \n" +
					"prior was set as "+prior_value+" ...");
			System.exit(1);
		}
		this.setPrior(prior_value);
		
//		double[] vector1 = this.hypothesis_cylinder.getV();
//		double[] vector2 = hyp_prev_estimate.getOrientation();
//		System.out.println(
//						"obtained prior value was: "+
//						prior_value+
//						" and angle_diff was "+angle_diff+" , cos("+
//						(vector1[0]*vector2[0]+vector1[1]*vector2[1]+vector1[2]*vector2[2])+
//						")  is it >1? "+((vector1[0]*vector2[0]+vector1[1]*vector2[1]+vector1[2]*vector2[2])>1));
		
	}
	
	public void calculatePrior(double radius_prev_estimate, double[] orientation_prev_estimate, double std_radius, double  std_orientation_angle_diff){
		double angle_diff = Transf.angle_between_vectors_3d_radians(
				this.getCylinder().getV(), 
				orientation_prev_estimate);
		double prior_value = 
				Math.exp(-Math.pow(angle_diff, 2)									 /(2*Math.pow(std_orientation_angle_diff, 2)))
				*
				Math.exp(-Math.pow((this.getNeuriteRadius()-radius_prev_estimate), 2)/(2*Math.pow(std_radius, 2)))
				;
		
		setPrior(prior_value);
		
	}
	
	public void setPrior(double prior_value){
		prior = prior_value;
		if(!isPriorCalculated()){
			priorCalculated = true;		
		}
	}
	
	public void exportAsStack(ImageStack in_stack, String out_name, int radius_res, int height_res){
		
		ImagePlus im_hyp = NewImage.createByteImage("name", radius_res, radius_res, height_res, NewImage.FILL_WHITE);
		
		double row_start 	= -this.getCylinder().getR();
		double row_end 		=  this.getCylinder().getR();
		double row_range	= 2 * this.getCylinder().getR();
		
		double lay_start 	= -this.getCylinder().getH()/2;
		double lay_end 		=  this.getCylinder().getH()/2;
		double lay_range	= this.getCylinder().getH();
		
		int index_row = -1;
		for (double row = row_start; row <= row_end; row+=row_range/(radius_res-1)) {
			
			index_row++;
			
			int index_col = -1;
			for (double col = row_start; col <= row_end; col+=row_range/(radius_res-1)) {
				
				index_col++;
				
				int index_lay = -1;
				for (double lay = lay_start; lay <= lay_end; lay+=lay_range/(radius_res-1)) {
					
					index_lay++;
					
					double[] c = new double[]{this.getOrientationX(), this.getOrientationY(), this.getOrientationZ()};
					double[] a = new double[3];
					double[] b = new double[3];
					
					Transf.cartesian(c[0], c[1], c[2], a, b);
					// make it global cartesian
					double x = this.getPositionX() + row*a[0] + col*b[0] + lay*c[0];
					double y = this.getPositionY() + row*a[1] + col*b[1] + lay*c[1];
					double z = this.getPositionZ() + row*a[2] + col*b[2] + lay*c[2];
					
					if(x>=0 && x<in_stack.getHeight() && y>=0 && y<in_stack.getWidth() && z>=0 && z<in_stack.getSize()){
						
						double taken_value 	= (new IntensityCalc(in_stack)).interpolateAt_new((float)x, (float)y, (float)z);
						im_hyp.getStack().setVoxel(index_col, index_row, index_lay, taken_value);
						
					}
					
				}
			}
			
		}
		
		(new FileSaver(im_hyp)).saveAsTiffStack(out_name);
		
	}
	
	/*
	 * LIKELIHOOD
	 */
	
	// TODO: substitute this one for speed reasons with calculateLikelihood_new
	public void calculateLikelihood(ImagePlus input_img, double dr, double dh){

		double[] avg_in_out_extracted_samples = hypothesis_cylinder.extractAvgInOut(input_img, (1/this.k), dr, dh);
		
		avgIn	 		= 	avg_in_out_extracted_samples[0];
		avgOut		 	= 	avg_in_out_extracted_samples[1];
		
		double likelihood_value = 0; 
		
		if(avg_in_out_extracted_samples[2]>0.5){
			double I_contrast 	= 1.0;
			int 	s 			= 1;
			likelihood_value = Math.pow((avgIn-avgOut)/I_contrast, s); //likelihood formula
		}
		
		setLikelihood(likelihood_value);
		
	}
	
	public void calculateLikelihood_new(IntensityCalc img_calc, double dr, double dh){

		double[] avg_in_out_extracted_samples = hypothesis_cylinder.extractAvgInOut_new(img_calc, (1/k), dr, dh);
		
		avgIn	 		= 	avg_in_out_extracted_samples[0];
		avgOut		 	= 	avg_in_out_extracted_samples[1];
		
		double likelihood_value = 0; 
		
		if(avg_in_out_extracted_samples[2]>0.5){
			double I_contrast 	= 1.0;
			int 	s 			= 1;
			likelihood_value = Math.pow((avgIn-avgOut)/I_contrast, s); //likelihood formula
		}
		
		setLikelihood(likelihood_value);
		
	}
	
	public void calculateLikelihood(ImageStack img_stack){ 
		
		////// NEW CALCULATION ///////
		
		/*
		 * here...!!!
		 */
//		double 	dr = 0.5;
//		double	dh = 0.5;
		
		////// NEW CALCULATION ///////
		
		/////// THIS IS WHERE THE LIKELIHOOD CAN BE (RE)DEsigINED //////
		
		int voxel_limit_nr = Sphere.numberOfVoxInSphere(
				(int)Math.ceil(Math.sqrt(
						Math.pow(this.getCylinder().getH()/2,2)+
						Math.pow(this.getCylinder().getR(),2))));
		// allocate as usual...
		int[][] roi_coord			= new int[3][voxel_limit_nr];
		int[] 	roi_vals			= new int[voxel_limit_nr];
		double[][] roi_x_y_crosssec	= new double[2][voxel_limit_nr];
		
		// extract the values from the image stack
		int extracted_vox_nr = this.hypothesis_cylinder.extractVox(
				img_stack, roi_coord, roi_vals, roi_x_y_crosssec);
		
		//System.out.println("extracted cyl: "+extracted_vox_nr);
		
		// average in
		int 	count_in  	= 0;
		avgIn   			= 0; // average intensity in
		// average out
		int 	count_out 	= 0;
		avgOut				= 0; // average intensity out
		
		double likelihood_value = 0;
		
		if(extracted_vox_nr>voxel_limit_nr/4){  
			// put higher amount here so that those that stick out of borders are 
			// taken with 0 likelihood - maybe set the thresholds wrt volume ratio sphere/cylinder
		
		for (int i = 0; i < extracted_vox_nr; i++) {
			if(
					Math.sqrt( 
							Math.pow(roi_x_y_crosssec[0][i],2) + 
							Math.pow(roi_x_y_crosssec[1][i],2) 
							) 
							<= hypothesis_cylinder.getR()/k){
				count_in++;		avgIn 		+= roi_vals[i];
				//count_in += roi_vals[i];	avgIn 		+= roi_vals[i]*roi_vals[i];
			}
			else{
				count_out++;	avgOut 		+= roi_vals[i];
				//count_out += roi_vals[i];	avgOut 		+= roi_vals[i]*roi_vals[i];
			}
		}
		
		avgIn = 	(count_in>0)?	(avgIn/count_in)	:avgIn;
		avgOut = 	(count_out>0)?	(avgOut/count_out)	:avgOut;
		
		// calculate likelihoood
		double I_contrast 	= 1.0;
		int 	s 			= 2;
				
		// likelihood
		likelihood_value = (avgIn>avgOut)? Math.pow((avgIn-avgOut)/I_contrast, s) : 0;
		
		}
		
		/////// THIS IS WHERE THE LIKELIHOOD CAN BE (RE)DEFINED //////
		
		setLikelihood(likelihood_value);
	}
	
	public void setLikelihood(double likelihood_value){
		likelihood 		= likelihood_value;
		if(!isLikelihoodCalculated()){
			likelihoodCalculated 	= true;
		}
	}
	
	/*
	 * POSTERIOR
	 */
	
	public void calculatePosterior(){
		if(isPriorCalculated() && isLikelihoodCalculated()){
			
			posterior = this.getPrior() * this.getLikelihood();
			
//			if(posterior<0){
//				posterior = 0;
//			}
			
			posteriorCalculated = true;
			
		}
		else{
			System.err.println("Hypothesis:calculatePosterior(): \n" +
					"cannot calculate posterior without having prior & likelihood!");
			return;
		}
	}
	
	public void scalePosterior(double coeff){
		if(isPosteriorCalculated()){
			this.posterior *= coeff;
		}
		else{
			System.err.println("Hypothesis:scalePosterior(): \n" +
					"cannot scale posterior without having it calculated first!");
			System.exit(1);
		}
	}

	public static double[] 		extractLikelihoods     (Hypothesis[] hyps){
		
		double[] out = new double[hyps.length];
		
		for (int i = 0; i < hyps.length; i++) {
			out[i] = hyps[i].getLikelihood();
		}
		
		return out;
		
	}
	
	public static boolean 		isTracable(Hypothesis[] hyps, double threshold){
		
		// enable tracing if the average likelihood probability is high enough 
		// if it is low for all - then it is likely that there is no branch to trace or
		// the branch has ended
		
		if(hyps.length<=0){
			// there were no hypotheses...
			System.err.println("Hypothesis:isTracable()\n0 hypotheses to check...");
			System.exit(1);
		}
		
		double[] l = extractLikelihoods(hyps);

		double l_norm = 0;
		// normalize them - to probabilities
		for (int i = 0; i < l.length; i++) {
			l_norm += l[i];
		}
		
		if(l_norm>0){
			l_norm = l_norm/l.length;
			return (l_norm>threshold)?true:false;
		}
		else{
			return false;//probably the end
		}
		
	}
	
	public static void 			drawOverColorImage(Hypothesis[] hypotheses_to_draw, ImagePlus template_image, int valueToDraw){
		for (int i = 0; i < hypotheses_to_draw.length; i++) {
			hypotheses_to_draw[i].getCylinder().drawOverColorImage(template_image, valueToDraw);
		}
	}
	
	public static void 			drawOverColorImage(Hypothesis hypothesis_to_draw, ImagePlus template_image, int valueToDraw){
			hypothesis_to_draw.getCylinder().drawOverColorImage(template_image, valueToDraw);
	}

	public static void 			drawOverColorImage(Hypothesis hypothesis_to_draw, ImagePlus template_image, int r, int g, int b){
		hypothesis_to_draw.getCylinder().drawOverColorImage(template_image, r , g, b);
}
	
	public static double[] 		angularDiff(Hypothesis[] hyps, Hypothesis ref_hyp){
		double[] diffs = new double[hyps.length];
		for (int i = 0; i < diffs.length; i++) {
			Transf.angle_between_vectors_3d_radians(
					hyps[i].getOrientation(), ref_hyp.getOrientation());
		}
		return diffs;
	}
	
	public static double 		angularDiff(Hypothesis hyp, Hypothesis ref_hyp){
		return Transf.angle_between_vectors_3d_radians(
				hyp.getOrientation(), 
				ref_hyp.getOrientation());
	}
	
	public void 				print(){
		System.out.format( "pos[%5.2f, %5.2f, %5.2f] ort[%5.2f, %5.2f, %5.2f] r.: %5.2f | neur r.: %5.2f   k:%5.2f | prior:%f  l'hood:%f \n", 
				getPositionX(),
				getPositionY(),
				getPositionZ(),
				getOrientationX(),
				getOrientationY(),
				getOrientationZ(),
				getHypothesisRadius(),
				getNeuriteRadius(),
				getK(),
				getPrior(),
				getLikelihood()
				);
	}

}

//public void calculateLikelihood(ImageStack img_stack, int nr_rad, int nr_angl, int nr_long){
//
////DebugExport f = new DebugExport("likelihood.csv");
//
//// average in
//int 	count_in  	= 0;
//avgIn   			= 0; // average intensity in
//// average out
//int 	count_out 	= 0;
//avgOut				= 0; // average intensity out
//		
//int cnt = 0;
//double likelihood_value = 0;
//
//for (double radial_dist = 0; 
//		radial_dist <= this.hypothesis_cylinder.getR(); 
//		radial_dist+=this.hypothesis_cylinder.getR()/(nr_rad-1)) {
//	
//	
//	for (double radial_angle = 0; 
//			radial_angle < 2*Math.PI; 
//			radial_angle+=(2*Math.PI)/nr_angl) {
//		
//		
//		for (double long_dist = -this.hypothesis_cylinder.getH()/2; 
//				long_dist <= this.hypothesis_cylinder.getH()/2; 
//				long_dist+=(this.hypothesis_cylinder.getH())/(nr_long-1)) {
//			
//			// get cartesian coords
//			double cs_x = radial_dist*Math.cos(radial_angle);
//			double cs_y = radial_dist*Math.sin(radial_angle);
//			double[] x_y_z = this.hypothesis_cylinder.cyl2global(cs_x, cs_y, long_dist); 
//			
//			if(
//					x_y_z[0]>=0 && x_y_z[0]<=(img_stack.getHeight()-1) 	&&
//					x_y_z[1]>=0 && x_y_z[1]<=(img_stack.getWidth()-1) 	&&
//					x_y_z[2]>=0 && x_y_z[2]<=(img_stack.getSize()-1)
//					){
//				
//				double curr_value = (new IntensityCalc(img_stack)).interpolateAt(x_y_z[0], x_y_z[1], x_y_z[2]);
//				
//				cnt++;
//				
//				//f.writeLine(String.format("%f, %f, %f", x_y_z[0], x_y_z[1], x_y_z[2]));
//				
//				if(radial_dist<this.getNeuriteRadius()){
//					
//					count_in++;
//					avgIn += curr_value;
//				}
//				else{
//					count_out++;
//					avgOut += curr_value;
//				}
//				
//			}
//			
//		}
//		
//	}
//	
//}
//
//
////f.closeDebug();
////System.exit(1);
//
//if(cnt>0.5*nr_angl*nr_long*nr_rad){
//	
//	avgIn = 	(count_in>0)?	(avgIn/count_in)	:avgIn;
//	avgOut = 	(count_out>0)?	(avgOut/count_out)	:avgOut;
//	
//	// calculate likelihoood
//	double I_contrast 	= 1.0;
//	int 	s 			= 1;
//			
//	// likelihood
//	likelihood_value = (avgIn>avgOut)? Math.pow((avgIn-avgOut)/I_contrast, s) : 0;
//	
//}
//
//setLikelihood(likelihood_value);
//
//}
/*
 * HYPOTHESIS CREATION
 */
//double radius_step, double radius_limit,

//public static void			predictHypotheses(
//		Hypothesis[] hyps, 
//		double[] point,
//		double[] orient,
//		double radius, 
//		int N_orientations, 
//		double[] radiuses, 
//		double k){
//	// form the points on sphere - they will be placed on the semi-sphere in the direction of tangent
//	double[][] points 			= new double[3][N_orientations];
//	double[][] point_orients	= new double[3][N_orientations];
//	
//	// form them at semi-sphere with radius of that semi-sphere depending on the neurite radius linearly with scale_prediction 
//	(new Sphere(point[0], point[1], point[2], radius)).generate3DSemiSpherePts(
//					N_orientations, 
//					orient[0], 
//					orient[1], 
//					orient[2], 
//					points, 
//					point_orients);
//	
//	int radius_count	= radiuses.length;
//	
//	// set the hypotheses... once more
//	int count_hypotheses = 0;
//	System.out.format("new hypotheses construct...");
//	for (int hypo_idx = 0; hypo_idx < N_orientations; hypo_idx++) {
//		for (int radius_idx = 0; radius_idx < radius_count; radius_idx++) {
//			
//			// (re)set the Hypothesis
//			hyps[count_hypotheses].setHypothesis(
//					points[0][hypo_idx], 
//					points[1][hypo_idx], 
//					points[2][hypo_idx], 
//					point_orients[0][hypo_idx], 
//					point_orients[1][hypo_idx], 
//					point_orients[2][hypo_idx], 
//					radiuses[radius_idx], 
//					k); //k*
//			
//			count_hypotheses++;
//			
//			//System.out.print(".");
//			
//		}
//	}
//	
//	System.out.println("done.");
//	
//}