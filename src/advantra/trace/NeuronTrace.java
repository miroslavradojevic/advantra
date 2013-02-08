package advantra.trace;

import java.util.ArrayList;

import advantra.processing.Moments;
import advantra.shapes.Sphere;
import advantra.swc.ExportSWC;
import advantra.trace.TinyBranchTrace.CharacteristicPoint;

import ij.ImagePlus;

public class NeuronTrace {

	ImagePlus 	traced_img;
	
	int 		max_branch_nr;
	int 		traced_branch_nr;			 
	int[]		branches_mother_idxs;	 
	
	TinyBranchTrace[] 			branches;
	ArrayList<Hypothesis>		branch_queue;
	ArrayList<Integer>			branch_queue_mother_idx;
	
	public NeuronTrace(ImagePlus im, int max_branch_nr){
		
		traced_img 	= im;
		
		this.max_branch_nr 	= max_branch_nr;
		traced_branch_nr 	= 0;
		
		branches_mother_idxs = new int[max_branch_nr];
		
		// set mother indexes to -1 
		for (int i = 0; i < branches_mother_idxs.length; i++) {
			branches_mother_idxs[i] = -1;
		}
		
		branch_queue		= new ArrayList<Hypothesis>();
		branch_queue_mother_idx	= new ArrayList<Integer>();
		
		branches = new TinyBranchTrace[max_branch_nr];
		for (int i = 0; i < max_branch_nr; i++) {
			branches[i] = new TinyBranchTrace(im);
		}
	}
	
	public void trace(double x, double y, double z){ // trace from the selected point
		
		boolean first = true;
		Hypothesis	take_hypothesis = new Hypothesis();
		Hypothesis  estimated_hyp 	= new Hypothesis();
		Hypothesis[] new_hypotheses = null;
		CharacteristicPoint p = CharacteristicPoint.BODY;
		TinyBranchTrace current_branch = new TinyBranchTrace();
		int sphere_res = 32;
		ImagePlus sphere_img = new ImagePlus();
		
		double[] start_xyz = new double[]{x, y, z};//refineStartPoint(new double[]{x, y, z}, 20); // refine 
		System.out.println("start at: "+start_xyz[0]+" , "+start_xyz[1]+" , "+start_xyz[2]);
		
		while((first || branch_queue.size()>0) && traced_branch_nr<max_branch_nr){
			
			System.out.println("looping...");//  queue size(>0)=%d, traced branches(<%d)=%d \n", branch_queue.size(), max_branch_nr, traced_branch_nr);
			
			if(first){
				System.out.println("manual press loop...");
				// detect best matching hypothesis at refined start point
				take_hypothesis = (new TinyBranchTrace(traced_img)).hyp_at(start_xyz);
				take_hypothesis.print();
				current_branch = getCurrentlyTracedBranch();
				current_branch.setStart(take_hypothesis); 
				
				//*** detection
				double scale = TinyBranch.check_bifurcations*current_branch.getCurrentTraceHypothesisRadius();
				Sphere sp = new Sphere(current_branch.getCurrentTraceHypothesisCenterpoint(), scale);
				sphere_img     = sp.extract(traced_img, sphere_res); 
				//sphere_img.show();
				current_branch.bifurcation_detection_MeanShift3D(sphere_img, sp, first); // will process new hypotheses
			
			}
			else{ 
				System.out.println("branch trace loop...");
				//take first hypothesis if not the first execution
				take_hypothesis = branch_queue.get(0);	
				branch_queue.remove(0);
				branches_mother_idxs[traced_branch_nr] = branch_queue_mother_idx.get(0);
				branch_queue_mother_idx.remove(0);
				
				current_branch = getCurrentlyTracedBranch();
				current_branch.setStart(take_hypothesis); 
				
				// trace from the point on
				new_hypotheses = null;
				int loop_count = 0;
				while(loop_count<TinyBranch.N-1){ 
					current_branch.predict();
					current_branch.calculatePriors();
					current_branch.calculateLikelihoods();
					double sum_posts = current_branch.calculatePosteriors();
					System.out.print(".");
					if(sum_posts!=0){
						
						estimated_hyp 	= current_branch.getMaximumPosteriorHypothesis();
						current_branch.storeHypothesis(estimated_hyp, false);
						
						//*** detection
						double scale = TinyBranch.check_bifurcations*current_branch.getCurrentTraceHypothesisRadius();
						Sphere sp = new Sphere(current_branch.getCurrentTraceHypothesisCenterpoint(), scale);
						sphere_img     = sp.extract(traced_img, sphere_res); 
						//sphere_img.show();
						boolean stop = current_branch.bifurcation_detection_MeanShift3D(sphere_img, sp, first); // will process new hypotheses
						if(stop) break;
					}
					else{
						// tracing itself resulted in "characteristic point"
						break;
					}
					
					if(
							p==CharacteristicPoint.TRACE_CANCELLED || 
							p==CharacteristicPoint.BRANCH ||
							p==CharacteristicPoint.END
							) break; //break while loop if the point is a stop point
					
					loop_count++;
				}
			}
			
			// done with the trace or the loop
			// store new hypotheses to the queue
			sphere_img.show();
			new_hypotheses = current_branch.getNewSeedHypotheses();
			
			if(new_hypotheses!=null){ 
				for (int b = 0; b < new_hypotheses.length; b++) {
					branch_queue.add(new_hypotheses[b]);
					if(first) 	branch_queue_mother_idx.add(traced_branch_nr);
					else 		branch_queue_mother_idx.add(-1);
				}
				System.out.println("\nadded "+new_hypotheses.length+" branches");
			}
			else{
				System.out.println("\nno added branches");
			}
			
//			System.out.format("\n------------------------------------------------------\n");
//			System.out.format("branch # %4d  |  %4d  remaining   |  limit: %4d \n", traced_branch_nr+1, branch_queue.size(), max_branch_nr);
			
			traced_branch_nr++; // count next one once this is finished
			
			if(first) first = false;
			
		}
		
		System.out.println("\n reconstructed "+traced_branch_nr+" branches");
		System.out.println("\n"+branch_queue.size()+" branches waiting in the queue");
		
	}		
	
	public void export_swc(String file_name){
		System.out.println("saving to: "+file_name);
		
		int[] last_line = new int[traced_branch_nr];
		
		ExportSWC swc_file = new ExportSWC(file_name);
		
		int line = 1;
		int prev_line;
		for (int i = 0; i < traced_branch_nr; i++) {
			boolean first = true;
			for (int j = 0; j < branches[i].getCount(); j++) {
				if(first){
					if(branches_mother_idxs[i]==(-1)){
						prev_line = -1;
					}
					else{
						prev_line = last_line[i];
					}
					first = false;
				}
				else{
					prev_line = line-1;
				}

				swc_file.writelnSWC(String.format("%5d %2d %6.2f %6.2f %6.2f %6.2f %5d", 
						line, 
						6,
						branches[i].centerlines[j][1],
						branches[i].centerlines[j][0],
						branches[i].centerlines[j][2],
						branches[i].radiuses[j],
						prev_line));
				line++;
				
			}
			
			last_line[i] = line-1;
			
		}
		
		swc_file.closeSWC();
		
	}
	
	private TinyBranchTrace getCurrentlyTracedBranch(){
		return branches[traced_branch_nr];
	}
	
//	private int				getCurrentlyTracedBranchNr(){
//		return traced_branch_nr;
//	}
	
//	private double[] 		getCurrentlyTracedBranchPos(){
//		return getCurrentlyTracedBranch().getCurrentPos();
//	}
	
//	private double 			getCurrentlyTracedBranchNeuriteRad(){
//		return getCurrentlyTracedBranch().getCurrentTraceHypothesisNeuriteRadius();
//	}	

//	private double[] AplusKtimesB(double[] A, double K, double[] B){
//		double[] c = new double[3];
//		c[0] = A[0]+K*B[0];
//		c[1] = A[1]+K*B[1];
//		c[2] = A[2]+K*B[2];
//		return c;
//	}
	
	private double[] refineStartPoint(double[] point3d, int range){
		
		// expand Sphere around the seed
		Sphere startSphere	= new Sphere(point3d[0], point3d[1], point3d[2], range);
		
		// take voxel values & locations
		int[][] roi_coord = new int[3][	Sphere.numberOfVoxInSphere(range)];	
		int[]   roi_vals  = new int[	Sphere.numberOfVoxInSphere(range)];
		
		int count_sphere_voxels = startSphere.extractVox(traced_img, roi_coord, roi_vals);
		
		double[] momts = new double[9]; // allocate space to store moments from extracted roi		
		// extract moments: momts[0], momts[1], momts[2] define CENTROID
		double sum_of_intensities = Moments.extract_moments_3D(roi_coord, roi_vals, count_sphere_voxels, momts);		
		
		double[] refined_point3d = new double[3];
		if(sum_of_intensities>0){
			refined_point3d[0] = momts[0];
			refined_point3d[1] = momts[1];
			refined_point3d[2] = momts[2];
		}
		else{
			System.out.println("Point was not refined... sum of the intensities was "+sum_of_intensities+" ...");
			refined_point3d[0] = point3d[0];
			refined_point3d[1] = point3d[1];
			refined_point3d[2] = point3d[2];
		}

		
		return refined_point3d;

	}

//	private double[][] calculateTraceDirections(double[] at_pos, double at_r, int at_resolution){
//		
//		Sphere sp = new Sphere(at_pos[0], at_pos[1], at_pos[2], at_r);
//		
//		ImagePlus output 		= sp.extract(traced_img, at_resolution);	
//		
//		double angle_range_deg 	= 20;
//		double angle_range_rad 	= (angle_range_deg/180)*Math.PI;
//		int 	N 				= 200;
//		int		N_th			= 20;
//		MeanShift3DSphere ms3dSph = new MeanShift3DSphere(output, angle_range_rad, N);
//		int max_iter 			= 200;
//		double epsilon 			= 0.00001;
//		ms3dSph.run(max_iter, epsilon);
//		ms3dSph.extractClusters(0.05, N_th); // result is stored in fields of the ms3dSph class
//		double[][] out_dirs = ms3dSph.getClusterDirs();
//		//double[][] cluster_seed = ms3dSph.getClusterSeed();
//		
//		//double[][] out_dirs = ms3dSph.extractDirections(0.05, N_th);
//		
//		(new FileSaver(output)).saveAsTiffStack(String.format("sphere_r_%f_%d_dirs.tif", at_r, out_dirs.length));
//		
//		return out_dirs;
//		
//	}
	
}