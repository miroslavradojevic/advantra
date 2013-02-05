package advantra.trace;

import java.util.ArrayList;

import advantra.processing.Moments;
import advantra.shapes.Sphere;
import advantra.swc.ExportSWC;
import advantra.tools.MeanShift3DSphere;
import advantra.trace.TinyBranchTrace.CharacteristicPoint;

import ij.ImagePlus;
import ij.io.FileSaver;

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
	
	public void trace(double x, double y, double z){ // trace from the selected point (coord.)
		
		int RES = 32;
		int MS_PTS = 100;
		int MS_PTS_TH = 25;
		int MS_MAX_ITER = 100;
		double MS_EPS = 0.001; 
		double MS_NEIGHBOUR = 0.01;
		double MS_ANGLE_RANGE_DEG = 30;
		
		
		// refine the point
		double[] start_xyz = refineStartPoint(new double[]{x, y, z}, 20); // 20 radius of refinement
		System.out.format("\nrefined start points at (%f,%f,%f)... \n" +
				"------------------------------------\n", 
				start_xyz[0],
				start_xyz[1],
				start_xyz[2]);
		
		// detect best matching hypothesis at refined start point
		System.out.println("calculating start hypothesis...");
		Hypothesis start_hyp = (new TinyBranchTrace(traced_img)).hyp_at(start_xyz);
		start_hyp.print();
		
		Sphere sph = new Sphere(start_xyz, start_hyp.getHypothesisRadius());
		ImagePlus sphere_img = sph.extract(traced_img, RES);
		
		MeanShift3DSphere ms3dSph = new MeanShift3DSphere(sphere_img, ((MS_ANGLE_RANGE_DEG/180)*Math.PI), MS_PTS);
		
		long t11 = System.currentTimeMillis();
		ms3dSph.run(MS_MAX_ITER, MS_EPS);
		double[][] out_dirs = ms3dSph.extractDirections(MS_NEIGHBOUR, MS_PTS_TH); // Mx3
		long t12 = System.currentTimeMillis();
		System.out.format(
				"at sphere radius R = %2.1f found %d dirs, elapsed %f sec. \n", 
				start_hyp.getHypothesisRadius(), out_dirs.length, (t12-t11)/1000f);
		
		if(true){return;}
		
		// detect directions
		double relative_dist_check_bifurcations = TinyBranch.check_bifurcations*start_hyp.getNeuriteRadius();
//		int 	spherical_image_resolution = 32;
//		double[][] start_directions = calculateTraceDirections(
//				start_xyz, 
//				relative_dist_check_bifurcations, 
//				RES);
		
		if(out_dirs==null || out_dirs.length==0){
			System.out.println("Trace stopped... " +
					"no directions detected at selected point.");
			return;
		}
		else{
			System.out.format("\nTrace starting at(%f,%f,%f)... %d directions detected at starting point:\n", 
					start_xyz[0],	start_xyz[1],	start_xyz[2],
					out_dirs.length);
			System.out.println("out_dirs: "+out_dirs.length+" x "+out_dirs[0].length);
			
		}
		
		// add them to hypothesis queue
		double[] branch_start_xyz = new double[3];
		for (int i = 0; i < out_dirs.length; i++) {
			branch_start_xyz[0] = start_xyz[0]+relative_dist_check_bifurcations*out_dirs[i][0];
			branch_start_xyz[1] = start_xyz[1]+relative_dist_check_bifurcations*out_dirs[i][1];
			branch_start_xyz[2] = start_xyz[2]+relative_dist_check_bifurcations*out_dirs[i][2];
			Hypothesis to_add = new Hypothesis(start_xyz, out_dirs[i], start_hyp.getNeuriteRadius(), start_hyp.getK());
			branch_queue.add(to_add);
			branch_queue_mother_idx.add(-1);
		}
		
		long t22 = System.currentTimeMillis();
		System.out.println("ms took : "+((t22-t11)/1000f)+" sec.");	
		System.out.println("returning..."); if(true){return;}
		
		/*
		 * loop
		 */
		
		while(branch_queue.size()>0 && traced_branch_nr<max_branch_nr){
			
			//take first hypothesis
			Hypothesis hyp_from_queue = branch_queue.get(0);
			branch_queue.remove(0);
			
			// take its mother idx
			branches_mother_idxs[traced_branch_nr] = branch_queue_mother_idx.get(0);
			branch_queue_mother_idx.remove(0);
			
			System.out.format("\n------------------------------------------------------\n");
			System.out.format("branch # %4d  |  %4d  remaining   |  limit: %4d \n", traced_branch_nr+1, branch_queue.size(), max_branch_nr);
			
			branches[traced_branch_nr].setStart(hyp_from_queue);
			
			/*
			 * MAIN STUFF: trace it
			 */
			
			CharacteristicPoint p = getCurrentlyTracedBranch().trace();// will calculate 'new_directions' 'new_seeds'
			System.out.println(p+" detected...");
			
			if(p==CharacteristicPoint.BRANCH){ // add them to queue
				
				// make hypotheses for further tracing
				for (int b = 0; b < getCurrentlyTracedBranch().numberOfNewSeeds(); b++) {
					branch_queue.add(new Hypothesis(getCurrentlyTracedBranch().getNewSeed(b),
													getCurrentlyTracedBranch().getNewDirection(b),
													getCurrentlyTracedBranchNeuriteRad(),
													getCurrentlyTracedBranch().getK()));
					branch_queue_mother_idx.add(traced_branch_nr);
				}
				System.out.println("\nadded "+getCurrentlyTracedBranch().numberOfNewSeeds()+" branches");
				
			}
			
			traced_branch_nr++;
			
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
	
	private double[] getCurrentlyTracedBranchPos(){
		return getCurrentlyTracedBranch().getCurrentPos();
	}
	
	private double getCurrentlyTracedBranchNeuriteRad(){
		return getCurrentlyTracedBranch().getCurrentRadius();
	}	

	private double[] AplusKtimesB(double[] A, double K, double[] B){
		double[] c = new double[3];
		c[0] = A[0]+K*B[0];
		c[1] = A[1]+K*B[1];
		c[2] = A[2]+K*B[2];
		return c;
	}
	
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

	private double[][] calculateTraceDirections(double[] at_pos, double at_r, int at_resolution){
		
		Sphere sp = new Sphere(at_pos[0], at_pos[1], at_pos[2], at_r);
		
		ImagePlus output 		= sp.extract(traced_img, at_resolution);	
		
		double angle_range_deg 	= 20;
		double angle_range_rad 	= (angle_range_deg/180)*Math.PI;
		int 	N 				= 200;
		int		N_th			= 20;
		MeanShift3DSphere ms3dSph = new MeanShift3DSphere(output, angle_range_rad, N);
		int max_iter 			= 200;
		double epsilon 			= 0.00001;
		ms3dSph.run(max_iter, epsilon);
		double[][] out_dirs = ms3dSph.extractDirections(0.05, N_th);
		
		(new FileSaver(output)).saveAsTiffStack(String.format("sphere_r_%f_%d_dirs.tif", at_r, out_dirs.length));
		
		return out_dirs;
		
	}
	
}