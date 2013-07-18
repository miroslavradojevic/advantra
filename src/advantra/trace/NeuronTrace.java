package advantra.trace;

import java.util.ArrayList;
import java.util.Vector;

import advantra.file.ExportSWC;
import advantra.general.ImageConversions;
import advantra.processing.IntensityCalc;
import advantra.shapes.Sphere;
import advantra.tools.MAP;
import advantra.tools.MeanShift3DSphere;

import ij.ImagePlus;

public class NeuronTrace {

	IntensityCalc img_calc;
	ImagePlus     img_rgb;
	
	int 		max_branch_nr;
	int 		traced_branch_nr;
	
	int			ms_null 					= 4;	// branch that went back to the trace	
	int			no_in 						= 1; 	// body branch
	int			body 						= 2; 	// body branch
	int			end 						= 3; 	// branch that finished with endpoint
	int 		looped_back					= 0;	// went back to already traced space
	int 		map_zero 					= 5;	// weird end
	int			branch						= 6;
	
	ImagePlus 			sphere_img; 
	public double[][] 	new_seeds; 					// br_dirs x3 (because they're 3d space coordinates)

	
	ArrayList<Integer>			branches_mother_idxs;	 
	ArrayList<BranchTrace> 		branches;
	ArrayList<Integer>			branches_types;
	
	ArrayList<Hypothesis>		branch_queue;
	ArrayList<Integer>			branch_queue_mother_idx;
	
	public NeuronTrace(ImagePlus im, int max_branches){
		
		img_calc = new IntensityCalc(im.getStack());
		img_rgb  = ImageConversions.ImagePlusToRGB(im);
		
		max_branch_nr 		= max_branches;
		traced_branch_nr 	= 0;
		
		branches_mother_idxs 	= new ArrayList<Integer>(max_branch_nr);
		branches_types			= new ArrayList<Integer>(max_branch_nr);
		branches 				= new ArrayList<BranchTrace>(max_branch_nr); 
		branch_queue			= new ArrayList<Hypothesis>();
		branch_queue_mother_idx	= new ArrayList<Integer>();
		
	}
	
	public void 			trace(double x, double y, double z){ // trace from the selected point
		
		boolean manual_start				= true;
		Hypothesis	take_hypothesis 		= new Hypothesis();
		Vector<Hypothesis> new_hypotheses 	= new Vector<Hypothesis>();
		int CPU_NR = Tracing.CPU_NR;
		double[] start_xyz = new double[]{x, y, z}; 
		double rr;	// sphere radius to check bifurcations at
		
		MAP.loadImage(img_calc);
		MAP map_jobs[] = new MAP[CPU_NR];
		int ITERS;
		int extend = 0;
		
		while((manual_start || branch_queue.size()>0) && traced_branch_nr<max_branch_nr+extend){
			
			if(manual_start){
				
				System.out.print("{");
				
				MAP.atPoint(start_xyz);
				ITERS  = MAP.total_hyps; // ITERS has to be multiple of 8
				for (int i = 0; i < map_jobs.length; i++) {
					map_jobs[i] = new MAP(i * ITERS/CPU_NR, (i+1) * ITERS/CPU_NR);
					map_jobs[i].start();
				}
				for (int i = 0; i < map_jobs.length; i++) {
					try {
						map_jobs[i].join();
					} catch (InterruptedException e) {
						e.printStackTrace();
					}
				}
				
				take_hypothesis = MAP.takeMap();
				
				if(take_hypothesis!=null){
					
					branches.add(traced_branch_nr, new BranchTrace());
					branches.get(traced_branch_nr).setStart(take_hypothesis);
					branches_mother_idxs.add(traced_branch_nr, -1);
					branches_types.add(traced_branch_nr, body);
					/* 
					 * detection of the new directions (new_hypotheses assigned!)
					 */
					
					rr = take_hypothesis.getNeuriteRadius();
					rr = (rr<=2)? (3*rr) : (3*rr);
					new_seeds = extractDirections(rr);
					
					if(new_seeds==null){ // ms did not converge
						System.out.print("MS_NULL");
						new_hypotheses.clear();
						branches_types.set(traced_branch_nr, ms_null);
						break;
					}
					
					new_hypotheses = extractSeeds(manual_start);
					
					manual_start = false;
					System.out.print("}");
					
				}
				else{ 
					// all posteriors were zero - terminate
					System.out.println("Sum of posteriors at start point was close to zero! Quitting trace... ");
					return;
					
				}

			}
			else{ // not a manual start
				
				System.out.print("\nbranch #"+traced_branch_nr+"/"+(max_branch_nr-1)+"+"+extend+"<");
				
				//take hypothesis from the queue
				take_hypothesis = branch_queue.get(0);	
				branch_queue.remove(0);
				
				BranchTrace to_add = new BranchTrace();
				to_add.setStart(take_hypothesis);
				branches.add(traced_branch_nr, to_add);
				
				branches_mother_idxs.add(traced_branch_nr, branch_queue_mother_idx.get(0));
				branch_queue_mother_idx.remove(0);
				
				branches_types.add(traced_branch_nr, body);
				
				// trace from the point on
				boolean enable_bif_detection = false;
				int loop_count = 0;
				
				while(loop_count<Tracing.N-1){ 
					
					System.out.print("=");
					
					/*
					 * one trace iteration
					 */
					
					MAP.atHyp(take_hypothesis);
					ITERS  = MAP.total_hyps; 
					for (int i = 0; i < map_jobs.length; i++) {
						map_jobs[i] = new MAP(i * ITERS/CPU_NR, (i+1) * ITERS/CPU_NR);
						map_jobs[i].start();
					}
					for (int i = 0; i < map_jobs.length; i++) {
						try {
							map_jobs[i].join();
						} catch (InterruptedException e) {
							e.printStackTrace();
						}
					}
					take_hypothesis = MAP.takeMap();
					
					// posterior almost zero
					if(take_hypothesis==null){
						System.out.print("MAP_ZERO");
						new_hypotheses.clear();
						branches_types.set(traced_branch_nr, map_zero);
						break;// tracing itself got 0 posterior sum
					}
					
					// check if it was traced
					if(isTraced(take_hypothesis.getPosition())){
						System.out.print("IS_TRACED");
						branches_types.set(traced_branch_nr, looped_back);
						extend++;
						// remove its mother index as well
						//branch_queue_mother_idx.remove(0);
						new_hypotheses.clear();
						break;
					}
					
					if(take_hypothesis!=null){
						/*
						 * save it
						 */
						branches.get(traced_branch_nr).storeHypothesis(take_hypothesis);
						
						double dx = branches.get(traced_branch_nr).getCurrentPosX()-
								branches.get(traced_branch_nr).getSeedPosX();
						double dy = branches.get(traced_branch_nr).getCurrentPosY()-
								branches.get(traced_branch_nr).getSeedPosY();
						double dz = branches.get(traced_branch_nr).getCurrentPosZ()-
								branches.get(traced_branch_nr).getSeedPosZ();
						
						rr = take_hypothesis.getNeuriteRadius();
						rr = (rr<=2)? (3*rr) : (3*rr);
						
						if((dx*dx+dy*dy+dz*dz>rr*rr) && !enable_bif_detection)	
							enable_bif_detection = true;
						
						if(enable_bif_detection){
														
							System.out.print("*");
							
							/*
							 *  check bifs
							 */	
							
							new_seeds = extractDirections(rr);
							
							if(new_seeds==null){ // ms did not converge
								System.out.print("MS_NULL");
								new_hypotheses.clear();
								branches_types.set(traced_branch_nr, ms_null);
								break;
							}
							
							new_hypotheses = extractSeeds(manual_start);
							
							if(new_hypotheses.size()==0){
								System.out.print("NO_NEW_HYPS");
								new_hypotheses.clear();
								branches_types.set(traced_branch_nr, no_in);
								break;
							}
							
							if(new_hypotheses.size()==1){
								System.out.print("END");
								new_hypotheses.clear();
								branches_types.set(traced_branch_nr, end);
								break;
							}
							
							//if(new_hypotheses.size()==2){
								
								//System.out.print("BODY");
								
							//}
							
							if(new_hypotheses.size()>2){
								System.out.print("BRANCH");
								branches_types.set(traced_branch_nr, branch);
								break;
							}
							
						}
						else{
							System.out.print("#");
						}

					}

					
					loop_count++;
				
				} // while loop
				
				if(loop_count==Tracing.N-1) System.out.print("LOOP_LIMIT_REACHED");
				
				System.out.print(">");
			
			//}
				
			} // not a manual start
			
			/*
			 *  align new seed hypotheses with respect to the image using MAP and add them
			 */
			
			for (int k = 0; k < new_hypotheses.size(); k++) {
				
				double[] new_hyp_pos = new_hypotheses.get(k).getPosition();
				double[] new_hyp_ort = new_hypotheses.get(k).getOrientation();
				
				MAP.atPoint(new_hyp_pos, new_hyp_ort);
				ITERS  = MAP.total_hyps;
				for (int i = 0; i < map_jobs.length; i++) {
					map_jobs[i] = new MAP(i * ITERS/CPU_NR, (i+1) * ITERS/CPU_NR);
					map_jobs[i].start();
				}
				for (int i = 0; i < map_jobs.length; i++) {
					try {
						map_jobs[i].join();
					} catch (InterruptedException e) {
						e.printStackTrace();
					}
				}
				
				Hypothesis aligned = MAP.takeMap();
				//new_hypotheses.setElementAt(, k);
				
				// add it to queue if != than null
				if(aligned!=null){
					branch_queue.add(aligned);
					branch_queue_mother_idx.add(traced_branch_nr);
				}
				
			}
			
			int nbr = new_hypotheses.size();
				
			System.out.println(nbr+" new branches added to queue, ("+branch_queue.size()+" currently pending)");
			
			traced_branch_nr++; // count next one once this is finished
			
		} // while loop
		
		System.out.print("reconstructed "+(traced_branch_nr-extend)+" branches, ");
		System.out.println("\n"+branch_queue.size()+" branches waiting in the queue");
		
	}
	
	private double[][] 		extractDirections(double rd){
		
		// extracts array with new directions detected by mean-shift
		double[] p 					= branches.get(traced_branch_nr).getCurrentTraceHypothesisCenterpoint();
		Sphere sphere_from_image 	= new Sphere(p, rd);
		sphere_img     				= sphere_from_image.extract(img_calc, Tracing.extract_sphere_resolution); 
		
		MeanShift3DSphere.load(
				sphere_img, 
				sphere_from_image, 
				Tracing.MS_ANGLE_RANGE_RAD, 
				Tracing.MS_PTS, 
				Tracing.MS_MAX_ITER, 
				Tracing.MS_EPS);
		MeanShift3DSphere ms_jobs[] = new MeanShift3DSphere[Tracing.CPU_NR];
		for (int i = 0; i < ms_jobs.length; i++) {
			ms_jobs[i] = new MeanShift3DSphere(i*Tracing.MS_PTS/Tracing.CPU_NR,  (i+1)*Tracing.MS_PTS/Tracing.CPU_NR);
			ms_jobs[i].start();
		}
		for (int i = 0; i < ms_jobs.length; i++) {
			try {
				ms_jobs[i].join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		MeanShift3DSphere.extractClusters(Tracing.MS_NEIGHBOUR_RAD, Tracing.MS_PTS_TH); 
		
		return MeanShift3DSphere.cluster_seed;
		
	}
	
	private Vector<Hypothesis>		extractSeeds(boolean isManualStart){
		
		Vector<Hypothesis> new_hyps = new Vector<Hypothesis>();

		//boolean expelled = false; 
		
		for (int i = 0; i < new_seeds.length; i++) {
			
			if(new_seeds[i]!=null){
				
				double[] S_vec = new double[3]; // seed orientation
				S_vec[0] = new_seeds[i][0]-branches.get(traced_branch_nr).getCurrentPosX();
				S_vec[1] = new_seeds[i][1]-branches.get(traced_branch_nr).getCurrentPosY();
				S_vec[2] = new_seeds[i][2]-branches.get(traced_branch_nr).getCurrentPosZ();
				
				// not necessary to normalize
				double S_vec_norm = Math.sqrt(S_vec[0]*S_vec[0]+S_vec[1]*S_vec[1]+S_vec[2]*S_vec[2]);
				S_vec[0] /= S_vec_norm;
				S_vec[1] /= S_vec_norm;
				S_vec[2] /= S_vec_norm;
				
				double r = branches.get(traced_branch_nr).getCurrentTraceHypothesisNeuriteRadius();
				
				Hypothesis hyp_to_add = new Hypothesis(new_seeds[i], S_vec, r, Tracing.k);
				
				new_hyps.add(hyp_to_add);
				
			}
			
		}
		
		return new_hyps;
		
	}
	
	
	public boolean 			isTraced(double[] pos){
		for (int i = 0; i < traced_branch_nr; i++) {
			for (int j = 0; j < branches.get(i).count; j++) {
				
				double c0 = branches.get(i).centerlines[j][0];
				double c1 = branches.get(i).centerlines[j][1];
				double c2 = branches.get(i).centerlines[j][2];
				
				double r  = branches.get(i).radiuses[j];
				
				double d2 = Math.pow(pos[0]-c0, 2) +
							Math.pow(pos[1]-c1, 2) +
							Math.pow(pos[2]-c2, 2);
				
				if(d2<=Tracing.RADIUS_SCALE_BORDER*Tracing.RADIUS_SCALE_BORDER*r*r){
					return true;
				}

			}
		}
		
		return false;
	}
	
	public void 			export_swc(String file_name){
		
		for (int i = 0; i < traced_branch_nr; i++) {
			System.out.println(
					"branch "+i+", size "+branches.get(i).count+", mother idx: "+branches_mother_idxs.get(i));
		}
		
		int[] last_line = new int[traced_branch_nr];
		
		ExportSWC swc_file = new ExportSWC(file_name);
		
		int line = 1; // current line in swc file
		int prev_line = -1;
		
		for (int i = 0; i < traced_branch_nr; i++) {
			boolean first = true;
			for (int j = 0; j < branches.get(i).count; j++) {
				
				if(first){
					if(branches_mother_idxs.get(i)==(-1)) 
						prev_line = -1;
					else								
						prev_line = last_line[branches_mother_idxs.get(i)];

					first = false;
				}

				int branch_color = 6;//branches_types.get(i);//i%8;
				
				if(branch_color!=looped_back && branch_color!=ms_null && branch_color!=map_zero){//&& branches.get(i).count>1
					swc_file.writelnSWC(String.format("%5d %2d %6.2f %6.2f %6.2f %6.2f %5d", 
							line, 
							branch_color,
							branches.get(i).centerlines[j][1]+10,
							branches.get(i).centerlines[j][0]+10,
							branches.get(i).centerlines[j][2],
							branches.get(i).radiuses[j],
							prev_line));

				}
				
				prev_line = line;
				line++;
				
			}
			
			last_line[i] = prev_line;
			
		}
		System.out.println("finished... closing");
		swc_file.closeSWC();
		
	}
	
}

//if(!isManualStart){
///* check if it belongs to the current branch already
// * chk checkings at max
// */
//
////int chk = 0;
//
//for (int k = branches.get(traced_branch_nr).count-1; k >=0; k--) {
//	
//	
//	
//	
//	// x1 and x2 will be centers
//	double[] x1 = branches.get(traced_branch_nr).centerlines[k];
//	double   r1 = branches.get(traced_branch_nr).radiuses[k];
////	double[] x2 = branches.get(traced_branch_nr).centerlines[k-1];
////	double   r2 = branches.get(traced_branch_nr).radiuses[k-1];
//	double[] x0 = new_seeds[i];
//	
//	double d2 = Math.pow(x0[0]-x1[0], 2)+Math.pow(x0[1]-x1[1], 2)+Math.pow(x0[2]-x1[2], 2);
//	
//	if(d2<=r1*r1){
//		
//	
//			
////	chk++;
////	
////	if(chk>Tracing.STEPS_BACKWARDS){
////		break;
////	}
//	
////	double in1 = 
////			(x0[0]-x1[0])*(x2[0]-x1[0])+
////			(x0[1]-x1[1])*(x2[1]-x1[1])+
////			(x0[2]-x1[2])*(x2[2]-x1[2]);
////	
////	double in2 = 
////			(x0[0]-x2[0])*(x1[0]-x2[0])+
////			(x0[1]-x2[1])*(x1[1]-x2[1])+
////			(x0[2]-x2[2])*(x1[2]-x2[2]);
//	
//	
////	double c1dist = 
////			Math.pow(x1[0]-branches.get(traced_branch_nr).getCurrentPosX(),2)+
////			Math.pow(x1[1]-branches.get(traced_branch_nr).getCurrentPosY(),2)+
////			Math.pow(x1[2]-branches.get(traced_branch_nr).getCurrentPosZ(),2);
////	
////	double c2dist = 
////			Math.pow(x2[0]-branches.get(traced_branch_nr).getCurrentPosX(),2)+
////			Math.pow(x2[1]-branches.get(traced_branch_nr).getCurrentPosY(),2)+
////			Math.pow(x2[2]-branches.get(traced_branch_nr).getCurrentPosZ(),2);
////	
////	double sddist = 
////			Math.pow(new_seeds[i][0]-branches.get(traced_branch_nr).getCurrentPosX(),2)+
////			Math.pow(new_seeds[i][1]-branches.get(traced_branch_nr).getCurrentPosY(),2)+
////			Math.pow(new_seeds[i][2]-branches.get(traced_branch_nr).getCurrentPosZ(),2);
////	
////	System.out.format("seed%d: %d:(%.2f, %.2f, %.2f)\n",i, k,c1dist,sddist,c2dist);
////	
////	//boolean x0_btw = ;
////	
////	boolean x0_btw = c2dist>sddist;
////	
////	if(x0_btw && ((in1>0)&&(in2>0))){
////		
////		//System.out.println("centers "+k+" and "+(k-1)+"are capturing seed!");
////		double x10_2 = Math.pow(x1[0]-x0[0], 2)+Math.pow(x1[1]-x0[1], 2)+Math.pow(x1[2]-x0[2], 2);
////		double x21_2 = Math.pow(x2[0]-x1[0], 2)+Math.pow(x2[1]-x1[1], 2)+Math.pow(x2[2]-x1[2], 2);
////		double x10_dot_x21 = 
////				(x1[0]-x0[0])*(x2[0]-x1[0])+
////				(x1[1]-x0[1])*(x2[1]-x1[1])+
////				(x1[2]-x0[2])*(x2[2]-x1[2]);
////		
////		double d2  = (x10_2*x21_2 - x10_dot_x21*x10_dot_x21)/x21_2;
////
////		if(d2<=Math.pow(Tracing.RADIUS_SCALE_BORDER,2)*Math.pow(r1,2)
////				|| d2<=Math.pow(Tracing.RADIUS_SCALE_BORDER,2)*Math.pow(r2,2)){
//		new_seeds[i] = null;
//			//System.out.println("seed "+i+" expelled");
//			//expelled = true;
//		break;
//	}
//}
//}