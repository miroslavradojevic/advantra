package advantra.trace;

import java.util.ArrayList;
import java.util.Vector;

import advantra.file.ExportSWC;
import advantra.general.ImageConversions;
import advantra.processing.IntensityCalc;
import advantra.tools.MAP;

import ij.ImagePlus;

public class NeuronTrace {

	IntensityCalc img_calc;
	ImagePlus     img_rgb;
	
	int 		max_branch_nr;
	int 		traced_branch_nr;
	
	int			looped_back 				= 0;	// branch that went back to the trace	
	int			default_branch_type 		= 2; 	// body branch
	int			terminal_branch 			= 3; 	// branch that finished with endpoint
	int 		end_with_posteriors_zero 	= 4;	// weird end
	
	ArrayList<Integer>			branches_mother_idxs;	 
	ArrayList<TinyBranchTrace> 	branches;
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
		branches 				= new ArrayList<TinyBranchTrace>(max_branch_nr); 
		branch_queue			= new ArrayList<Hypothesis>();
		branch_queue_mother_idx	= new ArrayList<Integer>();
		
	}
	
	public void 			trace(double x, double y, double z){ // trace from the selected point
		
		boolean manual_start				= true;
		Hypothesis	take_hypothesis 		= new Hypothesis();
		Vector<Hypothesis> new_hypotheses 	= new Vector<Hypothesis>();
		int CPU_NR = Tracing.CPU_NR;
		double[] start_xyz = new double[]{x, y, z}; 
		
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
					
					branches.add(traced_branch_nr, new TinyBranchTrace());
					//branches.set(traced_branch_nr, );
					branches.get(traced_branch_nr).setStart(take_hypothesis);
					branches_mother_idxs.add(traced_branch_nr, -1);//[traced_branch_nr] = -1;
					branches_types.add(traced_branch_nr, default_branch_type);
					/* 
					 * detection of the new directions (new_hypotheses assigned!)
					 */
					
					branches.get(traced_branch_nr).calculateNewSeeds(img_calc); // new_seeds calculated inside
					new_hypotheses = branches.get(traced_branch_nr).getNewSeedHypotheses(img_calc, manual_start); // will eliminate some seeds
					
					manual_start = false;
					System.out.print("}");
					
				}
				else{ // take_hypothesis is null
					
					System.out.println("Sum of posteriors at start point was close to zero! Quitting... ");
					return;
					
				}

			}
			else{ // not a manual start
				
				System.out.print("\nbranch #"+traced_branch_nr+"/"+(max_branch_nr-1)+"+"+extend+"<");
				
				//take first hypothesis if not the first execution
				take_hypothesis = branch_queue.get(0);	
				branch_queue.remove(0);
				
				//else{
				TinyBranchTrace to_add = new TinyBranchTrace();
				to_add.setStart(take_hypothesis);
				branches.add(traced_branch_nr, to_add);
				//branches.get(traced_branch_nr).setStart(take_hypothesis); 
				
				branches_mother_idxs.add(traced_branch_nr, branch_queue_mother_idx.get(0));//[traced_branch_nr] = branch_queue_mother_idx.get(0);
				branches_types.add(traced_branch_nr, default_branch_type);
				branch_queue_mother_idx.remove(0);
				
				// trace from the point on
				boolean enable_bif_detection = false;
				int loop_count = 0;
				while(loop_count<Tracing.N-1){ 
					
					System.out.print("=");
					
					MAP.atHyp(take_hypothesis);
					ITERS  = MAP.total_hyps; // ITERS has to be multiple of CPU_NR
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
						
						/*
						 * save it
						 */
						
						branches.get(traced_branch_nr).storeHypothesis(take_hypothesis);
						
						/*
						 * check 
						 */
						
//						double dst = 0.8*TinyBranch.check_bifurcations*current_branch.getCurrentTraceHypothesisRadius();
//						double dx = current_branch.getCurrentPosX()-current_branch.getSeedPosX();
//						double dy = current_branch.getCurrentPosY()-current_branch.getSeedPosY();
//						double dz = current_branch.getCurrentPosZ()-current_branch.getSeedPosZ();
						
						if(loop_count>3 && !enable_bif_detection)	// (dx*dx+dy*dy+dz*dz>dst*dst)
							enable_bif_detection = true;
						
						if(enable_bif_detection){
							
							// check if it was traced
							if(isTraced(take_hypothesis.getPosition())){
								
								System.out.print("SKIP_THE_REST(TRACED)");
								branches_types.set(traced_branch_nr, looped_back);
								extend++;
								// remove its mother index as well
								//branch_queue_mother_idx.remove(0);
								new_hypotheses.clear();
								break;
								
							}

							
							System.out.print("*");
							
							/*
							 *  detection
							 */	
							
							branches.get(traced_branch_nr).calculateNewSeeds(img_calc); // new_seeds calculated inside
							new_hypotheses = branches.get(traced_branch_nr).getNewSeedHypotheses(img_calc, manual_start); // will eliminate some seeds
							
							if(new_hypotheses.size()!=1){
								if(new_hypotheses.size()==1){System.out.println("SUSPICIOUS");}
								if(new_hypotheses.size()==0){branches_types.set(traced_branch_nr, terminal_branch);}
								break; // break while loop : 2-branch, 0-ms-gave suspicious result
							}
							
						}
						else{
							System.out.print("#");
						}

					}
					else{
						System.out.print("SUM_POSTS_ZERO");
						new_hypotheses.clear();
						branches_types.set(traced_branch_nr, end_with_posteriors_zero);
						break;// tracing itself got 0 posterior sum
					}
					
					loop_count++;
				
				} // while loop
				
				if(loop_count==Tracing.N-1) System.out.print("LOOP_LIMIT_REACHED");
				
				System.out.print(">");
			
			//}
				
			} // not a manual start
			
			/*
			 *  align new seed hypotheses with respect to the image using MAP
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

				new_hypotheses.setElementAt(MAP.takeMap(), k);
				
				// add it to queue
				branch_queue.add(new_hypotheses.get(k));
				branch_queue_mother_idx.add(traced_branch_nr);
				
			}
			
			int nbr = new_hypotheses.size();
				
			System.out.println(nbr+" new branches added to queue, ("+branch_queue.size()+" currently pending)");
			
			traced_branch_nr++; // count next one once this is finished
			
		} // while loop
		
		System.out.print("reconstructed "+traced_branch_nr+" branches, ");
		System.out.println("\n"+branch_queue.size()+" branches waiting in the queue");
		
	}		
	
	public boolean 			isTraced(double[] pos){
		for (int i = 0; i < traced_branch_nr; i++) {
			for (int j = 0; j < branches.get(i).getCount(); j++) {
				
				double c0 = branches.get(i).centerlines[j][0];
				double c1 = branches.get(i).centerlines[j][1];
				double c2 = branches.get(i).centerlines[j][2];
				
				double r  = branches.get(i).radiuses[j];
				
				double d2 = Math.pow(pos[0]-c0, 2) +
							Math.pow(pos[1]-c1, 2) +
							Math.pow(pos[2]-c2, 2);
				
				if(d2<=r*r){
					return true;
				}

			}
		}
		
		return false;
	}
	
	public void 			export_swc(String file_name){
		
		for (int i = 0; i < traced_branch_nr; i++) {
			System.out.println(
					"branch "+i+", size "+branches.get(i).getCount()+", mother idx: "+branches_mother_idxs.get(i));
		}
		
		int[] last_line = new int[traced_branch_nr];
		
		ExportSWC swc_file = new ExportSWC(file_name);
		
		int line = 1; // current line in swc file
		int prev_line = -1;
		
		for (int i = 0; i < traced_branch_nr; i++) {
			boolean first = true;
			for (int j = 0; j < branches.get(i).getCount(); j++) {
				
				if(first){
					if(branches_mother_idxs.get(i)==(-1)) 
						prev_line = -1;
					else								
						prev_line = last_line[branches_mother_idxs.get(i)];

					first = false;
				}

				int branch_color = branches_types.get(i);//i%8;
				
				//if(branch_color!=0){
					swc_file.writelnSWC(String.format("%5d %2d %6.2f %6.2f %6.2f %6.2f %5d", 
							line, 
							branch_color,
							branches.get(i).centerlines[j][1],
							branches.get(i).centerlines[j][0],
							branches.get(i).centerlines[j][2],
							branches.get(i).radiuses[j],
							prev_line));

				//}
				
				prev_line = line;
				line++;
				
			}
			
			last_line[i] = prev_line;
			
		}
		System.out.println("finished... closing");
		swc_file.closeSWC();
		
	}
	
}