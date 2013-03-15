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
	int[]		branches_mother_idxs;	 
	
	TinyBranchTrace[] 			branches;
	
	ArrayList<Hypothesis>		branch_queue;
	ArrayList<Integer>			branch_queue_mother_idx;
	
	public NeuronTrace(ImagePlus im, int max_branches){
		
		img_calc = new IntensityCalc(im.getStack());
		img_rgb  = ImageConversions.ImagePlusToRGB(im);
		
		max_branch_nr 	= max_branches;
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
			branches[i] = new TinyBranchTrace();
		}
		
	}
	
	public void 			trace(double x, double y, double z){ // trace from the selected point
		
		boolean manual_start				= true;
		Hypothesis	take_hypothesis 		= new Hypothesis();
		Vector<Hypothesis> new_hypotheses 	= new Vector<Hypothesis>();
		//TinyBranchTrace current_branch 		= new TinyBranchTrace();
		
		double[] start_xyz = new double[]{x, y, z}; 
		
		MAP.loadImage(img_calc);
		int CPU_NR = 8;
		MAP map_jobs[] = new MAP[CPU_NR];
		int ITERS;
		
		while((manual_start || branch_queue.size()>0) && traced_branch_nr<max_branch_nr){
			
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
					
					//current_branch = branches[traced_branch_nr];
					branches[traced_branch_nr].setStart(take_hypothesis);
					branches_mother_idxs[traced_branch_nr] = -1;
					
					/* 
					 * detection of the new directions (new_hypotheses assigned!)
					 */
					
					branches[traced_branch_nr].calculateNewSeeds(img_calc); // new_seeds calculated inside
					new_hypotheses = branches[traced_branch_nr].getNewSeedHypotheses(img_calc, manual_start); // will eliminate some seeds
					
					manual_start = false;
					System.out.print("}");
					
				}
				else{ // take_hypothesis is null
					
					System.out.println("Sum of posteriors at start point was close to zero! Quitting... ");
					return;
					
				}

			}
			else{ // not a manual start
				
				System.out.print("\nbranch #"+traced_branch_nr+"/"+(max_branch_nr-1)+"<");
				
				//take first hypothesis if not the first execution
				take_hypothesis = branch_queue.get(0);	
				branch_queue.remove(0);
				
				// check if it was traced
				if(isTraced(take_hypothesis.getPosition())){
					
					System.out.print("SKIP(TRACED)");
					// remove its mother index as well
					branch_queue_mother_idx.remove(0);
					new_hypotheses.clear();
					
				}
				else{
				
				branches[traced_branch_nr].setStart(take_hypothesis); 
				branches_mother_idxs[traced_branch_nr] = branch_queue_mother_idx.get(0);
				branch_queue_mother_idx.remove(0);
				
				// trace from the point on
				boolean enable_bif_detection = false;
				int loop_count = 0;
				while(loop_count<TinyBranch.N-1){ 
					
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
						
						branches[traced_branch_nr].storeHypothesis(take_hypothesis);
						
						/*
						 * check 
						 */
						
//						double dst = 0.8*TinyBranch.check_bifurcations*current_branch.getCurrentTraceHypothesisRadius();
//						double dx = current_branch.getCurrentPosX()-current_branch.getSeedPosX();
//						double dy = current_branch.getCurrentPosY()-current_branch.getSeedPosY();
//						double dz = current_branch.getCurrentPosZ()-current_branch.getSeedPosZ();
						
						if(loop_count>2 && !enable_bif_detection)	// (dx*dx+dy*dy+dz*dz>dst*dst)
							enable_bif_detection = true;
						
						if(enable_bif_detection){
							
							System.out.print("*");
							
							/*
							 *  detection
							 */	
							
							branches[traced_branch_nr].calculateNewSeeds(img_calc); // new_seeds calculated inside
							new_hypotheses = branches[traced_branch_nr].getNewSeedHypotheses(img_calc, manual_start); // will eliminate some seeds
							
							if(new_hypotheses.size()!=1){
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
						break;// tracing itself got 0 posterior sum
					}
					
					loop_count++;
				
				} // while loop
				
				if(loop_count==TinyBranch.N-1) System.out.print("LOOP_LIMIT_REACHED");
				
				System.out.print(">");
			
			}
				
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
			for (int j = 0; j < branches[i].getCount(); j++) {
				
				double c0 = branches[i].centerlines[j][0];
				double c1 = branches[i].centerlines[j][1];
				double c2 = branches[i].centerlines[j][2];
				
				double r  = branches[i].radiuses[j];
				
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
					"branch "+i+", size "+branches[i].getCount()+", mother idx: "+branches_mother_idxs[i]);
		}
		
		int[] last_line = new int[traced_branch_nr];
		
		ExportSWC swc_file = new ExportSWC(file_name);
		
		int line = 1; // current line in swc file
		int prev_line = -1;
		
		for (int i = 0; i < traced_branch_nr; i++) {
			boolean first = true;
			for (int j = 0; j < branches[i].getCount(); j++) {
				
				if(first){
					if(branches_mother_idxs[i]==(-1)) 
						prev_line = -1;
					else								
						prev_line = last_line[branches_mother_idxs[i]];

					first = false;
				}

				int branch_color = i%8;
				
				swc_file.writelnSWC(String.format("%5d %2d %6.2f %6.2f %6.2f %6.2f %5d", 
						line, 
						branch_color,
						branches[i].centerlines[j][1],
						branches[i].centerlines[j][0],
						branches[i].centerlines[j][2],
						branches[i].radiuses[j],
						prev_line));
				
				prev_line = line;
				line++;
				
			}
			
			last_line[i] = prev_line;
			
		}
		System.out.println("finished... closing");
		swc_file.closeSWC();
		
	}
	
}