package advantra.trace;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Vector;

import advantra.file.ExportSWC;
import advantra.processing.IntensityCalc;

import ij.ImagePlus;
import ij.io.FileSaver;

public class NeuronTrace {

	IntensityCalc img_calc;
	
	int 		max_branch_nr;
	int 		traced_branch_nr;			 
	int[]		branches_mother_idxs;	 
	
	TinyBranchTrace[] 			branches;
	
	ArrayList<Hypothesis>		branch_queue;
	ArrayList<Integer>			branch_queue_mother_idx;
	
	public NeuronTrace(ImagePlus im, int max_branches){
		
		img_calc = new IntensityCalc(im.getStack());
		
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
		
		boolean manual_start			= true;
		Hypothesis	take_hypothesis 	= new Hypothesis();
		Hypothesis  estimated_hyp 		= new Hypothesis();
		Vector<Hypothesis> new_hypotheses = new Vector<Hypothesis>();
		//double[][] detected_seeds       = null;
		TinyBranchTrace current_branch 	= new TinyBranchTrace();
		//ImagePlus sphere_img 			= new ImagePlus();
		
		double[] start_xyz = new double[]{x, y, z}; 
		
		while((manual_start || branch_queue.size()>0) && traced_branch_nr<max_branch_nr){
			
			if(manual_start){
				
				System.out.print("{CLICK!");
				// detect best matching hypothesis at refined start point
				take_hypothesis = current_branch.hyp_at(start_xyz, img_calc); 
				current_branch 	= getCurrentlyTracedBranch();
				current_branch.setStart(take_hypothesis);
				branches_mother_idxs[traced_branch_nr] = -1;
				
				/* 
				 * detection of the new directions (new_hypotheses assigned!)
				 */
				// new seeds calculated inside
				current_branch.calculateNewSeeds(img_calc);
				new_hypotheses = current_branch.getNewSeedHypotheses(img_calc); // will eliminate some seeds
				manual_start = false;
				System.out.print("}");
				
				
			}
			else{ // not a manual start
				
				System.out.print("\nbranch #"+traced_branch_nr+"/"+(max_branch_nr-1)+"<");
				//take first hypothesis if not the first execution
				take_hypothesis = branch_queue.get(0);	
				branch_queue.remove(0);
				
				current_branch = getCurrentlyTracedBranch();
				
				current_branch.setStart(take_hypothesis); 
				branches_mother_idxs[traced_branch_nr] = branch_queue_mother_idx.get(0);
				branch_queue_mother_idx.remove(0);
				
				// trace from the point on
				boolean enable_bif_detection = false;
				int loop_count = 0;
				while(loop_count<TinyBranch.N-1){ 
					
					current_branch.predict();
					current_branch.calculatePriors();
					current_branch.calculateLikelihoods(img_calc);
					
					double sum_posts = current_branch.calculatePosteriors();
					
					System.out.print("=");
					
					if(sum_posts>0){
						
						estimated_hyp 	= current_branch.getMaximumPosteriorHypothesis();
						
						/*
						 * save it
						 */
						
						current_branch.storeHypothesis(estimated_hyp, false);
						
						/*
						 * check if it made loop
						 */
						
						double dst = 0.8*TinyBranch.check_bifurcations*current_branch.getCurrentTraceHypothesisRadius();
						double dx = current_branch.getCurrentPosX()-current_branch.getSeedPosX();
						double dy = current_branch.getCurrentPosY()-current_branch.getSeedPosY();
						double dz = current_branch.getCurrentPosZ()-current_branch.getSeedPosZ();
						
						if((dx*dx+dy*dy+dz*dz>dst*dst) && !enable_bif_detection)	
							enable_bif_detection = true; // || loop_count>3
						
						if(enable_bif_detection){
							
							System.out.print("*");
							
							/*
							 *  detection
							 */							
							current_branch.calculateNewSeeds(img_calc);
							new_hypotheses.clear();
							new_hypotheses = current_branch.getNewSeedHypotheses(img_calc); // will eliminate some seeds
							
							if(new_hypotheses.size()>1) {
								
								
								
								/*
								 * save sphere image
								 */
								String file_name = String.format("b%dl%d.tif", traced_branch_nr, loop_count);
								FileSaver fs = new FileSaver(current_branch.getSphereImage());
								fs.saveAsTiffStack(file_name);
								/*
								 * save convergence positions
								 */
								FileWriter fw; 
								try{
									fw = new FileWriter(file_name+"_converg.csv");
									
									double[][] conv_points = current_branch.getTcartesian();
									if(conv_points!=null){
										for (int i = 0; i < conv_points.length; i++) {
										
											fw.write(String.format("%f, %f, %f\n", 
													conv_points[i][1]+1,
													conv_points[i][0]+1,
													conv_points[i][2]+1));
											
										}
									}
									fw.close();
									
									fw = new FileWriter(file_name+"_extracted_seeds.csv");
									double[][] seed_points = current_branch.getNewSeeds();
									
									if(seed_points!=null){
									for (int i = 0; i < seed_points.length; i++) {
										if(seed_points[i]!=null){
											for (int j = 0; j < seed_points[i].length; j++) {
												fw.write(""+seed_points[i][j]+", ");
											}
											fw.write("\n");
										}
									}
									}
									
									fw.close();
									
//									fw = new FileWriter(file_name+"_extracted_seeds_test.csv");
//									double[][] seed_points_test = current_branch.getNewSeedsTest();
//									
//									if(seed_points_test!=null){
//									for (int i = 0; i < seed_points_test.length; i++) {
//										if(seed_points_test[i]!=null){
//											for (int j = 0; j < seed_points_test[i].length; j++) {
//												fw.write(""+seed_points_test[i][j]+", ");
//											}
//											fw.write("\n");
//										}
//									}
//									}
//									
//									fw.close();
									
									
								}
								catch(IOException exIO){
									System.out.println("error writing.");
								}
								
								System.out.println("done writing.");
								
								break;
							}
						}
						else{
							System.out.print("#");
						}

					}
					else{
						System.out.print("SUM_POSTS_ZERO");
						break;// tracing itself got 0 posteriors
					}
					loop_count++;
				}
				
				if(loop_count==TinyBranch.N-1) System.out.print("LOOP_LIMIT_REACHED");
				
				System.out.print(">");
			
			}
			
			if(new_hypotheses!=null){ 
				int how_many_new_hyps = 0;
				for (int b = 0; b < new_hypotheses.size(); b++) {
					if(new_hypotheses.get(b)!=null){
						branch_queue.add(new_hypotheses.get(b));
						branch_queue_mother_idx.add(traced_branch_nr);
						how_many_new_hyps++;
					}
				}
				
				// current queue
				System.out.println(how_many_new_hyps+" new branches, ("+branch_queue.size()+" pending)");
			}
			
			traced_branch_nr++; // count next one once this is finished
			
			
		}
		
		System.out.print("reconstructed "+traced_branch_nr+" branches, ");
		System.out.println("\n"+branch_queue.size()+" branches waiting in the queue");
		
		
	}		
	
	public void 			export_swc(String file_name){
		
		System.out.println("*-------");
		for (int i = 0; i < traced_branch_nr; i++) {
			System.out.println(
					"branch "+i+", size "+branches[i].getCount()+", mother idx: "+branches_mother_idxs[i]);
			System.out.println("*-------");
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
	
	private TinyBranchTrace getCurrentlyTracedBranch(){
		return branches[traced_branch_nr];
	}
	
}