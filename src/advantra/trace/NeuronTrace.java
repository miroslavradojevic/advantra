package advantra.trace;

import java.util.ArrayList;

import advantra.file.ExportSWC;
import advantra.shapes.Sphere;
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
	
	public void 			trace(double x, double y, double z){ // trace from the selected point
		
		//ImagePlus img_traced_output;
		//img_traced_output 	= ImageConversions.ImagePlusToRGB(traced_img);
		
		//boolean first 					= true;
		boolean manual_start			= true;
		Hypothesis	take_hypothesis 	= new Hypothesis();
		Hypothesis  estimated_hyp 		= new Hypothesis();
		Hypothesis[] new_hypotheses 	= null;
		TinyBranchTrace current_branch 	= new TinyBranchTrace();
		ImagePlus sphere_img 			= new ImagePlus();
		
		double[] start_xyz = new double[]{x, y, z};//refineStartPoint(new double[]{x, y, z}, 20); // refine 
		
		while((manual_start || branch_queue.size()>0) && traced_branch_nr<max_branch_nr){
			
			if(manual_start){
				
				System.out.print("{");
				// detect best matching hypothesis at refined start point
				take_hypothesis = (new TinyBranchTrace(traced_img)).hyp_at(start_xyz);
				current_branch = getCurrentlyTracedBranch();
				current_branch.setStart(take_hypothesis);
				branches_mother_idxs[traced_branch_nr] = -1;
				
				//*** detection
				double scale = TinyBranch.check_bifurcations*current_branch.getCurrentTraceHypothesisRadius();
				Sphere sp = new Sphere(current_branch.getCurrentTraceHypothesisCenterpoint(), scale);
				sphere_img     = sp.extract(traced_img, TinyBranch.extract_sphere_resolution); 
				current_branch.bifurcation_detection_MeanShift3D(sphere_img, sp, manual_start); 
				// will process new hypotheses
				manual_start = false;
				System.out.print("}");
				
			}
			else{ // not a manual start
				
				System.out.print("\nbranch #"+traced_branch_nr+"/"+max_branch_nr+"<");
				//take first hypothesis if not the first execution
				take_hypothesis = branch_queue.get(0);	
				branch_queue.remove(0);
				current_branch = getCurrentlyTracedBranch();
				
				current_branch.setStart(take_hypothesis); 
				branches_mother_idxs[traced_branch_nr] = branch_queue_mother_idx.get(0);
				branch_queue_mother_idx.remove(0);
				
				// trace from the point on
				new_hypotheses = null;
				boolean enable_bif_detection = false;
				int loop_count = 0;
				while(loop_count<TinyBranch.N-1){ 
					
					current_branch.predict();
					current_branch.calculatePriors();
					current_branch.calculateLikelihoods();
					
					double sum_posts = current_branch.calculatePosteriors();
					
					System.out.print("=");
					
					if(sum_posts>0){
						
						estimated_hyp 	= current_branch.getMaximumPosteriorHypothesis();
						
						current_branch.storeHypothesis(estimated_hyp, false);
						
						//*** detection
						double scale = TinyBranch.check_bifurcations*current_branch.getCurrentTraceHypothesisRadius();
						double dx = current_branch.getCurrentPosX()-current_branch.getSeedPosX();
						double dy = current_branch.getCurrentPosY()-current_branch.getSeedPosY();
						double dz = current_branch.getCurrentPosZ()-current_branch.getSeedPosZ();
						
						if((dx*dx+dy*dy+dz*dz>scale*scale) && !enable_bif_detection)	
							enable_bif_detection = true; // || loop_count>3
						
						if(enable_bif_detection){
							Sphere sp 		= new Sphere(current_branch.getCurrentTraceHypothesisCenterpoint(), scale);
							sphere_img     	= sp.extract(traced_img, TinyBranch.extract_sphere_resolution); 
							System.out.print("*");
							
//							System.out.println("bif detection at ");
//							ArrayHandling.print1DArray(current_branch.getCurrentTraceHypothesisCenterpoint());
							
							boolean stop 	= 
									current_branch.bifurcation_detection_MeanShift3D(sphere_img, sp, manual_start); 
							// will process new hypotheses
							if(stop) break;
						}
//						else{
//							//Hypothesis.drawOverColorImage(estimated_hyp, img_traced_output, 255, 0, 0);
//						}
//						System.out.println(
//								"mean-shift detection says "+((stop)?"STOP":"CONTINUE"));
//						if(stop){
//							new FileSaver(sphere_img).saveAsTiffStack(String.format("branch_%d,mother_%d,loop_%d.tif", traced_branch_nr, branches_mother_idxs[traced_branch_nr], loop_count));
//						}
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
			
			new_hypotheses = current_branch.getNewSeedHypotheses();
			
			if(new_hypotheses!=null){ 
				int how_many_new_hyps = 0;
				for (int b = 0; b < new_hypotheses.length; b++) {
					if(new_hypotheses[b]!=null){
						//System.out.println("adding:");
						//Hypothesis.drawOverColorImage(new_hypotheses[b], img_traced_output, 255, 0, 0);
						//new_hypotheses[b].print();
						branch_queue.add(new_hypotheses[b]);
						branch_queue_mother_idx.add(traced_branch_nr);
						how_many_new_hyps++;
					}
				}
				
				// current queue
				System.out.println(how_many_new_hyps+" new branches");
			}
			
			traced_branch_nr++; // count next one once this is finished
			//if(first) first = false;
			
		}
		
		System.out.println("\n reconstructed "+traced_branch_nr+" branches");
		System.out.println("\n"+branch_queue.size()+" branches waiting in the queue");
		
		// debug
//		String name = String.format("seeds.tif");
//		img_traced_output.setTitle(name);
//		FileSaver fs = new FileSaver(img_traced_output);
//		fs.saveAsTiffStack(name);
		
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