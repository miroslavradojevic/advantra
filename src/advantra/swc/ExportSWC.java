package advantra.swc;

import java.io.*;
import java.util.ArrayList;

import advantra.trace.Hypothesis;
import advantra.trace.SeedPoint3D;
import advantra.trace.SeedPointQueue;

public class ExportSWC {
	
	String  			file_name;
	FileWriter 			fw; 
	
	SeedPointQueue		queue_seed_pts;			
	ArrayList<Integer>	queue_indexes;
	
	int					curr_entry_index;
	int					prev_entry_index;
	
	Hypothesis			next_hypothesis;
	
	boolean 			started;
	boolean 			queue_empty;
	
	public ExportSWC(String file_name){
		this.file_name 						= file_name;
		
		curr_entry_index 					= 	1;
		prev_entry_index 					= 	-1;
		
		queue_indexes	= new ArrayList<Integer>();
		queue_seed_pts	= new SeedPointQueue();
		
		started			= false;
		queue_empty		= false;
		
		next_hypothesis	= new Hypothesis();
		
		try{fw = new FileWriter(file_name);} 
		catch(IOException exIO){System.out.printf("SWCExport:SWCExport():Couldn't open SWC file for writing...");}
		
	}

	public void writeEntryBODY( 
			int branch_type, 
			Hypothesis input_hypothesis){
		if(!started){
				// add the point to the log
				try{fw.write(
				curr_entry_index+" "+
				branch_type+" "+
				input_hypothesis.getPositionY()	+" "+//I switched Y and Z so that it can be visualized in vaa3d
				input_hypothesis.getPositionX()	+" "+
				input_hypothesis.getPositionZ()	+" "+
				input_hypothesis.getNeuriteRadius()	+" "+
				prev_entry_index+"\n");} 
				catch(IOException exIO){System.out.println("Couldn't write to SWC.");}
				// append to queue
				queue_indexes.add(-1);
				queue_seed_pts.appendElement(new SeedPoint3D(
						input_hypothesis.getPositionX(),
						input_hypothesis.getPositionY(),
						input_hypothesis.getPositionZ(),
						-input_hypothesis.getOrientationX(),
						-input_hypothesis.getOrientationY(),
						-input_hypothesis.getOrientationZ(),
						input_hypothesis.getNeuriteRadius()
						));
				
				prev_entry_index=curr_entry_index;
				curr_entry_index++;
			}else{
				// add the point to the log
				try{fw.write(
				curr_entry_index+" "+
				branch_type+" "+
				input_hypothesis.getPositionY()	+" "+//I switched Y and Z so that it can be visualized in vaa3d
				input_hypothesis.getPositionX()	+" "+
				input_hypothesis.getPositionZ()	+" "+
				input_hypothesis.getNeuriteRadius()	+" "+
				prev_entry_index+"\n");} 
				catch(IOException exIO){System.out.println("Couldn't write to SWC.");}
				prev_entry_index=curr_entry_index;
				curr_entry_index++;
				
				
			}
			started = true;
			
		}
		
	public void writeEntryBRANCHING(int branch_type, 
				Hypothesis input_hypothesis,
				Hypothesis branch_hypothesis){
			
			// add the point to the log
			try{fw.write(
			curr_entry_index+" "+
			branch_type+" "+
			input_hypothesis.getPositionY()	+" "+ //I switched Y and Z so that it can be visualized in vaa3d
			input_hypothesis.getPositionX()	+" "+
			input_hypothesis.getPositionZ()	+" "+
			input_hypothesis.getNeuriteRadius()	+" "+
			prev_entry_index+"\n");} 
			catch(IOException exIO){System.out.println("Couldn't write to SWC.");}
			// append to queue
			queue_indexes.add(curr_entry_index);
			queue_seed_pts.appendElement(new SeedPoint3D(
					branch_hypothesis.getPositionX(),
					branch_hypothesis.getPositionY(),
					branch_hypothesis.getPositionZ(),
					branch_hypothesis.getOrientationX(),
					branch_hypothesis.getOrientationY(),
					branch_hypothesis.getOrientationZ(),
					branch_hypothesis.getNeuriteRadius()
					));
			
			prev_entry_index=curr_entry_index;
			curr_entry_index++;
			
			started = true;
			
	}
		
	public Hypothesis writeEntryEND(int branch_type, 
				Hypothesis input_hypothesis){
			// add the point to the log
			try{fw.write(
			curr_entry_index+" "+
			branch_type+" "+
			input_hypothesis.getPositionY()	+" "+ //I switched Y and Z so that it can be visualized in vaa3d
			input_hypothesis.getPositionX()	+" "+
			input_hypothesis.getPositionZ()	+" "+
			input_hypothesis.getNeuriteRadius()	+" "+
			prev_entry_index+"\n");} 
			catch(IOException exIO){System.out.println("Couldn't write to SWC.");}
			// take the next one from the queue
			if(queue_seed_pts.checkSize()<=0){
				System.out.println("no seed points to trace any more, empty queue.");
				queue_empty = true;
			}else{
				System.out.println("back to the point where started: prev="+prev_entry_index);
				prev_entry_index 	= takeSeedPointIndexFromQueue();
				curr_entry_index++;
				next_hypothesis		= takeSeedPointFromQueue();
			}
			
			started = true;
			
			return next_hypothesis;
	}
	
	public boolean isQueueEmpty(){
		return queue_empty;
	}
	
	public Hypothesis getNextSeedPointHypothesis(){
		Hypothesis hyp_out = new Hypothesis();
		hyp_out.setHypothesis(next_hypothesis);
		return hyp_out;
	}
	
	public Hypothesis takeSeedPointFromQueue(){
		Hypothesis hyp = new Hypothesis();
		SeedPoint3D p  = this.queue_seed_pts.takeElement();
		hyp.setHypothesis(
				p.getSeedPointPositionX(), 
				p.getSeedPointPositionY(), 
				p.getSeedPointPositionZ(), 
				p.getSeedPointOrientationX(), 
				p.getSeedPointOrientationY(), 
				p.getSeedPointOrientationZ(), 
				p.getRadius(), 2.0);
		return hyp;
	}
	
	public int takeSeedPointIndexFromQueue(){
		if(queue_indexes.size()<=0){
			System.err.println("SWCExport:takeSeedPointIndexFromQueue(): couldn't take the index from the list - it is empty!");
			System.exit(1);
		}
		int index_out = queue_indexes.remove(queue_indexes.size()-1);
		return index_out;

	}	
	/*
	public void writeBranchEntry(int branch_type, double x, double y, double z, double R){
		countEntries += 1;
		curr_entry = countEntries;
		prev_entry = parent_list.get(parent_list.size()-1);

		
		
		parent_list.remove(parent_list.size()-1);//remove the most recently added one
		

		parent_list.add(new Integer(curr_entry));
		
		if(EXPORT_DATA_TO_SWC_FILE){
			try{fw.write(countEntries+" "+branch_type+" "+x+" "+y+" "+z+" "+R+" "+prev_entry+"\n");} 
			catch(IOException exIO){System.out.printf("SWCExport:SWCExport():  \nCouldn't write to SWC.");}
		}
	}
	
	public void writeBranchEntryEnd(int T, double x, double y, double z, double R){
		
		branch_ended = true;
		countEntries+=1;
		curr_entry = countEntries;
		prev_entry = parent_list.get(0);
		parent_list.remove(parent_list.size()-1);//remove the most recently added one
		if(EXPORT_DATA_TO_SWC_FILE){
			try{fw.write(countEntries+" "+T+" "+x+" "+y+" "+z+" "+R+" "+prev_entry+"\n");} catch(IOException exIO){System.out.println("Couldn't write to SWC.");}
		}
	}	
	
	public void printCurrentParentList(){
		System.out.println("PARENT LIST["+parent_list.size()+"] for "+this.FILE_NAME+"    :   ");
		for (int i = 0; i < parent_list.size(); i++) {
			System.out.format("Parent[%02d]=%03d   ", i, parent_list.get(i));
			if(i%5==4){
				System.out.format("\n");
			}
		}
	}
	*/
	
	public void writelnSWC(String line){
		
			try{fw.write(line+"\n");} catch(IOException exIO){System.out.println("Couldn't write to SWC.");}
		
	}	

	public void closeSWC(){
		//if(EXPORT_DATA_TO_SWC_FILE){
			try{fw.close();} catch(IOException exIO){System.out.println("Couldn't write to SWC.");}
		//}
	}	
	
}
