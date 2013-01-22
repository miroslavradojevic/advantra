package advantra.trace;

import java.util.ArrayList;

import advantra.general.ArrayHandling;

public class SeedPointQueue {
	
	// list of seed points pending for trace
	// seed point = (point & orientation & radius)
	ArrayList<SeedPoint3D> 		queue_list;
	
	public SeedPointQueue(){
		
		queue_list = new ArrayList<SeedPoint3D>();
		
	}
	
	public int checkSize(){
		return queue_list.size();
	}
	
	public SeedPoint3D takeElement(){
		if(queue_list.size()<=0){
			System.err.println("SeedPointQueue:takeElement(): couldn't take the element from the list - it is empty!");
			System.exit(1);
		}
		return queue_list.remove(queue_list.size()-1);
	}
	
	public void appendElement(SeedPoint3D e){
		queue_list.add(e);// adds at the end of the list
	}
	
	public void printElements(){
		System.out.format("LIST SEED POINT QUEUE: \n");
		if(queue_list.size()<=0){
			System.out.format("***NONE***\n");
		}
		for (int i = 0; i < queue_list.size(); i++) {
			System.out.format("SEED_POINT[%4d]=  \n", i);
			ArrayHandling.print1DArray(queue_list.get(i).getSeedPointPosition());
			ArrayHandling.print1DArray(queue_list.get(i).getSeedPointOrientation());
			System.out.format("---------------------\n");
		}
	}

}
