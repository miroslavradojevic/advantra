package advantra.tools;

import java.util.ArrayList;

import ij.ImagePlus;

public class ConnetedRegsSegment {
	
	/*
	 * this class takes gray8 image (2d or 3d), 
	 * and labels connected regions that have exactly the same intensity
	 * wit hthe same label
	 * therefore it's convenient for images with just several 
	 * grey levels, mostly interesting for binary images to
	 * determine the connectivity of the landmarks
	 */
	
	ImagePlus img;
	boolean   imgIsStack;
	
	public ConnetedRegsSegment(ImagePlus img){
		
		this.img = img;
		imgIsStack = img.getStack().getSize()>1;
		ArrayList<int[]> seed_queue = new ArrayList<int[]>(); 
		seed_queue.add(new int[]{0, 0, 0}); // initial seed
		
	}
	
	public void run(){
		
		int H = img.getHeight();
		int W = img.getWidth();
		int L = img.getStack().getSize();
		
		byte[][] img_array = new byte[L][H*W];
		byte[][] img_labels = new byte[L][H*W];
		
		int total_labelled = 0;
		//int[] seed = new int[]{0, 0, 0};
		
		while(total_labelled<H*W*L){
			
		}
		
	}
	
	private int getNeighbours(int x, int y, int z, int[][] locations){
		
		//int[][] locations = new int[26][3];
		int nr_neighbours = 0;
		
		if(imgIsStack){
			// 26 neigh.
			for (int i = x-1; i <= x+1; i++) {
				for (int j = y-1; j <= y+1; j++) {
					for (int k = z-1; k < z+1; k++) {
						if(!(i==x && j==y && k==z) && i>=0 && j>=0 && k>=0 && i<img.getHeight() && j<img.getWidth() && k<img.getStack().getSize()){
							
							locations[nr_neighbours][0] = i;
							locations[nr_neighbours][1] = j;
							locations[nr_neighbours][2] = k;
							nr_neighbours++;
							
						}
					}
				}
			}
		}
		else{
			// 8 neigh.
			for (int i = x-1; i <= x+1; i++) {
				for (int j = y-1; j <= y+1; j++) {
					//for (int k = z-1; k < z+1; k++) {
					if(!(i==x && j==y) && i>=0 && j>=0 && i<img.getHeight() && j<img.getWidth()){
							
						locations[nr_neighbours][0] = i;
						locations[nr_neighbours][1] = j;
						//locations[nr_neighbours][2] = k;
						nr_neighbours++;
							
					}
					//}
				}
			}
		}
		
		return nr_neighbours; // nr. of extracted neighbours
		
	}

}
