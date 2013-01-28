package advantra.plugins;

import java.io.File;

import advantra.tools.Find_Connected_Regions;

import ij.ImagePlus;

public class Test {
	
	public static void main(String[] args){
		
		
		String file_name = args[0];
		
		int[] p1 = new int[]{ Integer.parseInt(args[1]),Integer.parseInt(args[2]),Integer.parseInt(args[3]) };
		int[] p2 = new int[]{ Integer.parseInt(args[4]),Integer.parseInt(args[5]),Integer.parseInt(args[6]) };
		
		file_name = (new File(file_name)).getAbsolutePath();
		if(new File(file_name).exists()){
			System.out.println("loaded: "+file_name);
		}
		else{
			return;
		}
		
		ImagePlus img = new ImagePlus(file_name);
		
		Find_Connected_Regions try_it = new Find_Connected_Regions(img, true);
		try_it.run("");
		
//		int[][] locs = try_it.getLocations(1);
//		for (int i = 0; i < locs.length; i++) {
//			System.out.println(""+locs[i][0]+" , "+locs[i][1]+" , "+locs[i][2]);
//		}
		System.out.println("\n---\n");
		System.out.println(try_it.getNrConnectedRegions()+" connected regions");
//		ImagePlus to_save = try_it.showLabels();
//		(new FileSaver(to_save)).saveAsTiff();
		
		System.out.format("(%d.%d,%d) and (%d,%d,%d) are connected:  %b\n", p1[0],p1[1],p1[2],p2[0],p2[1],p2[2],try_it.connected(p1, p2));
		System.out.format("(%d.%d,%d) belongs to biggest reg:  %b\n", p1[0],p1[1],p1[2],try_it.belongsToBiggestRegion(p1));
		System.out.format("(%d.%d,%d) belongs to biggest reg:  %b\n", p2[0],p2[1],p2[2],try_it.belongsToBiggestRegion(p2));
		System.out.format("(%d.%d,%d) rank:  %d \n",  p1[0],p1[1],p1[2],try_it.getRank(p1));
		System.out.format("(%d.%d,%d) rank:  %d \n",  p2[0],p2[1],p2[2],try_it.getRank(p2));
		
	}

}
