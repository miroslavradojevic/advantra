package advantra.tools;

import java.io.File;

import ij.ImagePlus;

public class Test {
	
	public static void main(String[] args){
		System.out.println("Hello? is it working?");
		
		String file_name = "image_when_having_1_results.tif";
		file_name = (new File(file_name)).getAbsolutePath();
		
		ImagePlus img = new ImagePlus(file_name);
		
		Get_Connected_Regions try_it = new Get_Connected_Regions(img);
		System.out.println("class initialized");
		
		try_it.run();
	}

}
