package advantra.commands;

import java.io.File;

import ij.io.FileSaver;
import advantra.general.CreateDirectory;
import advantra.tools.MeanShiftMasks;

public class DemoMsMasks {
	
	public static void main(String[] args){
		
		double ang_span = Math.PI/4; // can range from small arc: ~0 till full angular span 4*pi steradians
		int resolution 	= 32;
		
		MeanShiftMasks ms_masks = new MeanShiftMasks(resolution, ang_span);
		ms_masks.createMasks();
		
		String output_directory = 
				System.getProperty("user.dir")+File.separator+
				String.format("mask_angularDist_%f", ang_span)+File.separator;
		
		CreateDirectory.createOneDir(output_directory);
		
		for (int i = 0; i < resolution; i++) {
			for (int j = 0; j < 2*resolution; j++) {
				
				String file_name = String.format("mask_at%d_%d.tif", i, j);
				(new FileSaver(ms_masks.getMask(i, j))).saveAsTiff(String.format(output_directory+file_name));
				
			}
		}
		
	}

}
