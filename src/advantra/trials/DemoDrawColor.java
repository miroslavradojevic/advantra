package advantra.trials;

import advantra.general.ImageConversions;
import advantra.shapes.Sphere;
import ij.IJ;
import ij.ImagePlus;
import ij.plugin.PlugIn;

public class DemoDrawColor  implements PlugIn  {

	public void run(String arg0) {
		
		IJ.run("Close All");
		
		IJ.run("Image Sequence...", "open=[] number=[] starting=1 increment=1 scale=[] file=[] or=[] convert_to_rgb sort");
        
		ImagePlus img_in 		= new ImagePlus("LoadedSequence", IJ.getImage().getStack());
		
		int colour_intensity = 180; 
		// coords where to draw
		/*
		for (int i = 0; i < img_in.getWidth(); i++) {
			for (int j = 0; j < img_in.getHeight(); j++) {
				for (int j2 = 0; j2 < img_in.getStack().getSize(); j2++) {
					if(i==j){
						ImageConversions.setJet256ColorValue(img_in, i, j, j2, colour_intensity);
					}
				}
			}
		}
		*/
		int radius = 4;
		int[][] roi_coord = new int[3][Sphere.numberOfVoxInSphere(radius)];
		int nr = (new Sphere(8, 3, 5, radius)).extractCoords(img_in, roi_coord);
		
		int[][] to_extract = new int[3][nr];
		for (int i = 0; i < nr; i++) {
			to_extract[0][i] = roi_coord[0][i];
			to_extract[1][i] = roi_coord[1][i];
			to_extract[2][i] = roi_coord[2][i];
		}
		ImageConversions.setJet256ColorValues(img_in, to_extract[0], to_extract[1], to_extract[2], colour_intensity);
		
		img_in.show();
					
	}
	
}
