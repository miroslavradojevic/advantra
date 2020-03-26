package advantra.plugins;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.NewImage;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

public class StackMaxZ implements PlugInFilter {

	protected ImageStack stack;
	
	public int setup(String arg, ImagePlus imp) {
		if(imp==null){
			IJ.noImage();// Displays a "no images are open" dialog box.
			return DONE;
		}
		stack = imp.getStack();
		return DOES_8G+NO_CHANGES;
	}
	
	public void run(ImageProcessor ip) {
	
		byte[] pixels; // to hold the slice pix

		int stack_width 	= stack.getWidth();
		int stack_height 	= stack.getHeight();
		int dimension 		= stack_width * stack_height;
		int[] max = new int[dimension];
		
		ImagePlus img_maxZ 	= NewImage.createByteImage("Maximum along Z", stack_width, stack_height, 1, NewImage.FILL_BLACK);
		byte[] maxZ	= (byte[]) img_maxZ.getProcessor().getPixels();
		
		// set the first layer as max at the beginning
		pixels = (byte[]) stack.getPixels(1);
		for (int i = 0; i < dimension; i++) {
			max[i] = pixels[i] & 0xff;
		}
		
		for (int i=2;i<=stack.getSize();i++) {
			pixels = (byte[]) stack.getPixels(i);
			for (int j=0;j<dimension;j++) {
				if((pixels[j] & 0xff) > max[j]){
					max[j] = pixels[j] & 0xff;
				}
			}
		}
		
		for (int j=0;j<dimension;j++) {
			maxZ[j] = (byte) (max[j] & 0xff);
		}
		
		img_maxZ.show();
		img_maxZ.updateAndDraw();
		
	}

}