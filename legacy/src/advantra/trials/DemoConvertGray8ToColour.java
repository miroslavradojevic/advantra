package advantra.trials;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.PlugIn;
import ij.process.ColorProcessor;

public class DemoConvertGray8ToColour implements PlugIn {

	public void run(String arg0) {
		
		IJ.run("Close All");
        IJ.run("Image Sequence..."); // open image sequence
        
        ImagePlus img_in 			= new ImagePlus("StackToExamine", IJ.getImage().getStack());
        
        int img_height 				= img_in.getStack().getHeight();
		int img_width 				= img_in.getStack().getWidth();
		int img_len     			= img_in.getStack().getSize();
		
		ImageStack im_stack 		= new ImageStack(img_width, img_height); // container for each layer
		
		for (int i = 1; i <= img_len; i++) {
			im_stack.addSlice(new ColorProcessor(img_in.getStack().getProcessor(i).createImage()));
			//System.out.println("stack size: "+im_stack.getHeight()+" x "+im_stack.getWidth()+" x "+im_stack.getSize());
		}
		
		ImagePlus img_col 			= new ImagePlus("ConvertedImage", im_stack);
		
		img_col.show();
		
	}
	
}
