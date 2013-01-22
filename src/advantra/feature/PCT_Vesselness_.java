package advantra.feature;

import java.io.File;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

public class PCT_Vesselness_  implements PlugIn  {

	/*
	 * THIS PLUGIN WILL EXTRACT THE PCT VESSELNESS VALUES
	 * 
	 * ORIGINAL PAPER
	 * DOI: 10.1109/TIP.2012.2185938
	 * 
	 * INPUT IS B&W, BYTE IMAGE SELECTED BY USER
	 * FINAL OUTPUT IS PCT VESSELNESS IMAGE 
	 * 
	 * DEPENDING ON THE PLUGIN CALL VESSELNESS CAN BE 
	 * 1. DISPLAYED, (IMAGEJ PLUGIN CALL)
	 * 2. SAVED TO THE WORKING FOLDER (JAVA TERMINAL CALL)
	 * java -cp  
	 * ij.jar:
	 * imagescience.jar:
	 * Jama-1.0.2.jar:
	 * Advantra_.jar 
	 * advantra.feature.PCT_Vesselness_ 
	 * lena.png 
	 * "0.5, 1, 2, 3" 8 10 0.5 1
	 * 
	 * OR 
	 * 3. JUST EXPORTED AS PRODUCT OF CLASS METHOD FOR FURTHER USAGE (method getPC of PC_Extract class)
	 */
	
	private static final String plugInTitle = "pct vesselness v0.1";
	
	public void run(String arg) {
		
        ImagePlus img;
        
        IJ.run("Close All");
        img = IJ.openImage();
        img.show();
        
		// DEFINE WHICH DATA TYPE YOU WORK WITH - FLOAT
        ImageProcessor imp = img.getProcessor().convertToFloat();
		
		// DEFINE THE PARAMETERS THAT WILL BE USED
		// pc 	params
		String scales 	= "1";// "0.5, 1, 2, 3, 4"; //comma separated
		int norients 	= 6;
		double k		= 10.0;
		// "pc vesselness" params
		double b 		= 0.5;
		double c		= 0.5;
		
		// GENERIC DIALOG
		GenericDialog gd = new GenericDialog(plugInTitle, IJ.getInstance());
		
		gd.addMessage("Define scales:");
		gd.addStringField("scales:", scales);
		
		gd.addMessage("Number of orientations: ");
		gd.addNumericField("orientations per pi :", norients, 2, 5, "");

		gd.addMessage("k parameter for the noise : ");
		gd.addNumericField("Nose energy standard deviations number:", k, 2, 5, "");
		
		gd.addMessage("Vesselness parameters : ");
		gd.addNumericField( "b :", 	b, 2, 5, "" );

		gd.addMessage("Vesselness parameters : ");
		gd.addNumericField("c :", 	c, 2, 5, "");
		
		gd.showDialog();
		if (gd.wasCanceled()) {
			return;
		}
		
		// parameters 
		scales				= gd.getNextString();
		norients 			= (int)gd.getNextNumber();
		k					= (double)gd.getNextNumber();
		b					= (double)gd.getNextNumber();
		c 					= (double)gd.getNextNumber();
		
		// check the type of the image (only 8bit grey-level)
        try {
        	
            if (!isGray8(img)) {
                IJ.error("Input must be a 8-bit gray-scale image.\n");
                return;
            }
            
        } finally {
        	
            img.unlock();
            
        }
        
        System.out.println("extract vesselness...");
        
        (new PCT_Vesselness(imp, img.getShortTitle())).run(scales, norients, k, b, c, true, true);
		
		System.out.println("extracted! ");
		
	}

	public static void main(String[] args){ // terminal call
		
		ImagePlus img;
        ImageProcessor imp;
		
		if(args.length != 6){
			
			System.out.println("enter arguments:");
			System.out.println(
					"1 - image path				\n" +
					"2 - string with scales... \"0.5, 1, 2, 3\" \n" +
					"3 - nr_orientations...  8	\n" +
					"4 - k (pc parameter)... 10	\n" +
					"5 - b... 0.5				\n" +
					"6 - c... 1  				\n");
			return;
			
		}
		
		// path
		img = IJ.openImage((new File(args[0])).getAbsolutePath());
		
		// DEFINE WHICH DATA TYPE YOU WORK WITH - FLOAT
		imp = img.getProcessor().convertToFloat();
		
		String scales			= args[1];
		int norients			= (int)Integer.parseInt(args[2]);
		double k				= (double)Integer.parseInt(args[3]);
		double b				= (double)Float.parseFloat(args[4]);
		double c				= (double)Float.parseFloat(args[5]);
			
		System.out.println("loaded params: ");
		System.out.println(scales);
		System.out.println(norients);
		System.out.println(k);
		System.out.println(b);
		System.out.println(c);
			
		(new PCT_Vesselness(imp, img.getShortTitle())).run(scales, norients, k, b, c, false, true); 
		// save it only
			
		
	}
	
	private static boolean isGray8(ImagePlus image){
		return (image.getType() == ImagePlus.GRAY8);
	}
	
}
