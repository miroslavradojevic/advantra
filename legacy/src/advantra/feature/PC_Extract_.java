package advantra.feature;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

public class PC_Extract_ implements PlugIn {
	
	/*
	 * THIS PLUGIN WILL EXTRACT THE PHASE CONGRUENCY VALUES
	 * 
	 * ORIGINAL PAPER
	 * DOI: 10.1007/s004260000024
	 * 
	 * MATLAB IMPLEMENTATION AND FURTHER DETAILS AT:
	 * http://www.csse.uwa.edu.au/~pk/research/matlabfns/
	 * 
	 * INPUT IS B&W, BYTE IMAGE SELECTED BY USER
	 * FINAL OUTPUT IS PHASE CONGRUENCY FOE EACH PIXEL PER ORIENTATION 
	 * PC[ORIENTATION] WHERE PC ---> [0, 1]
	 * 
	 * DEPENDING ON THE PLUGIN CALL PC[ORIENT] CAN BE 
	 * 1. DISPLAYED, (IMAGEJ PLUGIN CALL)
	 * 2. SAVED TO THE WORKING FOLDER (JAVA TERMINAL CALL)
	 * 	java -cp  /home/miroslav/fiji/jars/ij.jar:/home/miroslav/fiji/jars/imagescience.jar:/home/miroslav/fiji/plugins/advantra.jars/flanagan.jar:/home/miroslav/fiji/plugins/Advantra_.jar advantra.feature.PC_Extract_ /home/miroslav/lena.png 4 12 3 2.1 0.55 10.0 0.5 10 -2 >> log.txt
	 * 
	 * OR 
	 * 3. JUST EXPORTED AS PRODUCT OF CLASS METHOD FOR FURTHER USAGE (method getPC of PC_Extract class)
	 */
		
		private static final String plugInTitle = "Phase Congruency v3.0";
		
		public void run(String arg){
			
	        // check if appropriate version of ImageJ is installed
	        if (IJ.versionLessThan("1.24t"))
	            return;

	        // check Java
	        if (System.getProperty("java.version").compareTo("1.6.0") < 0) {
	            IJ.error("This plugin has been developed and tested with Java, version 1.7 and higher.\n" + "Please upgrade your JVM.");
	            return;
	        }
	        
	        ImagePlus img;
	        ImageProcessor imp;
	        
	        IJ.run("Close All");
	        img = IJ.openImage();
	        img.show();
	        
			// DEFINE WHICH DATA TYPE YOU WORK WITH - FLOAT
			imp = img.getProcessor().convertToFloat();
			
			PC_Params params = new PC_Params(); // introduce class that will contain parameters
			
			// check opened windows with images
	        int[] wList = WindowManager.getIDList();
	        if (wList == null) {
	            IJ.noImage();
	            return;
	        }
	        
	        String[] titles = new String[wList.length];
	        for (int i = 0; i < wList.length; i++) {
	            ImagePlus opened_images = WindowManager.getImage(wList[i]);
	            if (opened_images != null)
	                titles[i] = opened_images.getTitle();
	            else
	                titles[i] = "";
	        }

			GenericDialog gd = new GenericDialog(plugInTitle, IJ.getInstance());
			
			gd.addMessage("Select image:");
			gd.addChoice("source image:", titles, titles[0]);
			
			// gd.addMessage("Number of wavelet scales: ");
			gd.addNumericField( "Number of wavelet scales: ", params.nscale, 0, 5, "" );

			// gd.addMessage("Number of filter orientations: ");
			gd.addNumericField( "Number of filter orientations: ", params.norient, 0, 5, "" );
			
			// gd.addMessage("Wavelength of smallest scale filter: ");
			gd.addNumericField( "Wavelength of smallest scale filter: ", params.minWvl, 0, 5, "px" );

			// gd.addMessage("Scaling factor between successive filters: ");
			gd.addNumericField( "Scaling factor between successive filters: ", params.mult, 2, 5, "" );

//			gd.addMessage("Ratio of the standard deviation of the Gaussian \n " +
//					"describing the log Gabor filter's transfer function \n" +
//					"in the frequency domain to the filter center frequency: ");
			
			gd.addNumericField( "std dev logGabor's Gaussian/center frequency :", params.sigmaOnf, 2, 5, "" );

			gd.addMessage("No of standard deviations of the noise energy beyond \n" +
					"the mean at which we set the noise threshold point: ");
			gd.addNumericField( "k :", params.k, 1, 5, "" );
			
			gd.addMessage("The fractional measure of frequency spread \n" +
					"below which phase congruency values get penalized:");
			gd.addNumericField( "cutOff :", params.cutOff, 2, 5, "" );
			
			gd.addMessage("Sharpness of the transition in the sigmoid function used to weight phase \n" +
					"congruency for frequency spread:");
			gd.addNumericField( "g :", params.g, 1, 5, "" );
			
			gd.addMessage(	"-1 use median of smallest scale filter responses \n" +
							"-2 use mode of smallest scale filter responses \n" +
							" 0+ use noiseMethod value as the fixed noise threshold");
			gd.addNumericField( "", params.noiseMethod, 0, 5, "" );
			
			gd.addCheckbox("save results as PNG", true);
			
			gd.showDialog();
			if (gd.wasCanceled()) {
				return;
			}
			
			// parameters 
			params.nscale 		= (int)gd.getNextNumber();
			params.norient 		= (int)gd.getNextNumber();
			params.minWvl 		= (int)gd.getNextNumber();
			params.mult 		= (float)gd.getNextNumber();
			params.sigmaOnf 	= (float)gd.getNextNumber();
			params.k 			= (float)gd.getNextNumber();
			params.cutOff 		= (float)gd.getNextNumber();
			params.g 			= (float)gd.getNextNumber();
			params.noiseMethod 	= (int)gd.getNextNumber();
			
			boolean saveIt		= gd.getNextBoolean();
			
			// check the type of the image (only 8bit grey-level)
	        try {
	            if ( !isGrayscale(img) || !isGrayscale(img) ) {
	                IJ.error("Input must be a 8-bit gray-scale image.\n");
	                return;
	            }

	            
	        } finally {
	            img.unlock();
	        }
	        
	        System.out.println("start extraction...");
	        System.out.println("save it: "+saveIt);
			PC_Extract.run(imp, img.getShortTitle(), params, true, saveIt);
			
			System.out.println("Extracted! ");
			
		}

		public static void main(String[] args){
			System.out.println("main()...");
	        ImagePlus img;
	        ImageProcessor imp;
			
	        // check Java
	        if (System.getProperty("java.version").compareTo("1.6.0") < 0) {
	            IJ.error("This plugin has been developed and tested with Java, version 1.7 and higher.\n" + "Please upgrade your JVM.");
	            return;
	        }
			
			if(args.length != 10){
				System.out.println("wrong number of arguments entered - unfortunately every argument has to be entered...");
			}
			
			if(args.length == 10){
				
				img = IJ.openImage(args[0]);
				
				// DEFINE WHICH DATA TYPE YOU WORK WITH - FLOAT
				imp = img.getProcessor().convertToFloat();

				PC_Params par 	= new PC_Params();
				
				par.nscale		= (int)Integer.parseInt(args[1]);
				par.norient		= (int)Integer.parseInt(args[2]);
				par.minWvl		= (int)Integer.parseInt(args[3]);

				par.mult		= (float)Float.parseFloat(args[4]);
				par.sigmaOnf	= (float)Float.parseFloat(args[5]);
				par.k			= (float)Float.parseFloat(args[6]);
				par.cutOff		= (float)Float.parseFloat(args[7]);
				par.g			= (float)Float.parseFloat(args[8]);
				
				par.noiseMethod	= (int)Integer.parseInt(args[9]);
				
				par.printAll();
				
				PC_Extract.run(imp, img.getShortTitle(), par, false, true);
				
				// this piece of code can be embedded in other packages
				// PC_Extract pc_ex = new PC_Extract(imp, img.getShortTitle(), par, false);
				// float[][] output = pc_ex.getPC();
				
			}

		}
		
		static public void showAbout() {
	        IJ.showMessage("About Phase Congruency...", "This plug-in has been written by M. Radojevic.\nThe code refers to phasecong3.m written by P. Kovesi.\n");
	    }
		
		private static boolean isGrayscale(ImagePlus image) {
	        return ((image.getType() == ImagePlus.GRAY8) ); //|| (image.getType() == ImagePlus.GRAY16) || (image.getType() == ImagePlus.GRAY32)
	    }
}
