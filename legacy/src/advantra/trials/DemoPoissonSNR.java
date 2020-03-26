package advantra.trials;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.NewImage;
import ij.plugin.PlugIn;
import imagescience.random.PoissonGenerator;
import imagescience.utility.Progressor;

public class DemoPoissonSNR  implements PlugIn {

	ImagePlus im;
	boolean doPoisson;
	byte[] pixels; // to store the values
	
	public void run(String arg0) {
		
		System.out.println("DEMO: add Poisson noise with set SNR...");
		
		IJ.run("Close All");
		
		doPoisson = true;
		
		//// dialog window ////
		String menu_title = "Demo Poisson: SNR...";
		byte bg = 10;
		double SNR = 5;
		
		GenericDialog gd = new GenericDialog(menu_title, IJ.getInstance());
		
		gd.addMessage("Background ");
		gd.addNumericField( "intensity level :", bg, 3, 5, "brightness" );

		gd.addMessage("SNR");
		gd.addNumericField( "ratio :", SNR, 2, 5, "mean/variance" );// label, defaultVal, digits, columns, units
		
		gd.showDialog();
		
		if (gd.wasCanceled()) {
			return;
		}
		
		bg 			= (byte)gd.getNextNumber();
		SNR			= (double)gd.getNextNumber();

		
		// define the signal level alphaSNR using equation alphaSNR/sqrt(alphaSNR+bg)=snr
		byte alphaSNR = (byte)((Math.pow(SNR, 2) + SNR * Math.sqrt(Math.pow(SNR, 2) + 4 * bg)) / 2.0 ); 
		
  		ImagePlus output_img = NewImage.createByteImage("poisson_noise", 512, 512, 1, NewImage.FILL_BLACK);
  		// arguments: 									java.lang.String title, int width, int height, int slices, int options
		
		int rows = output_img.getHeight();
		int cols = output_img.getWidth();

		pixels = new byte[rows*cols];
		pixels = (byte[])output_img.getProcessor().getPixels(); // pixels refers to the pixels in image stack now
		
  		// add the background 
  		for (int i = 0; i < rows*cols; i++) {
  				
			pixels[i] = (byte)bg; // add the background level
		
  		}
  		
  		// create synthetic circle
  		for (int i = 0; i < pixels.length; i++) {
  			if( Math.pow(i/cols-rows/2, 2)+Math.pow(i%cols-cols/2, 2) <= Math.pow(Math.min(rows, cols)/4, 2) ){
  				// it belongs to the circle - set it to logic true -> alphaSNR
  				pixels[i] += alphaSNR;// add it to the background
  			}
  		}
  	  	
  		// add poisson noise, define noise generator
	  	PoissonGenerator poissonGen = new PoissonGenerator();
	  		
	  	System.out.println("Adding Poisson noise...");
	  		
		Progressor pgs = new Progressor(); // progress bar
		pgs.display(true);
		pgs.steps(rows*cols); 
		pgs.start();
			
	  		for (int i = 0; i < rows*cols; i++) {
				int value = pixels[i] & 0xff;// read the byte value from the image
			  	pixels[i] = (byte)poissonGen.next(value);
				pgs.step();
			}
	  		
	  	pgs.stop();
  		
  		output_img.show();
  		
	}
	
}
