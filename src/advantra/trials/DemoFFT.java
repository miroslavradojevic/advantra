package advantra.trials;

import advantra.filter.Lowpass;
import advantra.filter.Template;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.NewImage;
import ij.plugin.PlugIn;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;

import imagescience.fourier.FFT;
import imagescience.image.Axes;
import imagescience.image.FloatImage;

// this plug in is intended as a demo for imagescience FFT function usage
// fourier transform of the image is obtained and the image is filtered
// with defined lowpass filter

public class DemoFFT implements PlugIn {
	
	ImageConverter conv;
	ImagePlus im;
	ImageProcessor imp;

	public void run(String arg0) {
		
		System.out.println("run()...");
		IJ.run("Close All");
		im = IJ.openImage();
		
		im.show();
		
		// DEFINE WHICH DATA TYPE YOU WORK WITH - FLOAT
		imp = im.getProcessor().convertToFloat();
		
		int rows = im.getHeight();
		int cols = im.getWidth();
		
		FloatImage real  	= new FloatImage(new ImagePlus("REAL FFT", imp));
		FloatImage imag 	= new FloatImage(NewImage.createFloatImage("IMAG FFT", 
				cols, rows, 1, NewImage.FILL_BLACK));

		FFT fft = new FFT();
		fft.forward(real, imag, new Axes(true, true));
		
		// show it
		real.imageplus().show();
		imag.imageplus().show();
		
		String menu_title = "Demo FFT: lowpass filter...";
		float cutOff = 0.4f;
		float n = 10;
		GenericDialog gd = new GenericDialog(menu_title, IJ.getInstance());
		
		gd.addMessage("Threshold");
		gd.addNumericField( "cut-off :", cutOff, 2, 5, "0.0-0.5" );

		gd.addMessage("Slope degree");
		gd.addNumericField( "n :", n, 2, 5, "~10" );// label, defaultVal, digits, columns, units
		
		gd.showDialog();
		if (gd.wasCanceled()) {
			return;
		}
		
		cutOff 		= (float)gd.getNextNumber();
		n 			= (float)gd.getNextNumber();
		
		
		float[] radius  	= Template.radius(rows, cols);
		float[] filter_fq 	= Lowpass.createRow(radius, cutOff, 5);
		

		
		float[] real_fft 	= (float[])real.imageplus().getProcessor().getPixels();
		float[] imag_fft 	= (float[])imag.imageplus().getProcessor().getPixels();
		
		// filtering
		
		for (int i = 0; i < rows*cols; i++) {
			real_fft[i] = real_fft[i] * filter_fq[i];
			imag_fft[i] = imag_fft[i] * filter_fq[i];
		}
		
		// not necessary because they are already linked
		//real.imageplus().getProcessor().setPixels(real_fft);
		//imag.imageplus().getProcessor().setPixels(imag_fft);
		
		fft.inverse(real, imag, new Axes(true, true));
		
		for (int i = 0; i < filter_fq.length; i++) {
		System.out.print(filter_fq[i]+"   ");
		if(i%cols==cols-1){
			System.out.println("");
		}
	}
		
		// show the inverse transform
		ImagePlus im_out = real.imageplus();
		im_out.setTitle("AFTER FILTERING");
		im_out.show();
		
	}
}