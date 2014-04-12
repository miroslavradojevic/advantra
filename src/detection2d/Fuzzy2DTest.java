package detection2d;

import detection2d.Fuzzy2D;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

import java.util.Arrays;

/**
 * Created by miroslav on 3/25/14.
 * terminal call, using imagej jar
 * java -Xmx4096m -jar ~/jarlib/ij.jar -ijpath ~/ImageJ/plugins/  -run "Fuzzy2DTest"
 */
public class Fuzzy2DTest implements PlugIn {

	public void run(String s) {demo();}

	public static void demo()
	{

		float ncc_high_mean                  = (float)   Prefs.get("critpoint.detection2d.ncc_high_mean", 1f);
		float ncc_high_sigma 				= (float)   Prefs.get("critpoint.detection2d.ncc_high_sigma", 0.15f);
		float ncc_low_mean					= (float) 	Prefs.get("critpoint.detection2d.ncc_low_mean", 0.7f);
		float ncc_low_sigma					= (float) 	Prefs.get("critpoint.detection2d.ncc_low_sigma", 0.15f);

		float likelihood_high_mean           = (float)   Prefs.get("critpoint.detection2d.likelihood_high_mean", 0.9f);
		float likelihood_high_sigma 			= (float)   Prefs.get("critpoint.detection2d.likelihood_high_sigma", 0.2f);
		float likelihood_low_mean			= (float) 	Prefs.get("critpoint.detection2d.likelihood_low_mean", 0.0f);
		float likelihood_low_sigma			= (float) 	Prefs.get("critpoint.detection2d.likelihood_low_sigma", 0.2f);

		float output_sigma					= (float) 	Prefs.get("critpoint.detection2d.output_sigma", 0.4f);

		GenericDialog gd = new GenericDialog("TEst");
		gd.addNumericField("NCC HIGH", 	        ncc_high_mean, 			1,  5, "MEAN");
		gd.addNumericField("NCC HIGH",        	ncc_high_sigma, 	    1,  5, "SIGMA");
		gd.addNumericField("NCC LOW",           ncc_low_mean, 		    1,  5, "MEAN");
		gd.addNumericField("NCC LOW",           ncc_low_sigma, 		    1,  5, "SIGMA");

		gd.addNumericField("LIKELIHOOD HIGH", 	likelihood_high_mean, 			1,  10, "MEAN");
		gd.addNumericField("LIKELIHOOD HIGH",   likelihood_high_sigma, 	    	1,  10, "SIGMA");
		gd.addNumericField("LIKELIHOOD LOW",    likelihood_low_mean, 		    1,  10, "MEAN");
		gd.addNumericField("LIKELIHOOD LOW",    likelihood_low_sigma, 		    1,  10, "SIGMA");

		gd.addNumericField("OUT ",    output_sigma, 		    1,  10, "SIGMA");

		gd.showDialog();
		if (gd.wasCanceled()) return;

		ncc_high_mean   = (float) gd.getNextNumber();   	Prefs.set("critpoint.detection2d.ncc_high_mean", ncc_high_mean);
		ncc_high_sigma  = (float) gd.getNextNumber();  		Prefs.set("critpoint.detection2d.ncc_high_sigma", ncc_high_sigma);
		ncc_low_mean    = (float) gd.getNextNumber();  		Prefs.set("critpoint.detection2d.ncc_low_mean", ncc_low_mean);
		ncc_low_sigma	= (float) gd.getNextNumber(); 		Prefs.set("critpoint.detection2d.ncc_low_sigma", ncc_low_sigma);
		likelihood_high_mean  	= (float) gd.getNextNumber(); Prefs.set("critpoint.detection2d.likelihood_high_mean", 	likelihood_high_mean);
		likelihood_high_sigma 	= (float) gd.getNextNumber(); Prefs.set("critpoint.detection2d.likelihood_high_sigma", 	likelihood_high_sigma);
		likelihood_low_mean 	= (float) gd.getNextNumber(); Prefs.set("critpoint.detection2d.likelihood_low_mean", 	likelihood_low_mean);
		likelihood_low_sigma 	= (float) gd.getNextNumber(); Prefs.set("critpoint.detection2d.likelihood_low_sigma", 	likelihood_low_sigma);

		output_sigma 	= (float) gd.getNextNumber(); Prefs.set("critpoint.detection2d.output_sigma", 	output_sigma);



		Fuzzy2D test = new Fuzzy2D(
													  // ncc
													  ncc_high_mean, ncc_high_sigma,      // high
													  ncc_low_mean, ncc_low_sigma,   // low
													  // lhood
													  likelihood_high_mean, likelihood_high_sigma, // high
													  likelihood_low_mean, likelihood_low_sigma,  // low
													  output_sigma  // std output membership functions - defines separation
		);
		test.showFuzzification();
		test.showDefuzzification();
		test.showDefuzzificationSurface();

		float[] temp = new float[2];
		System.out.print("(1, 1) -> " + test.branchStrengthDefuzzified(1f, 1f));
		test.branchStrengthFuzzified(1f, 1f, temp);
		System.out.println(" -> " + Arrays.toString(temp));
		test.showAgg();

		System.out.print("(.9, .9) -> " + test.branchStrengthDefuzzified(.9f, .9f));
		test.branchStrengthFuzzified(.9f, .9f, temp);
		System.out.println(" -> " + Arrays.toString(temp));
		test.showAgg();

        System.out.print("endpoint on membership (.9, .9, .5, .3) -> " + test.endpointFuzzified(.9f, .9f, .5f, .3f));

	}

}
