package detection2d;

import detection2d.Fuzzy2D;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.NonBlockingGenericDialog;
import ij.plugin.PlugIn;

import java.util.Arrays;

/**
 * Created by miroslav on 3/25/14.
 * terminal call
 * java -jar ~/jarlib/ij.jar -ijpath ~/ImageJ/plugins/  -run "Fuzzy2DTest"
 */
public class Fuzzy2DTest implements PlugIn {

	NonBlockingGenericDialog gd = new NonBlockingGenericDialog("TEST");

	public void run(String s) {demo();}

	public void demo()
	{

		float ncc_high                  = (float)   Prefs.get("critpoint.detection2d.ncc_high", 0.9f);
		float ncc_low					= (float) 	Prefs.get("critpoint.detection2d.ncc_low", 0.6f);

		float likelihood_high           = (float)   Prefs.get("critpoint.detection2d.likelihood_high", 0.6f);
		float likelihood_low			= (float) 	Prefs.get("critpoint.detection2d.likelihood_low", 0.2f);

        float smoothness_high           = (float)   Prefs.get("critpoint.detection2d.smoothness_high", 10f);
        float smoothness_low			= (float) 	Prefs.get("critpoint.detection2d.smoothness_low", 15f);

		float output_sigma					= (float) 	Prefs.get("critpoint.detection2d.output_sigma", 0.4f);

//		NonBlockingGenericDialog gd1 = new NonBlockingGenericDialog("my_title");
//		gd1.addMessage("Hehehhe");
//		gd1.showDialog();
//		if (gd1.wasCanceled()) return;

		gd.addNumericField("NCC HIGH", 	        ncc_high, 			        2,  5, "");
		gd.addNumericField("NCC LOW",           ncc_low, 		            2,  5, "");
		gd.addNumericField("LIKELIHOOD HIGH", 	likelihood_high, 			2,  10, "");
		gd.addNumericField("LIKELIHOOD LOW",    likelihood_low, 		    2,  10, "");
        gd.addNumericField("SMOOTHNESS HIGH", 	smoothness_high, 			2,  10, "");
        gd.addNumericField("SMOOTHNESS LOW",    smoothness_low, 		    2,  10, "");

		gd.addNumericField("OUT SIGMA",    		output_sigma, 		    	1,  10, "");

		gd.showDialog();
		if (gd.wasCanceled()) return;

		ncc_high   = (float) gd.getNextNumber();   	        Prefs.set("critpoint.detection2d.ncc_high",         ncc_high);
		ncc_low    = (float) gd.getNextNumber();  		    Prefs.set("critpoint.detection2d.ncc_low",          ncc_low);
		likelihood_high  	= (float) gd.getNextNumber();   Prefs.set("critpoint.detection2d.likelihood_high", 	likelihood_high);
		likelihood_low 	= (float) gd.getNextNumber();       Prefs.set("critpoint.detection2d.likelihood_low", 	likelihood_low);
        smoothness_high  	= (float) gd.getNextNumber();   Prefs.set("critpoint.detection2d.smoothness_high", 	smoothness_high);
        smoothness_low 	= (float) gd.getNextNumber();       Prefs.set("critpoint.detection2d.smoothness_low", 	smoothness_low);

		output_sigma 	= (float) gd.getNextNumber(); 		Prefs.set("critpoint.detection2d.output_sigma", 	output_sigma);



		Fuzzy2D test = new Fuzzy2D(
													  // ncc
													  ncc_high,
													  ncc_low,
													  // lhood
													  likelihood_high,
													  likelihood_low,
                                                      // smoothness
                                                      smoothness_high,
                                                      smoothness_low,
													  output_sigma  // std output membership functions - defines separation
		);
		test.showFuzzification();
		test.showDefuzzification();

		// test when having one input
		float ncc1 = .9f;
		float lhood1 = .9f;
		float smooth1 = 5f;
		float[] out = new float[3];    // critpoint
		float[] out1 = new float[4];   // on/off/none for each branch
		test.verbose = true;
		test.critpointScore(ncc1, lhood1, smooth1, out, out1);
		test.clearLog();
		test.critpointScore(
                ncc1, lhood1, smooth1,
                ncc1, lhood1, smooth1,
                ncc1, lhood1, smooth1,
                out, out1);
		new ImagePlus("", test.fls_steps).show();
		System.out.println(Arrays.toString(out));

	}

}
