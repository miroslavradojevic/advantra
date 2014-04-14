package generate;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.io.FileSaver;
import ij.plugin.PlugIn;

import java.io.File;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 10/23/13
 * Time: 4:49 PM
 * a wrapper plugin class that would use .swc files to create 2d or 3d
 * synthetic images of neurons with different snr's and ground truth
 * critical points inferred from the swc itself
 * input : path to swc, output: tif image stack or image, critical point swc, end point swc
 */
public class GeneratorSwcDemo implements PlugIn {

    /*
        terminal call :
        java -cp "$HOME/critpoint/*:$HOME/jarlib/*" generate.GeneratorSwcDemo
     */
    public static void main(String args[]){

        float   k;            	// scales the gaussian sigma - drop of the intensity as you go from the cone centerline, correlated to radius read from swc
        float   snr;            // noise ratio
        File    fileSWC;        // swc input file
        String  pathInSwc;      // path to input swc file
		boolean is2D;			// whether it will be generated in 2d

        /*
            check argument nr.
         */
        if (args.length!=4) {
            System.out.println("# GENERATE SYNTHETIC NEURON IMAGE FROM SWC #");
            System.out.println("USAGE:\t\tpath_to_SWC\tk\tSNR\tis2D");
            return;
        }

        /*
            check if the swc input file exists
         */
        pathInSwc = new File(args[0]).getAbsolutePath();
        fileSWC = new File(pathInSwc);

        if (!fileSWC.exists()) {
            System.out.println("file "+pathInSwc+" does not exist!");
            return;
        }

        k       = Float.valueOf(args[1]);
        snr     = Float.valueOf(args[2]);
		is2D	= Boolean.valueOf(args[3]);

        /*
         set paths for outputs: same folder, keep the name of the original, with prefixes added
         */
        String parentDir = fileSWC.getParent() + File.separator;          // mother directory
        String name = fileSWC.getName().substring(0, fileSWC.getName().length()-4);

        // new names for outputs
		String pathOutSwc = parentDir + "REC_" + name + ".swc";
		String pathOutTif = parentDir + "IMG_" + name + ".tif";
		String pathOutBif = parentDir + "BIF_" + name + ".swc";
		String pathOutEnd = parentDir + "END_" + name + ".swc";

//		System.out.println("generating... ");
//
//		System.out.println("k    = "+k);
//		System.out.println("snr  = "+snr);
//		System.out.println("is2D = "+is2D);
//
//		System.out.println(pathOutSwc);
//		System.out.println(pathOutTif);
//		System.out.println(pathOutBif);
//		System.out.println(pathOutEnd);
//		System.out.println("done.");

        /*
        generate image
         */
		GeneratorSwc neuronGenerator = new GeneratorSwc();
		ImageStack isOut = neuronGenerator.swc2image(pathInSwc, is2D, k, snr, pathOutSwc, pathOutBif, pathOutEnd);

		//name
//		String outName = new File(pathInSwc).getName();
//		outName = outName.substring(outName.length()-4);

		ImagePlus imOut = new ImagePlus("syn_"+name+"_snr_"+IJ.d2s(snr,1)+"_k_"+IJ.d2s(k,1)+"_is2d_"+is2D, isOut);
		imOut.show();
		FileSaver fs = new FileSaver(imOut);
		if (is2D)
			fs.saveAsTiff(pathOutTif);
		else
			fs.saveAsTiffStack(pathOutTif);
		System.out.println(pathOutTif+" exported");

	}

    /*
        IJ call
     */
    public void run(String s) {

		String pathInSwc =         Prefs.get("critpoint.generate.pathInSwc", System.getProperty("user.home"));
		float k 		 = (float) Prefs.get("critpoint.generate.k", 1);
		float snr 		 = (float) Prefs.get("critpoint.generate.snr", 3);
		boolean is2D	 = Prefs.get("critpoint.generate.is2D", true);

        // parameters
		GenericDialog gd = new GenericDialog("GENERATE 3D NEURON");
		gd.addStringField("swc                  :", pathInSwc, 100);
		gd.addNumericField("k (gauss sigma=k*r) :", k,  2);
		gd.addNumericField("snr                 :", snr,  2);
		gd.addCheckbox("2D neuron", 					is2D);

		gd.showDialog();
		if (gd.wasCanceled()) return;

		pathInSwc 	= gd.getNextString();
        k       	= (float) gd.getNextNumber();
        snr     	= (float) gd.getNextNumber();
		is2D 		= gd.getNextBoolean();

		Prefs.set("critpoint.generate.pathInSwc", 	pathInSwc);
		Prefs.set("critpoint.generate.k", 			k);
		Prefs.set("critpoint.generate.snr", 		snr);
		Prefs.set("critpoint.generate.is2D", 		is2D);

		/*
            check if the swc input file exists
         */
		//pathInSwc = new File(pathInSwc).getAbsolutePath();
		File fileSWC = new File(pathInSwc);

		if (!fileSWC.exists()) {
			IJ.log("file " + pathInSwc + " does not exist!");
			return;
		}

		/*
         set paths to outputs, same folder, keep the name with prefixes added
         */
		String parentDir = fileSWC.getParent() + File.separator;
		String name = fileSWC.getName().substring(0, fileSWC.getName().length()-4);

		// new names for outputs
		String pathOutSwc = parentDir + "REC_" + name + ".swc";
		String pathOutTif = parentDir + "IMG_" + name + ".tif";
		String pathOutBif = parentDir + "BIF_" + name + ".swc";
		String pathOutEnd = parentDir + "END_" + name + ".swc";

		/*
        generate image
         */

        GeneratorSwc neuronGenerator = new GeneratorSwc();
        ImageStack isOut = neuronGenerator.swc2image(pathInSwc, is2D, k, snr, pathOutSwc, pathOutBif, pathOutEnd);

		//name
		//String outName = new File(pathInSwc).getName();
		//outName = outName.substring(outName.length()-4);

		ImagePlus imOut = new ImagePlus("syn_"+name+"_snr_"+IJ.d2s(snr,1)+"_k_"+IJ.d2s(k,1)+"_is2d_"+is2D, isOut);
		imOut.show();

		FileSaver fs = new FileSaver(imOut);
		if (is2D)
			fs.saveAsTiff(pathOutTif);
		else
			fs.saveAsTiffStack(pathOutTif);
		System.out.println(pathOutTif+" exported");


    }

//
//
//        FileSaver fs = new FileSaver(imOut);
//        String outPath = System.getProperty("user.home")+File.separator+"test.tif";
//        fs.saveAsTiffStack(outPath);

}
