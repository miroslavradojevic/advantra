package generate;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.io.FileSaver;
import ij.io.OpenDialog;
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

	String      swc_path;
	String		out_dir;

    public void run(String s) {

		/*
        load the swc through the menu
         */
		String in_folder = Prefs.get("id.folder", System.getProperty("user.home"));
		OpenDialog.setDefaultDirectory(in_folder);
		OpenDialog dc = new OpenDialog("Select file");
		in_folder = dc.getDirectory();
		swc_path = dc.getPath();
		if (swc_path==null || swc_path.substring(swc_path.length()-4, swc_path.length()).equals(".swc")) return;
		Prefs.set("id.folder", in_folder);

		File fileSWC = new File(swc_path);
		if(fileSWC==null) return;

		float k 		 = (float) Prefs.get("critpoint.generate.k", 1);
		float snr 		 = (float) Prefs.get("critpoint.generate.snr", 3);
		boolean is2D	 = Prefs.get("critpoint.generate.is2D", true);
		String pathInSwc =         Prefs.get("critpoint.generate.pathInSwc", System.getProperty("user.home"));

		GenericDialog gd = new GenericDialog("GENERATE 3D NEURON");
		gd.addStringField("swc                  :", pathInSwc, 100);
		gd.addNumericField("k (gauss sigma=k*r) :", k,  2);
		gd.addNumericField("SNR", snr,  2);
		gd.addCheckbox("2D", 					is2D);

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
