package generate;

import aux.Tools;
import ij.*;
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
	float 		SNR;
	boolean		is2D;

	// some constants
	static float 		DMIN = 3;
	static float 		K = 2.5f;   // scaling the gaussian to match radius

    public void run(String s) {

		// load the swc file
		String in_folder = Prefs.get("id.folder", System.getProperty("user.home"));
		OpenDialog.setDefaultDirectory(in_folder);
		OpenDialog dc = new OpenDialog("Select SWC file");
		in_folder = dc.getDirectory();
		swc_path = dc.getPath();
		if (swc_path==null || (!swc_path.substring(swc_path.length()-4, swc_path.length()).equals(".swc"))) return;
		Prefs.set("id.folder", in_folder);

		File swc_file = new File(swc_path);
		if(swc_file==null) return;

		if (Macro.getOptions()==null) {
			// there was no macro call
			SNR 		= (float) Prefs.get("critpoint.generate.snr", 3);
			is2D	 	=         Prefs.get("critpoint.generate.is2D", true);
			out_dir		= 		  Prefs.get("critpoint.generate.out_dir", System.getProperty("user.home"));

			GenericDialog gd = new GenericDialog("GENERATE NEURON FROM SWC");

			gd.addNumericField("SNR", 				SNR,  2);
			gd.addCheckbox("2D", 					is2D);
			gd.addStringField("out_dir", 			out_dir, 50);

			gd.showDialog();
			if (gd.wasCanceled()) return;

			SNR     	= (float) gd.getNextNumber();
			is2D 		= gd.getNextBoolean();
			out_dir  	= gd.getNextString();

			Prefs.set("critpoint.generate.snr", 		SNR);
			Prefs.set("critpoint.generate.is2D", 		is2D);
			Prefs.set("critpoint.generate.out_dir",		out_dir);
		}
		else {
			// macro is called, parameters are there
			SNR     = Float.valueOf(Macro.getValue(Macro.getOptions(), "SNR", String.valueOf(3)));
			is2D 	= Boolean.valueOf(Macro.getValue(Macro.getOptions(), "is2D", String.valueOf(true)));
			out_dir  = Macro.getValue(Macro.getOptions(), "out_dir", System.getProperty("user.home"));
		}

//		float k 		 = (float) Prefs.get("critpoint.generate.k", 1);
//		String pathInSwc =         Prefs.get("critpoint.generate.pathInSwc", System.getProperty("user.home"));
//		gd.addStringField("swc                  :", pathInSwc, 100);
//		gd.addNumericField("k (gauss sigma=k*r) :", k,  2);
//		pathInSwc 	= gd.getNextString();
//        k       	= (float) gd.getNextNumber();
//		Prefs.set("critpoint.generate.pathInSwc", 	pathInSwc);
//		Prefs.set("critpoint.generate.k", 			k);

		//	set paths to outputs, same folder, keep the name with prefixes added
		//String parentDir = fileSWC.getParent() + File.separator;
		String name = swc_file.getName().substring(0, swc_file.getName().length()-4);
		// replace '_' with '-' in order not to make confusion later on (_ separates parameters stored in the name - according to nomenclature)
		char[] name_chars = name.toCharArray();
		for (int ii=0; ii<name_chars.length; ii++) if (name_chars[ii] == '_') name_chars[ii] = '-';
		name = String.valueOf(name_chars);

		// output folder name
		boolean ends_with_sep =  out_dir.substring(out_dir.length()-1, out_dir.length()).equals(File.separator);
		out_dir  += ((ends_with_sep)?File.separator:"")+
						   "swcgen.SNR_"+
						   String.format("%.1f", SNR)+
						   File.separator;
		Tools.createDir(out_dir);

		// new names for outputs
		String pathOutSwc = out_dir + String.format("%s.Dmed_%.1f", name, );// name + "" + ".swc";

		String pathOutTif = out_dir + "IMG_" + name + ".tif";

		String pathOutBif = parentDir + "BIF_" + name + ".swc";

		String pathOutEnd = parentDir + "END_" + name + ".swc";

		/*
        generate image
         */

        GeneratorSwc neuronGenerator = new GeneratorSwc();
        ImageStack isOut = neuronGenerator.swc2image(pathInSwc, is2D, K, SNR, pathOutSwc, pathOutBif, pathOutEnd);
        ImageStack isOut = neuronGenerator.swc2image(swc_file, is2D, K, DMIN, SNR, pathOutSwc, pathOutBif, pathOutEnd);

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
