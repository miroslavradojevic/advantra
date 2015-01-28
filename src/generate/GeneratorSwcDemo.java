package generate;

import aux.ReadSWC;
import aux.Stat;
import aux.Tools;
import ij.*;
import ij.gui.GenericDialog;
import ij.io.FileSaver;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 10/23/13
 * Time: 4:49 PM
 * plugin that calls GeneratorSwc class that would use .swc files to create 2d or 3d image stacks (8bit, tif)
 * synthetic images of neurons with possibility of having different snr
 * ground truth critical points are inferred from the swc itself automatically
 * input : path to swc,
 * output: tif image stack or image, critical point swc, end point swc
 */
public class GeneratorSwcDemo implements PlugIn {

	String      swc_path;   	// input swc path
	String		out_dir;    	// output directory
	float 		SNR;        	// snr level of the generated file
	boolean		is2D;       	// generate 2d or 3d image stack
	float 		um_pix_xy;  	// size of one pixel in micro meters in xy (used to convert micro-meter values into pixels) -- will be used for radius as well (?!)
	float 		um_pix_z;   	// size of one pixel in micro meters in z
	boolean     auto_um_to_pix; // select the pixel size automatically so that it fits NpixPerDiam per diameter

	// some constants
	static float 			pixPerDiam = 3f; 	// used in automatic selection of the pixel size, generation constant - at least this many pixels are used to cover the (average/median) diameter

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

		ReadSWC readerSwc = new ReadSWC(swc_path); // pulled out so that expected diameter can be known in advance
		float average_diam = readerSwc.averageDiameter(); //IJ.log("a.d. "+average_diam+"");
		float median_diam = readerSwc.medianDiameter();

		// take the parameters
		if (Macro.getOptions()==null) {

			SNR 		= (float) Prefs.get("critpoint.generate.snr", 3);
			is2D	 	=         Prefs.get("critpoint.generate.is2D", true);
			out_dir		= 		  Prefs.get("critpoint.generate.out_dir", System.getProperty("user.home"));
			um_pix_xy   = (float) Prefs.get("critpoint.generate.um_pix_xy", 0.1f);
			um_pix_z	= (float) Prefs.get("critpoint.generate.um_pix_z",  0.5f);
			auto_um_to_pix =      Prefs.get("critpoint.generate.auto_um_to_pix", true);

			GenericDialog gd = new GenericDialog("GENERATE NEURON FROM SWC");

			gd.addNumericField("SNR", 				SNR,  		2);
			gd.addCheckbox("2D", 			is2D);
			gd.addStringField("out_dir", 	out_dir, 50);
			gd.addMessage("average/median neurite diameter: " + IJ.d2s(average_diam,2)+"/"+IJ.d2s(median_diam,2)+"[um]");
			gd.addNumericField("x,y:", 		um_pix_xy,  2, 	4, "[um/pix]");
			gd.addNumericField("z:  ", 	 	um_pix_z, 	2,	4, "[um/pix]");
			gd.addCheckbox("automatic ("+IJ.d2s(pixPerDiam,1)+"  pix/avg.diameter)", 	auto_um_to_pix);

			gd.showDialog();
			if (gd.wasCanceled()) return;

			SNR     	= (float) gd.getNextNumber();
			is2D 		= gd.getNextBoolean();
			out_dir  	= gd.getNextString();
			auto_um_to_pix = gd.getNextBoolean();

			Prefs.set("critpoint.generate.snr", 		SNR);
			Prefs.set("critpoint.generate.is2D", 		is2D);
			Prefs.set("critpoint.generate.out_dir",		out_dir);
			Prefs.set("critpoint.generate.um_pix_xy", 	um_pix_xy);
			Prefs.set("critpoint.generate.um_pix_z", 	um_pix_z);
			Prefs.set("critpoint.generate.auto_um_to_pix", auto_um_to_pix);
		}
		else {

			SNR     	= Float.valueOf(Macro.getValue(Macro.getOptions(), "SNR", String.valueOf(3)));
			is2D 		= Boolean.valueOf(Macro.getValue(Macro.getOptions(), "is2D", String.valueOf(true)));
			out_dir  	= Macro.getValue(Macro.getOptions(), "out_dir", System.getProperty("user.home"));
			um_pix_xy 	= Float.valueOf(Macro.getValue(Macro.getOptions(), "um_pix_xy", String.valueOf(0.1f)));
			um_pix_z 	= Float.valueOf(Macro.getValue(Macro.getOptions(), "um_pix_z", String.valueOf(0.5f)));
			auto_um_to_pix = Boolean.valueOf(Macro.getValue(Macro.getOptions(), "auto_um_to_pix", String.valueOf(true)));

		}

		//	set paths to outputs, same folder, keep the name with prefixes added
		//String parentDir = fileSWC.getParent() + File.separator;
		String name = swc_file.getName().substring(0, swc_file.getName().length()-4);
		// replace '_' with '-' in order not to make confusion later on (_ separates parameters stored in the name - according to nomenclature)
		char[] name_chars = name.toCharArray();
		for (int ii=0; ii<name_chars.length; ii++) if (name_chars[ii] == '_') name_chars[ii] = '-';
		name = String.valueOf(name_chars);

		// output folder name
		boolean ends_with_sep =  out_dir.substring(out_dir.length()-1, out_dir.length()).equals(File.separator);
		out_dir  += ((!ends_with_sep)? File.separator : "")+
						   "SWCGEN.SNR_"+
						   String.format("%.1f", SNR)+
						   File.separator;
		if ( !(new File(out_dir).exists()) ) {
			Tools.createDir(out_dir); // create directory if it does not exist, otherwise just append in
		}

        // swc radius correction (swc has irregular radiuses and also generating tiny radiuses does not make sense)
//		float[] all_radiuses = new float[readerSwc.nodes.size()];
//        for (int ii=0; ii<readerSwc.nodes.size(); ii++)
//            all_radiuses[ii] = RMIN + (readerSwc.nodes.get(ii)[readerSwc.RADIUS] - readerSwc.minR);
//        float diam_est = 2 * Stat.average(all_radiuses);

		// new names for outputs
        String out_name = String.format("%s", name);//, diam_est);
//        File gnd_tth_end = new File(out_dir+out_name+".end");
        File rec_swc = new File(out_dir+out_name+"_rec"+".swc");
//        File gnd_tth_non = new File(out_dir+out_name+".non");
        File image_path  = new File(out_dir+out_name+".tif");
        File gnd_tth_swc = new File(out_dir+out_name+".swc"); // "_gndtth"

		System.out.println("\n" + swc_path);

		// NOTE: x,y,z and r are in micrometers!
		// all the values stored in readerSWC are in um, while generating works in pix
		// use constant um/pix to convert those into pixels - because pixels are used as metric when generating

		if (auto_um_to_pix) {
			um_pix_xy = median_diam / pixPerDiam; // redefine, override old pixel size
			um_pix_z  = um_pix_xy; // z gets the same pixel size - maybe not ok!!!
		}

		readerSwc.umToPix(um_pix_xy, um_pix_z); // convert it before plugged into generator, radiuses are scaled as well
		IJ.log(" " + IJ.d2s(um_pix_xy, 2) + " um/pix");

		readerSwc.lowRadiusBoundary(1.5f); //

		GeneratorSwc neuronGenerator = new GeneratorSwc();

		// clean output files if they exist
		if(gnd_tth_swc.exists()) try {
			new PrintWriter(gnd_tth_swc.getAbsolutePath()).close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		if(rec_swc.exists()) try {
			new PrintWriter(rec_swc.getAbsolutePath()).close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		neuronGenerator.swc2image(readerSwc, is2D, SNR,
				rec_swc,
				gnd_tth_swc,
				image_path
		);
		// rec_swc: output reconstruction, for the generated image
		// gnd_tth_swc: output swc with critical points
		// image_path: output image path

    }

}
