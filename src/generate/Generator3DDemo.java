package generate;

import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.io.FileSaver;
import ij.plugin.PlugIn;

import java.io.File;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 10/23/13
 * Time: 4:49 PM
 */
public class Generator3DDemo implements PlugIn {

    /*
        terminal call
     */
    public static void main(String args[]){

        float   k;
        float   snr;
        File    fileSWC;
        String  pathInSwc;

        /*
            check argument nr.
         */
        if (args.length!=3) {
            System.out.println("# GENERATE SYNTHETIC IMAGE FROM SWC #");
            System.out.println("usage: path_to_SWC k SNR");
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
        Generator3D g3D = new Generator3D();
        ImageStack isOut = g3D.swc2Stack(pathInSwc, k, snr, pathOutSwc, pathOutBif, pathOutEnd, pathOutTif);

        FileSaver fs = new FileSaver(new ImagePlus(name, isOut));
        fs.saveAsTiffStack(pathOutTif);

    }

    /*
        IJ call
     */
    public void run(String s) {

        // parameters
		GenericDialog gd = new GenericDialog("GENERATE 3D NEURON");
		gd.addStringField("swc                  :", "path-to-swc", 50);
		gd.addNumericField("k (gauss sigma=k*r) :", 1, 1);
		gd.addNumericField("snr                 :", 3,  1);

		gd.showDialog();
		if (gd.wasCanceled()) return;

		String pathInSwc 	= gd.getNextString();
        float   k       	= (float) gd.getNextNumber();
        float   snr     	= (float) gd.getNextNumber();

		/*
            check if the swc input file exists
         */
		pathInSwc = new File(pathInSwc).getAbsolutePath();
		File fileSWC = new File(pathInSwc);

		if (!fileSWC.exists()) {
			System.out.println("file "+pathInSwc+" does not exist!");
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
        Generator3D g3D = new Generator3D();
        ImageStack isOut = g3D.swc2Stack(pathInSwc, k, snr, pathOutSwc, pathOutBif, pathOutEnd, pathOutTif);


		FileSaver fs = new FileSaver(new ImagePlus(name, isOut));
		fs.saveAsTiffStack(pathOutTif);

//        String outName = new File(pathInSwc).getName();
//        outName = outName.substring(outName.length()-4);
//        FileSaver fs = new FileSaver(imOut);
//        String outPath = System.getProperty("user.home")+File.separator+"test.tif";
//        fs.saveAsTiffStack(outPath);

    }
}
