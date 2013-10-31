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
        String  pathOutSwc="";
        String  pathOutTif="";

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

        System.out.println("loaded: \n"+pathInSwc+"\n"+k+"\n"+snr);


        /*
         set paths to outputs, same folder, keep the name with prefixes added
         */
        String parentDir = fileSWC.getParent() + File.separator;
        String name = fileSWC.getName().substring(0, fileSWC.getName().length()-4);
        // new names for outputs
        pathOutSwc = parentDir + "REC_" + name + ".swc";
        pathOutTif = parentDir + "IMG_" + name + ".tif";
        //System.out.println("outputs = \n" + pathOutSwc + "\n" + pathOutTif);

        /*
        generate image
         */
        Generator3D g3D = new Generator3D();
        ImageStack isOut = g3D.Swc2Stack(pathInSwc, k, snr, pathOutSwc, pathOutTif);

        FileSaver fs = new FileSaver(new ImagePlus(name, isOut));
        fs.saveAsTiffStack(pathOutTif);

    }


    /*
        IJ call
     */
    public void run(String s) {

        // parameters
		//GenericDialog gd = new GenericDialog("GENERATE 3D NEURON");

        float   k       = 1;
        float   snr     = 3;
        File    fileSWC;
        String  pathInSwc = "/home/miroslav/test.swc";
        String  pathOutSwc="";
        String  pathOutTif="";

        Generator3D g3D = new Generator3D();
        ImageStack isOut = g3D.Swc2Stack(pathInSwc, k, snr, pathOutSwc, pathOutTif);

        String outName = new File(pathInSwc).getName();
        outName = outName.substring(outName.length()-4);

        System.out.println(" file name: "+outName);

        ImagePlus imOut = new ImagePlus(outName, isOut);

        // save with appropriate name and export cooresponding swc
        FileSaver fs = new FileSaver(imOut);
        String outPath = System.getProperty("user.home")+File.separator+"test.tif";
        fs.saveAsTiffStack(outPath);
        System.out.println(outPath+" saved.");


    }
}
