package generate;

import ij.ImagePlus;
import ij.plugin.PlugIn;

import java.io.File;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 10/23/13
 * Time: 4:49 PM
 */
public class Generator3DDemo implements PlugIn {

    public static void main(String args[]){

        float   k;
        float   snr;
        File    fileSWC;
        String  pathSWC;

        if (args.length!=3) {
            System.out.println("# GENERATE SYNTHETIC IMAGE FROM SWC #");
            System.out.println("usage: pathSWC k SNR");
            return;
        }

        pathSWC = new File(args[0]).getAbsolutePath();
        fileSWC = new File(pathSWC);

        if (!fileSWC.exists()) {
            System.out.println("file "+pathSWC+" does not exist!");
            return;
        }

        k       = Float.valueOf(args[1]);
        snr     = Float.valueOf(args[2]);

        System.out.println("loaded: \n"+pathSWC+"\n"+k+"\n"+snr);
        System.out.println("creating image...");

        Generator3D g3D = new Generator3D();
        ImagePlus formedImage = g3D.fromSWC(pathSWC, k, snr);
		formedImage.show();

    }

    public void run(String s) {

        // parameters

        String  pathSWC = "/home/miroslav/Copy/set3d/nmorpho/hippoc1.swc";
        float   k       = 1;
        float   snr     = 3;

        Generator3D g3D = new Generator3D();
        ImagePlus formedImage = g3D.fromSWC(pathSWC, k, snr);
        formedImage.show();


    }
}
