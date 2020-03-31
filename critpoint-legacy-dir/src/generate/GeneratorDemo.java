package generate;

import aux.Tools;
import ij.IJ;
import ij.ImagePlus;
import ij.Macro;
import ij.gui.GenericDialog;
import ij.io.FileSaver;
import ij.plugin.PlugIn;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 9/6/13
 * Time: 1:50 PM
 * generates 2D images with synthetic bifurcations using given parameters
 * outputs are images (parameter0000.tif) and annotations (parameter0000.swc)
 */
public class GeneratorDemo implements PlugIn {

	float p1, p2, p3, SNR;
	int N;
    static int nrP1 = 7;   // number of points that are stored along one dimension
	static String outDir=System.getProperty("user.home");

	public void run(String s) {

        if (Macro.getOptions()==null) {
            GenericDialog gd = new GenericDialog("Generate junctions");

			gd.addNumericField("SNR:   ", 3, 0);
            gd.addNumericField("p1:    ", 1, 1);
            gd.addNumericField("p2:    ", 1, 1);
            gd.addNumericField("p3:    ", 1, 1);
            //gd.addNumericField("scale: ", 8, 1);
            //gd.addNumericField("Dmin:  ", 3, 0);
            gd.addNumericField("nr_imgs:", 1, 0);
            gd.addStringField("out_dir:", outDir, 50);

            gd.showDialog();
            if (gd.wasCanceled()) return;

            SNR     = (float) gd.getNextNumber();
            p1 	    = (float) gd.getNextNumber();
            p2 	    = (float) gd.getNextNumber();
            p3 	    = (float) gd.getNextNumber();
            //sc 	    = (float) gd.getNextNumber();
            //Dmin 	= (float) gd.getNextNumber();
            N 		= (int) gd.getNextNumber();
            outDir  = gd.getNextString();
        }
        else {
            SNR     = Float.valueOf(Macro.getValue(Macro.getOptions(), "snr", String.valueOf(3)));
            p1 	    = Float.valueOf(Macro.getValue(Macro.getOptions(), "p1", String.valueOf(1)));
            p2 	    = Float.valueOf(Macro.getValue(Macro.getOptions(), "p2", String.valueOf(1)));
            p3 	    = Float.valueOf(Macro.getValue(Macro.getOptions(), "p3", String.valueOf(1)));
            //sc 	    = Float.valueOf(Macro.getValue(Macro.getOptions(), "scale", String.valueOf(8)));
            //Dmin 	= Float.valueOf(Macro.getValue(Macro.getOptions(), "Dmin", String.valueOf(3)));
            N 		= Integer.valueOf(Macro.getValue(Macro.getOptions(), "nr_imgs", String.valueOf(1)));
            outDir  = Macro.getValue(Macro.getOptions(), "out_dir", System.getProperty("user.home"));
        }

		float p = p1 + p2 + p3;
        p1 = p1 / p; // normalize them
		p2 = p2 / p;
        p3 = p3 / p;

        // generate only one category defined with p1,p2,p3 and Dmin -> Dmax can be inferred - essential parameter to set in detection
        // will correspond to one p-ty distribution p1,p2,p3, generate number of images, each with nrP2 junctions
        // and export the annotations for bif- and end- points and non-points

        // Dmin is lowest possible diameter that makes sense
        float Dmin = 3; // fixed minimal radius

        float p_min = Math.min(p1, Math.min(p2, p3));
        float d1 = Dmin*(p1/p_min);
        float d2 = Dmin*(p2/p_min);
        float d3 = Dmin*(p3/p_min);

        float d_max = Math.max(d1, Math.max(d2, d3));
        // output folder/image name, values to be later extracted are separated with _
        String name_pattern = "SYNGEN.Dmax.SNR.p1.p2.p3_" + String.format("%.1f_%.1f_%.2f_%.2f_%.2f", d_max, SNR, p1, p2, p3);

        if (!outDir.endsWith(File.separator)) outDir  += File.separator;// + name_pattern + File.separator;
        if (!new File(outDir).exists()) Tools.createDir(outDir);

        outDir += name_pattern + File.separator;
        if (!new File(outDir).exists()) Tools.createDir(outDir);

        Generator imageGenerator = new Generator(nrP1);

        for (int cnt=0; cnt<N; cnt++) {

            File gnd_tth_swc = new File(outDir+File.separator+name_pattern+String.format("_%04d", cnt)+".swc");
            File image_path  = new File(outDir+File.separator+name_pattern+String.format("_%04d", cnt)+".tif");
            imageGenerator.runDisProportional(SNR, d1, d2, d3, gnd_tth_swc, image_path);

        }

	}

}
