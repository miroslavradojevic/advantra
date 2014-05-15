package generate;

import ij.IJ;
import ij.ImagePlus;
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
 * outputs are images (0000.tif) and annotations (0000.bif, 0000.end, 0000.non)
 */
public class GeneratorDemo implements PlugIn {

	float p1, p2, p3, sc, Dmin, SNR;
	int N;
    static int nrP1 = 5;   // number of points that are stored along one dimension
	static String outDir=System.getProperty("user.home");

	public void run(String s) {

		GenericDialog gd = new GenericDialog("Generate junctions");
		gd.addNumericField("SNR:   ", 3, 0);
		gd.addNumericField("p1:    ", 1, 1);
		gd.addNumericField("p2:    ", 1, 1);
		gd.addNumericField("p3:    ", 1, 1);
		gd.addNumericField("scale: ", 8, 1);
		gd.addNumericField("Dmin:  ", 3, 0);
		gd.addNumericField("# imgs:", 1, 0);
//        gd.addMessage(""+nrP1+"x"+nrP1+" junctions per img");
		gd.addStringField("out dir:", outDir, 50);

		gd.showDialog();
		if (gd.wasCanceled()) return;

		SNR     = (float) gd.getNextNumber();
		p1 	    = (float) gd.getNextNumber();
		p2 	    = (float) gd.getNextNumber();
		p3 	    = (float) gd.getNextNumber();
		sc 	    = (float) gd.getNextNumber();
		Dmin 	= (float) gd.getNextNumber();
		N 		= (int) gd.getNextNumber();
		outDir  = gd.getNextString();

		float p = p1 + p2 + p3;
        p1 = p1 / p; // normalize them
		p2 = p2 / p;
        p3 = p3 / p;

		synthetize(SNR, p1, p2, p3, sc, N, Dmin);

	}

	private static void synthetize(float snr, float p1, float p2, float p3, float scale, int N, float Dmin)
	{

		// generate only one category defined with p1,p2,p3 and Dmin -> Dmax can be inferred - essential parameter to set in detection
		// will correspond to one p-ty distribution p1,p2,p3, generate number of images, each with nrP2 junctions
		// and export the annotations for bif- and end- points and non-points

		float p_min = Math.min(p1, Math.min(p2, p3));
		float d1 = Dmin*(p1/p_min);
		float d2 = Dmin*(p2/p_min);
		float d3 = Dmin*(p3/p_min);

		float d_max = Math.max(d1, Math.max(d2, d3));

        // output folder name, values to be later extracted are separated with _
		outDir  += File.separator+
						   "synjun(Dmax,SNR,s,p1,p2,p3)_"+
						   String.format("%d_%d_%.1f_%.2f_%.2f_%.2f", (int)Math.ceil(d_max), Math.round(snr), scale, p1, p2, p3)+
						   File.separator;
        createDir(outDir);

        Generator imageGenerator = new Generator(nrP1);

        for (int cnt=0; cnt<N; cnt++) {

            File gnd_tth_end = new File(outDir+File.separator+String.format("%04d", cnt)+".end");
            File gnd_tth_bif = new File(outDir+File.separator+String.format("%04d", cnt)+".bif");
            File gnd_tth_non = new File(outDir+File.separator+String.format("%04d", cnt)+".non");
            File image_path  = new File(outDir+File.separator+String.format("%04d", cnt)+".tif");
            File out_log     = new File(outDir+File.separator+String.format("params")+".csv");

            imageGenerator.runDisProportional(snr, d1, d2, d3, scale, gnd_tth_bif, gnd_tth_end, gnd_tth_non, out_log, image_path);

        }

	}

	public static String createDir(String dir_name){

		File f = new File(dir_name);

		if(!f.exists()){  //

			try{
				// Create one directory
				boolean success = f.mkdirs();
				if (!success) {
					System.out.println("Error: Directory: " + dir_name + "    .... was NOT created");
				}

			}
			catch (Exception e){//Catch exception if any
				System.err.println("Error: " + e.getMessage());
			}

		}
		else {



            // empty it before deleting
            String[] entries = f.list();

            for(String s: entries){
                File currentFile = new File(f.getPath(),s);
                currentFile.delete();
            }

            f.delete();

			try{
				// Create one directory
				boolean success = f.mkdirs();
				if (!success) {
					System.out.println("Error: Directory: " + dir_name + "    .... was NOT created");
				}

			}
			catch (Exception e){//Catch exception if any
				System.err.println("Error: " + e.getMessage());
			}

		}

		return f.getAbsolutePath();

	}

	public static void createDirs(String dirs_names){
		try{
			// Create multiple directories
			boolean success = (new File(dirs_names)).mkdirs();
			if (success) {
				System.out.println("Directories: " + dirs_names + " created");
			}
		}
		catch (Exception e){//Catch exception if any
			System.err.println("Error: " + e.getMessage());
		}
	}

}
