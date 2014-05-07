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
 * java -cp critpoint_.jar:/home/miroslav/jarlib/ij.jar -run "GeneratorDemo"
 * generates 2D images with synthetic bifurcations using given parameters
 * images and annotations are stored in separate Ci folders where i marks one configuration
 */
public class GeneratorDemo implements PlugIn {

	float p1, p2, p3, sc, Dmin, SNR;
	int N;
    static int nrP1 = 5;
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
        gd.addMessage(""+nrP1+"x"+nrP1+" junctions per img");
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

        p1 = p1 / p;
		p2 = p2 / p;
        p3 = p3 / p;

		float p_min = Math.min(p1, Math.min(p2, p3));

		// Dmin corresponds to min(p1,p2,p3)
		float D1 = Dmin*(p1/p_min);
		float D2 = Dmin*(p2/p_min);
		float D3 = Dmin*(p3/p_min);

		synthetize(SNR, D1, D2, D3, sc, N);

	}

	private static void synthetize(float snr, float d1, float d2, float d3, float scale, int N)
	{

		// generate only one category defined with p1,p2,p3 and Dmin -> resulting in Dmax - parameter to set in detection
		// will correspond to one p-ty distribution p1,p2,p3, generate number of images, each with 100 junctions
		// and export the annotaions for bif- and end- points

		float d_max = Math.max(d1, Math.max(d2, d3));

		outDir  += File.separator+"synthetic_junctions(Dmax_SNR_s),"+IJ.d2s(d_max,1)+","+snr+","+IJ.d2s(scale,1)+File.separator;
        createDir(outDir);
//		IJ.log(outDir);
        Generator imageGenerator = new Generator(nrP1);

        for (int cnt=0; cnt<N; cnt++) {

            File gnd_tth_end = new File(outDir+File.separator+String.format("%04d", cnt)+".end");
            File gnd_tth_bif = new File(outDir+File.separator+String.format("%04d", cnt)+".bif");
            File gnd_tth_non = new File(outDir+File.separator+String.format("%04d", cnt)+".non");
            File image_path  = new File(outDir+File.separator+String.format("%04d", cnt)+".tif");
            File out_log     = new File(outDir+File.separator+String.format("params")+".csv");

            imageGenerator.runDisProportional(snr, d1, d2, d3, scale, gnd_tth_bif, gnd_tth_end, gnd_tth_non, out_log, image_path);

        }

//		for (int snrIdx = 0; snrIdx < snr.length; snrIdx++) {
//
//			/****************************************************
//			N synthetic images, category 1 (D1)
//			 ***************************************************/
//			System.out.println("category 1. 3x"+p1+", snr "+snr[snrIdx]);
//
//			dirName = outDir+File.separator+"C1";//"SNR_"+snr[snrIdx]+",D["+D1+"]";
//			createDir(dirName);
//
//			for (int cnt=0; cnt<N; cnt++) {
//
//				csvFile = new File(dirName+File.separator+String.format("%04d", cnt)+".csv");
//				// empty file
//				try {
//					PrintWriter writer = new PrintWriter(csvFile);
//					writer.print("# GROUND TRUTH");
//					writer.close();
//				}
//				catch (FileNotFoundException ex) {}
//
//				ImagePlus imp = imageGenerator.runProportional(snr[snrIdx], D1, csvFile);
//				//imp.show();
//				FileSaver fs = new FileSaver(imp);
//				imagePath 	= new File(dirName+File.separator+String.format("%04d", cnt)+".tif").getAbsolutePath();
//				fs.saveAsTiff(imagePath);
//				System.out.println(imagePath);
//
//			}
//
//			/****************************************************
//			N synthetic images, category 3 (D2)
//			 ***************************************************/
//			System.out.println("category 3. 3x"+D2+", snr "+snr[snrIdx]);
//
//			dirName = outDir+File.separator+"C3";//"SNR_"+snr[snrIdx]+",D["+D2+"]";
//			createDir(dirName);
//
//			for (int cnt=0; cnt<N; cnt++) {
//
//				csvFile = new File(dirName+File.separator+String.format("%04d", cnt)+".csv");
//				// empty file
//				try {
//					PrintWriter writer = new PrintWriter(csvFile);
//					writer.print("#");
//					writer.close();
//				}
//				catch (FileNotFoundException ex) {}
//
//				ImagePlus imp = imageGenerator.runProportional(snr[snrIdx], D2, csvFile);
//				//imp.show();
//				FileSaver fs = new FileSaver(imp);
//				imagePath 	= new File(dirName+File.separator+String.format("%04d", cnt)+".tif").getAbsolutePath();
//				fs.saveAsTiff(imagePath);
//				System.out.println(imagePath);
//
//			}
//
//			/****************************************************
//			 N synthetic images, category 7 (D3)
//			 ***************************************************/
//			System.out.println("category 7. 3x"+D3+", snr "+snr[snrIdx]);
//
//			dirName = outDir+File.separator+"C7";//"SNR_"+snr[snrIdx]+",D["+D3+"]";
//			createDir(dirName);
//
//			for (int cnt=0; cnt<N; cnt++) {
//
//				csvFile = new File(dirName+File.separator+String.format("%04d", cnt)+".csv");
//				// empty file
//				try {
//					PrintWriter writer = new PrintWriter(csvFile);
//					writer.print("#");
//					writer.close();
//				}
//				catch (FileNotFoundException ex) {}
//
//				ImagePlus imp = imageGenerator.runProportional(snr[snrIdx], D3, csvFile);
//				//imp.show();
//				FileSaver fs = new FileSaver(imp);
//				imagePath 	= new File(dirName+File.separator+String.format("%04d", cnt)+".tif").getAbsolutePath();
//				fs.saveAsTiff(imagePath);
//				System.out.println(imagePath);
//
//			}
//
//			/****************************************************
//			 N synthetic images, category 2 (D1, D2)
//			 ***************************************************/
//			System.out.println("category 2. "+D1+", "+D2+" snr "+snr[snrIdx]);
//
//			dirName = outDir+File.separator+"C2";//"SNR_"+snr[snrIdx]+",D["+D1+","+D2+"]";
//			createDir(dirName);
//
//			for (int cnt=0; cnt<N; cnt++) {
//
//				csvFile = new File(dirName+File.separator+String.format("%04d", cnt)+".csv");
//				// empty file
//				try {
//					PrintWriter writer = new PrintWriter(csvFile);
//					writer.print("#");
//					writer.close();
//				}
//				catch (FileNotFoundException ex) {}
//
//				ImagePlus imp = imageGenerator.runDisProportional(snr[snrIdx], D1, D2, csvFile);
//				//imp.show();
//				FileSaver fs = new FileSaver(imp);
//				imagePath 	= new File(dirName+File.separator+String.format("%04d", cnt)+".tif").getAbsolutePath();
//				fs.saveAsTiff(imagePath);
//				System.out.println(imagePath);
//
//			}
//
//			/****************************************************
//			 N synthetic images, category 6 (D2, D3)
//			 ***************************************************/
//			System.out.println("category 6. "+D2+", "+D3+" snr "+snr[snrIdx]);
//
//			dirName = outDir+File.separator+"C6";//"SNR_"+snr[snrIdx]+",D["+D2+","+D3+"]";
//			createDir(dirName);
//
//			for (int cnt=0; cnt<N; cnt++) {
//
//				csvFile = new File(dirName+File.separator+String.format("%04d", cnt)+".csv");
//				// empty file
//				try {
//					PrintWriter writer = new PrintWriter(csvFile);
//					writer.print("#");
//					writer.close();
//				}
//				catch (FileNotFoundException ex) {}
//
//				ImagePlus imp = imageGenerator.runDisProportional(snr[snrIdx], D2, D3, csvFile);
//				//imp.show();
//				FileSaver fs = new FileSaver(imp);
//				imagePath 	= new File(dirName+File.separator+String.format("%04d", cnt)+".tif").getAbsolutePath();
//				fs.saveAsTiff(imagePath);
//				System.out.println(imagePath);
//
//			}
//
//			/****************************************************
//			 N synthetic images, category 5 (D1, D3)
//			 ***************************************************/
//			System.out.println("category 5. "+D1+", "+D3+" snr "+snr[snrIdx]);
//
//			dirName = outDir+File.separator+"C5";//"SNR_"+snr[snrIdx]+",D["+D1+","+D3+"]";
//			createDir(dirName);
//
//			for (int cnt=0; cnt<N; cnt++) {
//
//				csvFile = new File(dirName+File.separator+String.format("%04d", cnt)+".csv");
//				// empty file
//				try {
//					PrintWriter writer = new PrintWriter(csvFile);
//					writer.print("#");
//					writer.close();
//				}
//				catch (FileNotFoundException ex) {}
//
//				ImagePlus imp = imageGenerator.runDisProportional(snr[snrIdx], D1, D3, csvFile);
//				//imp.show();
//				FileSaver fs = new FileSaver(imp);
//				imagePath 	= new File(dirName+File.separator+String.format("%04d", cnt)+".tif").getAbsolutePath();
//				fs.saveAsTiff(imagePath);
//				System.out.println(imagePath);
//
//			}
//
//			/****************************************************
//			 N synthetic images, category 4 (D1, D2, D3)
//			 ***************************************************/
//			System.out.println("category 4. "+D1+", "+D2+", "+D3+" snr "+snr[snrIdx]);
//
//			dirName = outDir+File.separator+"C4";//"SNR_"+snr[snrIdx]+",D["+D1+","+D2+","+D3+"]";
//			createDir(dirName);
//
//			for (int cnt=0; cnt<N; cnt++) {
//
//				csvFile = new File(dirName+File.separator+String.format("%04d", cnt)+".csv");
//				// empty file
//				try {
//					PrintWriter writer = new PrintWriter(csvFile);
//					writer.print("#");
//					writer.close();
//				}
//				catch (FileNotFoundException ex) {}
//
//				ImagePlus imp = imageGenerator.runDisProportional(snr[snrIdx], D1, D2, D3, csvFile);
//				//imp.show();
//				FileSaver fs = new FileSaver(imp);
//				imagePath 	= new File(dirName+File.separator+String.format("%04d", cnt)+".tif").getAbsolutePath();
//				fs.saveAsTiff(imagePath);
//				System.out.println(imagePath);
//
//			}
//
//		} // loop snrs

//		String syntheticImgsDir = new File(outDir).getAbsolutePath();

		//Prefs.set("advantra.critpont.synthetic_data_dir", syntheticImgsDir);

//		System.out.println("finished! it's here: \n" + syntheticImgsDir + "\n");

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
