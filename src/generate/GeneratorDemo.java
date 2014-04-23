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
 */
public class GeneratorDemo implements PlugIn {

	float p1, p2, p3, Dmin, SNR;
	int N;

	// java -cp critpoint_.jar:/home/miroslav/jarlib/ij.jar -run "GeneratorDemo"
    // generates 2D images with synthetic bifurcations using given parameters
    // images and annotations are stored in separate Ci folders where i marks one configuration

//	public static void main(String args[]){
//
//		/*
//		System.out.println("main() "+args.length);
//		for (int i = 0; i < args.length; i++) System.out.println("ARGS["+i+"]="+args[i]);
//		*/
//
//		if(args.length!=5) {
//			System.out.println("usage:\nSNR\t\tD1\t\tD2\t\tD3\t\tNrImages\n");
//			return;
//		}
//
//		int     SNR     = Integer.valueOf(args[0]);
//		float   D1 	    = Float.valueOf(args[1]);
//		float   D2 	    = Float.valueOf(args[2]);
//		float   D3 	    = Float.valueOf(args[3]);
//		int	    N       = Integer.valueOf(args[4]);
//
//		float[] snr = new float[]{SNR};
//
//		synthetize(snr, D1, D2, D3, N);
//
//	}

	public void run(String s) {

		GenericDialog gd = new GenericDialog("Generate junctions");
		gd.addNumericField("SNR ", 3, 0);
		gd.addNumericField("p1 ", 1, 1);
		gd.addNumericField("p2 ", 1, 1);
		gd.addNumericField("p3 ", 1, 1);
		gd.addNumericField("Dmin ", 3, 0);
		gd.addNumericField("N imgs ", 5, 0);

		gd.showDialog();
		if (gd.wasCanceled()) return;

		SNR     = (float) gd.getNextNumber();
		   p1 	    = (float) gd.getNextNumber();
		   p2 	    = (float) gd.getNextNumber();
		   p3 	    = (float) gd.getNextNumber();
		   Dmin 	    = (float) gd.getNextNumber();
		     N 		= (int) gd.getNextNumber();

		float p = p1 + p2 + p3;

		p1 = p1/p;//(p1<0.2f)? 0.2f : p1;
		p2 = p2 / p;//(p1<0.2f)? 0.2f : p2;
		p3 = p3/p;

//				IJ.log("synthetizing...");

		IJ.log("p="+p);
		IJ.log("p (after)="+(p1+p2+p3));

		float p_min = Math.min(p1, Math.min(p2, p3));
		// Dmin corresponds to min(p1,p2,p3)
		float D1 = Dmin*(p1/p_min);
		float D2 = Dmin*(p2/p_min);
		float D3 = Dmin*(p3/p_min);

		IJ.log("D1 = "+D1);
		IJ.log("D2 = "+D2);
		IJ.log("D3 = "+D3);

		synthetize(SNR, D1, D2, D3, N);



	}

	private static void synthetize(float snr, float d1, float d2, float d3, int N)
	{

		// generate only one category defined with p1,p2,p3 and Dmin -> resulting in Dmax - parameter to set in detection
		// will correspond to one p-ty distribution p1,p2,p3, generate number of images, each with 100 junctions
		// and export the annotaions for bif- and end- points

		String 			timestamp = "";
		DateFormat dateFormat = new SimpleDateFormat("yyyy_MM_dd__HH_mm_ss");
		timestamp =  dateFormat.format(new Date());

		float d_max = Math.max(d1, Math.max(d2, d3));

		String outDir  = "synthetic_junctions(Dmax_SNR),"+IJ.d2s(d_max,1)+","+snr+File.separator;
		IJ.log("synthetize bifurcations... \n" + createDir(outDir));		

		// parameters used fo generate
		try {
			PrintWriter writer = new PrintWriter(outDir+"params.dat");
			writer.println("SNR\t= "+snr);
			writer.println("D1\t= "+d1);
			writer.println("D2\t= "+d2);
			writer.println("D3\t= "+d3);
			writer.println("# imgages nr. \t= "+N);
			writer.close();
		}
		catch (FileNotFoundException ex) {}


        Generator imageGenerator = new Generator();

        for (int cnt=0; cnt<N; cnt++) {

            File gnd_tth_end = new File(outDir+File.separator+String.format("%04d", cnt)+".end");
            File gnd_tth_bif = new File(outDir+File.separator+String.format("%04d", cnt)+".bif");
            File image_path  = new File(outDir+File.separator+String.format("%04d", cnt)+".tif");
            ImagePlus imp = imageGenerator.runDisProportional(snr, d1, d2, d3, gnd_tth_bif);
            FileSaver fs = new FileSaver(imp);
            fs.saveAsTiff(image_path.getAbsolutePath());
            IJ.log(image_path.getAbsolutePath() + " created.");

        }

		/****************************************************
		 * D1<D2<D3
		 * SNR1:
		 * C1 -> D1
		 * C2 -> D1 D2
		 * C3 -> D2
		 * C4 -> D1 D2 D3
		 * C5 -> D1 D3
		 * C6 -> D2 D3
		 * C7 -> D3
		 ****************************************************/

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
				boolean success = f.mkdir();
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
				boolean success = f.mkdir();
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
