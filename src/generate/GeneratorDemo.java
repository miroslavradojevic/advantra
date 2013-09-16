package generate;

import ij.ImagePlus;
import ij.Prefs;
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

	// java -cp critpoint_.jar:/home/miroslav/jarlib/ij.jar generate.GeneratorDemo arguments (will be called from generate_bifs)

	public static void main(String args[]){

		/*
		System.out.println("main() "+args.length);
		for (int i = 0; i < args.length; i++) System.out.println("ARGS["+i+"]="+args[i]);
		*/

		if(args.length!=5) {
			System.out.println("usage:\nSNR\t\tD1\t\tD2\t\tD3\t\tNrImages\n");
			return;
		}

		int min_SNR = Integer.valueOf(args[0]);
//		int max_SNR = Integer.valueOf(args[1]);
		float D1 	= Float.valueOf(args[1]);
		float D2 	= Float.valueOf(args[2]);
		float D3 	= Float.valueOf(args[3]);
		int	  N     = Integer.valueOf(args[4]);

//		float[] snr = new float[max_SNR-min_SNR+1];
//		for (int i = 0; i < snr.length; i++) snr[i] = min_SNR + i;

		float[] snr = new float[]{min_SNR};

		synthetize(snr, D1, D2, D3, N);

	}

	public void run(String s) {

		GenericDialog gd = new GenericDialog("Generate bifurcations");
		gd.addNumericField("min SNR ", 2, 1);
//		gd.addNumericField("max SNR ", 2, 1);
		gd.addNumericField("D1 ", 2, 1);
		gd.addNumericField("D2 ", 2, 1);
		gd.addNumericField("D3 ", 2, 1);
		gd.addNumericField("N imgs ", 0, 1);

		gd.showDialog();
		if (gd.wasCanceled()) return;

		int min_SNR = (int) gd.getNextNumber();
//		int max_SNR = (int) gd.getNextNumber();
		float D1 	= (float) gd.getNextNumber();
		float D2 	= (float) gd.getNextNumber();
		float D3 	= (float) gd.getNextNumber();
		int N 		= (int) gd.getNextNumber();

//		float[] snr = new float[max_SNR-min_SNR+1];
//		for (int i = 0; i < snr.length; i++) snr[i] = min_SNR + i;
		float[] snr = new float[]{min_SNR};

		synthetize(snr, D1, D2, D3, N);

	}

	private static void synthetize(float[] snr, float D1, float D2, float D3, int N)
	{

		String 			timestamp = "";
		DateFormat dateFormat = new SimpleDateFormat("yyyy_MM_dd__HH_mm_ss");
		timestamp =  dateFormat.format(new Date());

		String outDir  = "synth_bifs_"+timestamp+File.separator;
		System.out.println("synthetize bifurcations...   "+createDir(outDir));

		// parameters used fo generate
		try {
			PrintWriter writer = new PrintWriter(outDir+"params.dat");
			writer.println("SNR\t= "+snr[0]);
			writer.println("D1\t= "+D1);
			writer.println("D2\t= "+D2);
			writer.println("D3\t= "+D3);
			writer.println("# img\t= "+N);
			writer.println("# FOLDER NAMES:");
			writer.println("# C1 -> D1");
			writer.println("# C2 -> D1 D2");
			writer.println("# C3 -> D2");
			writer.println("# C4 -> D1 D2 D3");
			writer.println("# C5 -> D1 D3");
			writer.println("# C6 -> D2 D3");
			writer.println("# C7 -> D3");
			writer.close();
		}
		catch (FileNotFoundException ex) {}

		Generator imageGenerator = new Generator();
		String imagePath, dirName;
		File csvFile;

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

		for (int snrIdx = 0; snrIdx < snr.length; snrIdx++) {

			/****************************************************
			N synthetic images, category 1 (D1)
			 ***************************************************/
			System.out.println("category 1. 3x"+D1+", snr "+snr[snrIdx]);

			dirName = outDir+File.separator+"C1";//"SNR_"+snr[snrIdx]+",D["+D1+"]";
			createDir(dirName);

			for (int cnt=0; cnt<N; cnt++) {

				csvFile = new File(dirName+File.separator+String.format("%04d", cnt)+".csv");
				// empty file
				try {
					PrintWriter writer = new PrintWriter(csvFile);
					writer.print("# GROUND TRUTH");
					writer.close();
				}
				catch (FileNotFoundException ex) {}

				ImagePlus imp = imageGenerator.runProportional(snr[snrIdx], D1, csvFile);
				//imp.show();
				FileSaver fs = new FileSaver(imp);
				imagePath 	= new File(dirName+File.separator+String.format("%04d", cnt)+".tif").getAbsolutePath();
				fs.saveAsTiff(imagePath);
				System.out.println(imagePath);

			}

			/****************************************************
			N synthetic images, category 3 (D2)
			 ***************************************************/
			System.out.println("category 3. 3x"+D2+", snr "+snr[snrIdx]);

			dirName = outDir+File.separator+"C3";//"SNR_"+snr[snrIdx]+",D["+D2+"]";
			createDir(dirName);

			for (int cnt=0; cnt<N; cnt++) {

				csvFile = new File(dirName+File.separator+String.format("%04d", cnt)+".csv");
				// empty file
				try {
					PrintWriter writer = new PrintWriter(csvFile);
					writer.print("#");
					writer.close();
				}
				catch (FileNotFoundException ex) {}

				ImagePlus imp = imageGenerator.runProportional(snr[snrIdx], D2, csvFile);
				//imp.show();
				FileSaver fs = new FileSaver(imp);
				imagePath 	= new File(dirName+File.separator+String.format("%04d", cnt)+".tif").getAbsolutePath();
				fs.saveAsTiff(imagePath);
				System.out.println(imagePath);

			}

			/****************************************************
			 N synthetic images, category 7 (D3)
			 ***************************************************/
			System.out.println("category 7. 3x"+D3+", snr "+snr[snrIdx]);

			dirName = outDir+File.separator+"C7";//"SNR_"+snr[snrIdx]+",D["+D3+"]";
			createDir(dirName);

			for (int cnt=0; cnt<N; cnt++) {

				csvFile = new File(dirName+File.separator+String.format("%04d", cnt)+".csv");
				// empty file
				try {
					PrintWriter writer = new PrintWriter(csvFile);
					writer.print("#");
					writer.close();
				}
				catch (FileNotFoundException ex) {}

				ImagePlus imp = imageGenerator.runProportional(snr[snrIdx], D3, csvFile);
				//imp.show();
				FileSaver fs = new FileSaver(imp);
				imagePath 	= new File(dirName+File.separator+String.format("%04d", cnt)+".tif").getAbsolutePath();
				fs.saveAsTiff(imagePath);
				System.out.println(imagePath);

			}

			/****************************************************
			 N synthetic images, category 2 (D1, D2)
			 ***************************************************/
			System.out.println("category 2. "+D1+", "+D2+" snr "+snr[snrIdx]);

			dirName = outDir+File.separator+"C2";//"SNR_"+snr[snrIdx]+",D["+D1+","+D2+"]";
			createDir(dirName);

			for (int cnt=0; cnt<N; cnt++) {

				csvFile = new File(dirName+File.separator+String.format("%04d", cnt)+".csv");
				// empty file
				try {
					PrintWriter writer = new PrintWriter(csvFile);
					writer.print("#");
					writer.close();
				}
				catch (FileNotFoundException ex) {}

				ImagePlus imp = imageGenerator.runDisProportional(snr[snrIdx], D1, D2, csvFile);
				//imp.show();
				FileSaver fs = new FileSaver(imp);
				imagePath 	= new File(dirName+File.separator+String.format("%04d", cnt)+".tif").getAbsolutePath();
				fs.saveAsTiff(imagePath);
				System.out.println(imagePath);

			}

			/****************************************************
			 N synthetic images, category 6 (D2, D3)
			 ***************************************************/
			System.out.println("category 6. "+D2+", "+D3+" snr "+snr[snrIdx]);

			dirName = outDir+File.separator+"C6";//"SNR_"+snr[snrIdx]+",D["+D2+","+D3+"]";
			createDir(dirName);

			for (int cnt=0; cnt<N; cnt++) {

				csvFile = new File(dirName+File.separator+String.format("%04d", cnt)+".csv");
				// empty file
				try {
					PrintWriter writer = new PrintWriter(csvFile);
					writer.print("#");
					writer.close();
				}
				catch (FileNotFoundException ex) {}

				ImagePlus imp = imageGenerator.runDisProportional(snr[snrIdx], D2, D3, csvFile);
				//imp.show();
				FileSaver fs = new FileSaver(imp);
				imagePath 	= new File(dirName+File.separator+String.format("%04d", cnt)+".tif").getAbsolutePath();
				fs.saveAsTiff(imagePath);
				System.out.println(imagePath);

			}

			/****************************************************
			 N synthetic images, category 5 (D1, D3)
			 ***************************************************/
			System.out.println("category 5. "+D1+", "+D3+" snr "+snr[snrIdx]);

			dirName = outDir+File.separator+"C5";//"SNR_"+snr[snrIdx]+",D["+D1+","+D3+"]";
			createDir(dirName);

			for (int cnt=0; cnt<N; cnt++) {

				csvFile = new File(dirName+File.separator+String.format("%04d", cnt)+".csv");
				// empty file
				try {
					PrintWriter writer = new PrintWriter(csvFile);
					writer.print("#");
					writer.close();
				}
				catch (FileNotFoundException ex) {}

				ImagePlus imp = imageGenerator.runDisProportional(snr[snrIdx], D1, D3, csvFile);
				//imp.show();
				FileSaver fs = new FileSaver(imp);
				imagePath 	= new File(dirName+File.separator+String.format("%04d", cnt)+".tif").getAbsolutePath();
				fs.saveAsTiff(imagePath);
				System.out.println(imagePath);

			}

			/****************************************************
			 N synthetic images, category 4 (D1, D2, D3)
			 ***************************************************/
			System.out.println("category 4. "+D1+", "+D2+", "+D3+" snr "+snr[snrIdx]);

			dirName = outDir+File.separator+"C4";//"SNR_"+snr[snrIdx]+",D["+D1+","+D2+","+D3+"]";
			createDir(dirName);

			for (int cnt=0; cnt<N; cnt++) {

				csvFile = new File(dirName+File.separator+String.format("%04d", cnt)+".csv");
				// empty file
				try {
					PrintWriter writer = new PrintWriter(csvFile);
					writer.print("#");
					writer.close();
				}
				catch (FileNotFoundException ex) {}

				ImagePlus imp = imageGenerator.runDisProportional(snr[snrIdx], D1, D2, D3, csvFile);
				//imp.show();
				FileSaver fs = new FileSaver(imp);
				imagePath 	= new File(dirName+File.separator+String.format("%04d", cnt)+".tif").getAbsolutePath();
				fs.saveAsTiff(imagePath);
				System.out.println(imagePath);

			}

		} // loop snrs

		String syntheticImgsDir = new File(outDir).getAbsolutePath();

		//Prefs.set("advantra.critpont.synthetic_data_dir", syntheticImgsDir);

		System.out.println("finished! it's here: \n" + syntheticImgsDir + "\n");

	}

	public static String createDir(String dir_name){

		File f = new File(dir_name);

		if(!f.exists()){

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
