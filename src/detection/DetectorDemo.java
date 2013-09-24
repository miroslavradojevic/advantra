package detection;

import aux.AnalyzeCSV;
import generate.GeneratorDemo;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.io.FileSaver;
import ij.plugin.PlugIn;

import java.awt.*;
import java.io.*;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.Vector;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 9/10/13
 * Time: 5:35 PM
 */
public class DetectorDemo implements PlugIn {

	// ok to be static because it's not necessary to instantiate DetectorDemo, it's meant to be one class anyway
	static String 		inDir;
	static String 		outDir;

	static float		neuronD;
	static float		iDiff;


	static float 		minCosAngle 		= .8f;
	static float 		minFuzzyScore 		= .6f;
	static float 		scatterDistSquared 	=  5f;

	static int 			wStdRatioToD 		= 3;
	static int 	 		MIN_SIZE 			= 2;

	static float 	 	LOCATION_TOLERANCE_SCALE 	= 1.5f;

	static String 		logFile 	= "legend.dat"; // string (dir. name) for each category
	static PrintWriter 	logWriter;

	static String 		outFile 	= "data.dat";	// table with data
	static PrintWriter	outWriter;

	static String 		parFile		= "params.data";// file with params used for the detection
	static PrintWriter	parWriter;

	public static void main(String args[]){

		/***** read source directory *****/
		if(args.length!=3) {
			System.out.println("usage:\tinDir\tneuronDiameter\tiDiff");
			return;
		}
		inDir 		= args[0];
		neuronD 	= Float.valueOf(args[1]);
		iDiff 		= Float.valueOf(args[2]);

		/***** check source directory *****/
		File fcheck = new File(inDir);
		if (!fcheck.exists()) {
			System.err.println("Directory does not exist");
			return;
		}
		inDir = fcheck.getAbsolutePath();
		Prefs.set("advantra.critpont.synthetic_data_dir", inDir);
		System.out.println("reading from " + inDir);

		/***** out directory ****/
		File f = new File(inDir);
		if (!f.exists()) {
			System.err.println(inDir+ "     directory does not exist!");
			return;
		}
		else {
			outDir  = f.getParent() + File.separator + "DET_"+ f.getName() + File.separator;
			GeneratorDemo.createDir(outDir);
		}

		startLog(outDir);
		startOut(outDir);
		startPar(outDir);

		detect(inDir);

		finishLog();
		finishOut();
		finishPar();

	}

	public void run(String s) {

		/***** read source directory *****/
		GenericDialog gd = new GenericDialog("DetectorDemo");
		gd.addStringField("folder with synthetic data", Prefs.get("advantra.critpont.synthetic_data_dir", ""), 50);
		gd.showDialog();
		if (gd.wasCanceled()) return;
		inDir = gd.getNextString();

		/***** check source directory *****/
		File fcheck = new File(inDir);
		if (!fcheck.exists()) {
			System.err.println("Directory does not exist");
			return;
		}
		inDir = fcheck.getAbsolutePath();
		Prefs.set("advantra.critpont.synthetic_data_dir", inDir);
		IJ.log("reading from " + inDir);

		/***** out directory ****/
		File f = new File(inDir);
		if (!f.exists()) {
			System.err.println(inDir+ "     directory does not exist!");
			return;
		}
		else {
			outDir  = f.getParent() + File.separator + "DET_"+ f.getName() + File.separator;
			GeneratorDemo.createDir(outDir);
		}

		startLog(outDir);
		startOut(outDir);
		startPar(outDir);

		detect(inDir);

		finishLog();
		finishOut();
		finishPar();
	}

	private static void startPar(String dirPath)
	{

		String outFilePath = dirPath+parFile;

		try {
			PrintWriter writer = new PrintWriter(outFilePath);
			writer.print("");
			writer.close();
		}
		catch (FileNotFoundException ex) {}
		// initialize par file
		try {
			parWriter = new PrintWriter(new BufferedWriter(new FileWriter(outFilePath, true)));
		} catch (IOException e) {}

	}

	private static void startLog(String dirPath)
	{

		String outFilePath = dirPath+logFile;

		try {
			PrintWriter writer = new PrintWriter(outFilePath);
			writer.print("");
			writer.close();
		}
		catch (FileNotFoundException ex) {}

		// initialize log file
		try {
			logWriter = new PrintWriter(new BufferedWriter(new FileWriter(outFilePath, true)));
		} catch (IOException e) {}

	}

	private static void startOut(String dirPath)
	{

		String outFilePath = dirPath+outFile;

		try {
			PrintWriter writer = new PrintWriter(outFilePath);
			writer.print("");
			writer.close();
		}
		catch (FileNotFoundException ex) {}

		// initialize log file
		try {
			outWriter = new PrintWriter(new BufferedWriter(new FileWriter(outFilePath, true)));
		} catch (IOException e) {}

	}

	private static void finishLog()
	{

		logWriter.close();

	}

	private static void finishOut()
	{

		outWriter.close();

	}

	private static void finishPar()
	{

		parWriter.close();

	}

	private static String timestamp()
	{
		// timestamp for the folder name
		DateFormat dateFormat = new SimpleDateFormat("yyyy_MM_dd__HH_mm_ss");
		String timestamp =  dateFormat.format(new Date());
		return timestamp;
	}

	private static void detect(String inDir)
	{

		//System.out.println(outDir);
		//System.out.println(inDir);

		File f = new File(inDir);
		String[] listDirs 	= f.list();
		Arrays.sort(listDirs);
		int nrCategories	= 0;//listDirs.length; // equal to nr different FOLDERS here
		for (int aa=0; aa<listDirs.length; aa++) {
			File fdummy = new File(inDir+File.separator+listDirs[aa]);
			if (fdummy.isDirectory()) {
				nrCategories++;
			}
			//System.out.println("\n" + fdummy.getAbsolutePath()+" -> " + nrCategories);
		}

		String[] listFiles =  new File(f.getAbsolutePath()+File.separator+listDirs[0]).list();
		int nrImgsPerCategory = 0;
		for (int aa=0; aa<listFiles.length; aa++)
			if (extension(listFiles[aa]).equalsIgnoreCase("TIF"))
				nrImgsPerCategory++;

		int[][] performanceTable = new int[nrImgsPerCategory][3*nrCategories];

		System.out.println("NR IMGS FROM: "+f.getAbsolutePath()+File.separator+listDirs[0]);

		/***** initialize Detector *****/
		//Detector det = new Detector();    			// default
		Detector det = new Detector(minCosAngle, minFuzzyScore, scatterDistSquared, wStdRatioToD);    // with params

		/***** store legend ******/
		for (int ii=0; ii<listDirs.length; ii++) {
			if (new File(inDir+File.separator+listDirs[ii]).isDirectory())
				logWriter.print(listDirs[ii]+" ");
		}
		logWriter.println();

		/***** store parameters used ******/
		parWriter.println("MIN_COS_ANG\t\t= " 			+ Detector.MIN_COS_ANG);
		parWriter.println("MIN_FUZZY_SCORE\t\t= " 		+ Detector.MIN_FUZZY_SCORE);
		parWriter.println("SCATTER_DISR_SQUARED\t\t= " 	+ Detector.SCATTER_D2);
		parWriter.println("NR.CATEGORIES\t\t= " 		+ nrCategories);
		parWriter.println("NR.IMGS PER CATEGORY\t\t= " 	+ nrImgsPerCategory);
		parWriter.println("D\t\t= " 					+ neuronD);
		parWriter.println("iDiff\t\t= " 				+ iDiff);
		parWriter.println("MIN_SIZE\t\t= " 				+ MIN_SIZE);
		parWriter.println("LOCATION_TOLERANCE_SCALE\t\t= " + LOCATION_TOLERANCE_SCALE*neuronD);

		int imgIdx,catIdx=0;  // for the table with results

		for (int ii=0; ii<listDirs.length; ii++) { // ii -> category

			imgIdx=0; // for the table with results

			System.out.println("\n" + listDirs[ii] + "\n");
			//logWriter.println("#CATEGORY " + catIdx + "\t\t= " + listDirs[ii]);

			File ff = new File(inDir+File.separator+listDirs[ii]);

			if (!ff.isDirectory()) {System.out.println(ff.getAbsolutePath() + " was not directory skipping... catIdx: "+catIdx+", imgIdx:"+imgIdx); continue;}

			String[] listDirsDirs = ff.list();

			for (int iii=0; iii<listDirsDirs.length; iii++) { // iii -> file index

				String fileName 	= listDirsDirs[iii];
				String fileNameExt 	= extension(fileName);
				String fileNameTitle= title(fileName);

				if (fileNameExt.equalsIgnoreCase("TIF")) {

					// loop the rest of the folder to find corresponding .csv
					for (int iiii = 0; iiii < listDirsDirs.length; iiii++) {

						if (iiii!=iii) {

							String csvFileName = listDirsDirs[iiii];
							String csvFileNameExt = extension(csvFileName);
							String csvFileNameTitle = title(csvFileName);

							if (csvFileNameTitle.equalsIgnoreCase(fileNameTitle) && csvFileNameExt.equalsIgnoreCase("CSV")) {

								// load image
								File imgFile = new File(inDir+File.separator+listDirs[ii]+File.separator+listDirsDirs[iii]);
								File csvFile = new File(inDir+File.separator+listDirs[ii]+File.separator+listDirsDirs[iiii]);

								System.out.println(imgFile);
								System.out.println(csvFile);

								ImagePlus imp = new ImagePlus(imgFile.getAbsolutePath());
								// load annotation
								AnalyzeCSV anCsv = new AnalyzeCSV(csvFile.getAbsolutePath());
								int ln = anCsv.getLinesNr();
								float[][] annot = anCsv.readLn(2);

								/*
								bif. detection
					 			*/

								ArrayList<ArrayList<int[]>> detRegions = det.run(imp, neuronD, iDiff);  // connected regions

								Vector<float[]> detList = det.formDetectionList(detRegions, MIN_SIZE);  // list of discs

								int[] labels = det.clustering(detList, LOCATION_TOLERANCE_SCALE*neuronD);  								// prune detection list, list of discs

								Vector<float[]> detLstPruned = Detector.groupClusters(labels, detList); // group clusters using labels

								System.out.println(detList.size() + "detections, " + detLstPruned.size() + "pruned detections");

								//Overlay overlayDet = det.formPointOverlay(detLstPruned, Color.YELLOW);
								//imp.setOverlay(overlayDet);

								// allocate evaluation
								float[] eval = new float[3];
								Overlay overlayRes = det.classificationPerformanceOverlay(detLstPruned, annot, LOCATION_TOLERANCE_SCALE*neuronD, eval);
								imp.setOverlay(overlayRes);

								// save the image with detections/detection results
								FileSaver fs = new FileSaver(imp);
								GeneratorDemo.createDirs(outDir + listDirs[ii]);
								fs.saveAsTiff(outDir+listDirs[ii]+File.separator+"DET_"+fileName);

								// evaluate
								//float[] eval = Detector.classificationPerformance(detLstPruned, annot, LOCATION_TOLERANCE_SCALE*neuronD);
								System.out.println("tp: " + eval[0] + ", fp: " + eval[1] + ", fn: " + eval[2]);

								// save to the table
								System.out.println("save to table "+performanceTable.length+" x "+performanceTable[0].length+" catIdx: "+catIdx+", imgIdx:"+imgIdx);
								performanceTable[imgIdx][3*catIdx+0] = (int) eval[0];
								performanceTable[imgIdx][3*catIdx+1] = (int) eval[1];
								performanceTable[imgIdx][3*catIdx+2] = (int) eval[2];
								imgIdx++; // for the table with results

							}

						}

					}

				}

			}

			catIdx++; // for the table with results

		}

		// export table to res
		for (int i=0; i<performanceTable.length; i++) {
			for (int j=0; j<performanceTable[0].length; j++) {
				outWriter.print(String.format("%3d  ", performanceTable[i][j]));
			}
			outWriter.println();
		}

	}

	private static String extension(String in)
	{
		return in.substring(in.indexOf(".")+1);
	}

	private static String title(String in)
	{
		return in.substring(0, in.indexOf("."));
	}





}
