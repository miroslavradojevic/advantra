package detection;

import aux.Interpolator;
import aux.Tools;
import conn.Find_Connected_Regions;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.*;
import ij.measure.ResultsTable;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 8/7/13
 * Time: 9:45 AM
 */
public class JunctionDet_v11 implements PlugInFilter, MouseListener, MouseMotionListener, KeyListener {

    ImagePlus 		imp;
    ImageCanvas 	canvas;
	String 			inimgTitle;

	FloatProcessor 	fuzzyScores;

	String 			timestamp = "";
	String			logFile 	= "";
	String 			profileFile = "JunctionDet_v11.log";
	PrintWriter 	logWriter;

    float 	scale;
    double 	D;
	float 	iDiff;
	float 	r;
	float   resolDeg;

	ArrayList<ArrayList<float[]>>   						profiles;
    ArrayList<int[]>   										centersXY;
	ArrayList<ArrayList<ArrayList<Float>>>					anglesRad;
	ArrayList<ArrayList<ArrayList<float[]>>>				peaksXY;
    ArrayList<ArrayList<ArrayList<ArrayList<float[]>>>>   	topologyXY;
	ArrayList<float[][]> 									theta;

	PlotWindow pwP_A;//, pwI_A;

    Overlay     currOvl = new Overlay();  // mouse click overlay points
    Overlay 	detectionOverlay = new Overlay();
	ResultsTable resTab = new ResultsTable();

    // used for export variables (used by file export (buttonPressed) and mouse click plots (mouseCLicked))
    float[] ang_A, ang_B;  // they cover [0, 360) but with different sampling rate
    float[] i_A, i_B;
    float[] iTh_A, iTh_B;
    float[] peakAng_A, peakAng_B;
    float[] peakIdx_A, peakIdx_B;
    float[] profile_A = null, profile_B = null;

	Fuzzy fz;

    int CPU_NR;

    private static float 	Deg2Rad = (float) (Math.PI/180f);
    private static float 	Rad2Deg = (float) (180f/Math.PI);
    private static float 	MIN_COS_ANG = .6f;
	private static int 	 	MIN_SIZE = 3;
	private static float 	MIN_FUZZY_SCORE = .6f;
	private static float 	SCATTER_D2 = 5;
	private static int		PROFILE_STD_RATIO_TO_D = 3; // std = neuronDIam/PROFILE_W_STD
	private static boolean useMax = true; // to estimate theta (inupt to fuzzy)
	private static boolean showAll   = true;

	public int setup(String s, ImagePlus imagePlus)
	{
		if(imagePlus==null) return DONE;
		imp = Tools.convertToFloatImage(imagePlus);
		imp.setTitle("inimg");
		canvas = imagePlus.getCanvas();
		inimgTitle = imagePlus.getTitle();
		return DOES_8G+DOES_32+NO_CHANGES;
	}

	private Color getColor(int clusterIdx)
	{
		if (clusterIdx<=0) return Color.RED;
		else if (clusterIdx==1) return Color.GREEN;
		else if (clusterIdx==2) return Color.BLUE;
		else return Color.MAGENTA;
	}

	public void run(ImageProcessor imageProcessor)
    {
		/***********************************************************/
		DateFormat dateFormat = new SimpleDateFormat("yyyy_MM_dd__HH_mm_ss");
		timestamp =  dateFormat.format(new Date());
		logFile = "JunctionDet_v11_" + timestamp + ".log";

		int 	totalJobs; // used in paralelization log
		long 	t1, t2;
		/***********************************************************/


		/***********************************************************/
		D       			=  			Prefs.get("advantra.critpoint.D", 		3);
		scale               = (float) 	Prefs.get("advantra.critpoint.scale", 	1.5);
		iDiff				= (float) 	Prefs.get("advantra.critpoint.iDiff", 	5);

		GenericDialog gd = new GenericDialog("JUNCTION DET.");
		gd.addNumericField("neuronD ",  D, 		0, 5, " pix.");
		gd.addNumericField("scale ",  	scale, 	1, 5, " x(neuronD)");
		gd.addNumericField("iDiff ",  	iDiff, 	1, 5, " ");

		gd.showDialog();
		if (gd.wasCanceled()) return;

		D       =  gd.getNextNumber();
		Prefs.set("advantra.critpoint.D",   D);

		scale               = (float) gd.getNextNumber();
		Prefs.set("advantra.critpoint.scale",   scale);

		iDiff               = (float) gd.getNextNumber();
		Prefs.set("advantra.critpoint.iDiff",   iDiff);

		fz = new Fuzzy(iDiff);

		CPU_NR = Runtime.getRuntime().availableProcessors();

		r = (float) (scale * D);
		resolDeg = Profiler.getResolDeg(scale);

		// empty the log file
		try {
			PrintWriter writer = new PrintWriter(logFile);
			writer.print("");
			writer.close();
		}
		catch (FileNotFoundException ex) {}
		/***********************************************************/

        /***********************************************************/
		//ImageProcessor nness = MyHessian.nness(imp, new double[]{D/2, D, 2*D});
		//new ImagePlus("nness", nness).show();
		/***********************************************************/


		/***********************************************************/
		int neighbourhoodR = (int) Math.ceil(1.5f*scale*D);
		IJ.log("extracting background... neigh. radius = "+neighbourhoodR + " pixels, iDiff = "+iDiff);
		t1 = System.currentTimeMillis();
		Masker.loadTemplate(imp.getProcessor(), neighbourhoodR, iDiff);
		totalJobs = imp.getHeight()*imp.getWidth();
		Masker masker_jobs[] = new Masker[CPU_NR];
		for (int i = 0; i < masker_jobs.length; i++) {
			masker_jobs[i] = new Masker(i*totalJobs/CPU_NR,  (i+1)*totalJobs/CPU_NR);
			masker_jobs[i].start();
		}
		for (int i = 0; i < masker_jobs.length; i++) {
			try {
				masker_jobs[i].join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		t2 = System.currentTimeMillis();
		ImagePlus inmask = new ImagePlus("inmask", Masker.mask);
		inmask.show();
			//imp.show();
			//IJ.selectWindow("inimg");
        	//IJ.run("Add Image...", "image=inmask x="+0+" y="+0+" opacity=30");  detectionOverlay = imp.getOverlay();
			//inmask.close();
		int nr = Masker.getMaskLocationsTotal();
		float ratio = 100*(float)nr/(totalJobs);
		IJ.log("done "+((t2-t1)/1000f)+" sec. total "+nr+" given by Masker ("+ratio+" % kept) "+Profiler.getProfilerLocationsTotal());
		GenericDialog cont = new GenericDialog("Continue?");
		cont.addMessage("Is background eliminated? ("+ratio+" % kept)");
		cont.enableYesNoCancel();
		cont.showDialog();
		if (cont.wasCanceled()||!cont.wasOKed()) {inmask.close(); return; }
		inmask.close();
		/***********************************************************/


		/***********************************************************/
		IJ.log("calculating profiles... ");
		t1 = System.currentTimeMillis();
		// set the template first
		Profiler.loadTemplate(imp.getProcessor(), Masker.mask, D, scale, PROFILE_STD_RATIO_TO_D, false);
		totalJobs = Profiler.offsets.size();
		Profiler[] profiler_jobs;

		profiler_jobs = new Profiler[CPU_NR];
		for (int i = 0; i < profiler_jobs.length; i++) {
			profiler_jobs[i] = new Profiler(i*totalJobs/CPU_NR,  (i+1)*totalJobs/CPU_NR);
			profiler_jobs[i].start();
		}
		for (int i = 0; i < profiler_jobs.length; i++) {
			try {
				profiler_jobs[i].join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		profiles     	= new ArrayList<ArrayList<float[]>>(nr);
		addProfilerList(profiles, Profiler.profiles);
		t2 = System.currentTimeMillis();
		IJ.log("done "+((t2-t1)/1000f)+" sec.");
		//Profiler.getLocation2IndexMapping().show(); don't show index mapping
		/***********************************************************/



		/***********************************************************/
		IJ.log("calculating peaks for selected location profiles...");
		t1 = System.currentTimeMillis();

		Analyzer.loadProfiles(profiles);
		totalJobs = Analyzer.profiles.size();

		Analyzer analyzer_jobs[] = new Analyzer[CPU_NR];
		for (int i = 0; i < analyzer_jobs.length; i++) {
			analyzer_jobs[i] = new Analyzer(i*totalJobs/CPU_NR,  (i+1)*totalJobs/CPU_NR);
			analyzer_jobs[i].start();
		}
		for (int i = 0; i < analyzer_jobs.length; i++) {
			try {
				analyzer_jobs[i].join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}

		t2 = System.currentTimeMillis();
		IJ.log("done extracting peaks "+((t2-t1)/1000f)+" sec.");
		/***********************************************************/



		/***********************************************************/
		IJ.log("detection...");
		fuzzyScores = new FloatProcessor(imp.getWidth(), imp.getHeight());
		t1 = System.currentTimeMillis();

		centersXY 	= new ArrayList<int[]>(nr); // from Profiler
		for (int qq=0; qq<nr; qq++) centersXY.add(qq, new int[]{Profiler.locations[qq][0], Profiler.locations[qq][1]});

		anglesRad = Analyzer.exportPeakIdx();

//		// DBG
//		for (int q=0; q<Analyzer.peakIdx.size(); q++) {
//			for (int q1=0; q1<Analyzer.peakIdx.get(q).size(); q1++) {
//				for (int q2=0; q2<Analyzer.peakIdx.get(q).get(q1).size(); q2++) {
//					System.out.print(" (" + Analyzer.peakIdx.get(q).get(q1).get(q2) + ", " +  ") "); // anglesRad.get(q).get(q1).get(q2)[1] +
//				}
//				System.out.print(" | ");
//			}
//			System.out.println();
//		}

		for (int ii=0; ii<anglesRad.size(); ii++) {
			for (int jj=0; jj<anglesRad.get(ii).size(); jj++) {
				for (int kk=0; kk<anglesRad.get(ii).get(jj).size(); kk++) {
					float newValue = anglesRad.get(ii).get(jj).get(kk) * resolDeg * Deg2Rad;
					anglesRad.get(ii).get(jj).set(kk, newValue);
				}
			}
		}

//		// DBG
//		for (int q=0; q<anglesRad.size(); q++) {
//			for (int q1=0; q1<anglesRad.get(q).size(); q1++) {
//				for (int q2=0; q2<anglesRad.get(q).get(q1).size(); q2++) {
//					System.out.print(" (" + anglesRad.get(q).get(q1).get(q2) + ", " +  ") "); // anglesRad.get(q).get(q1).get(q2)[1] +
//				}
//				System.out.print(" | ");
//			}
//			System.out.println();
//		}

		peaksXY		= new ArrayList<ArrayList<ArrayList<float[]>>>(nr);
		for (int ii=0; ii<nr; ii++) {
			ArrayList<ArrayList<float[]>> peaksXY1 = new ArrayList<ArrayList<float[]>>(anglesRad.get(ii).size());
			for (int jj=0; jj<anglesRad.get(ii).size(); jj++) {
				ArrayList<float[]> peaksXY2 = new ArrayList<float[]>(anglesRad.get(ii).get(jj).size());
				for (int kk=0; kk<anglesRad.get(ii).get(jj).size(); kk++){
					float[] pkXY = new float[2];
					// use anglesDeg.get(ii).get(jj).get(kk), centersXY.get(ii)
					float currAngle = anglesRad.get(ii).get(jj).get(kk);
					pkXY[0] = (float) (centersXY.get(ii)[0] + r * Math.cos( currAngle ));
					pkXY[1] = (float) (centersXY.get(ii)[1] - r * Math.sin( currAngle ));
//					// DBG
//					System.out.print(" (" + pkXY[0] + ", " + pkXY[1] + "):(" + centersXY.get(ii)[0] + ", " + centersXY.get(ii)[1] + "):(" + r + " )");
					peaksXY2.add(kk, pkXY);
				}
				peaksXY1.add(jj, peaksXY2);
			}
			peaksXY.add(ii, peaksXY1);
//			// DBG
//			System.out.println();
		}

//		// DBG
//		for (int q=0; q<peaksXY.size(); q++) {
//			for (int q1=0; q1<peaksXY.get(q).size(); q1++) {
//				for (int q2=0; q2<peaksXY.get(q).get(q1).size(); q2++) {
//					System.out.print(" (" + peaksXY.get(q).get(q1).get(q2)[0] + ", " + peaksXY.get(q).get(q1).get(q2)[1] + ") "); //
//				}
//				System.out.print(" | ");
//			}
//			System.out.println();
//		}

		// initialize detection log file
		try {
			logWriter = new PrintWriter(new BufferedWriter(new FileWriter(logFile, true)));
			logWriter.println("BIFURCATION DETECTION LOG\n"+timestamp); // dateFormat.format(date)
		} catch (IOException e) {}


		// extract locations and fuzzy scores

		topologyXY	= new ArrayList<ArrayList<ArrayList<ArrayList<float[]>>>>(nr);
		theta		= new ArrayList<float[][]>(nr);   							// parameters describing intensity
		int[][] currentPeakPoint  	= new int[4][2];
		float[]	currentPeakVals		= new float[4];
		float[] nextPeakPoint		= new float[2];

		for (int ii=0; ii<nr; ii++) { // peaksXY.size()

			int[] currentCenter = centersXY.get(ii);
			float iBackgC = Masker.backgr.getf(currentCenter[0], currentCenter[1]);


			logWriter.println(IJ.d2s(ii,0)+": ("+IJ.d2s(currentCenter[0], 0)+","+IJ.d2s(currentCenter[1], 0)+") -> ");

			// check how many clusters there could be
			ArrayList<ArrayList<ArrayList<float[]>>> locationTopology = new ArrayList<ArrayList<ArrayList<float[]>>>(4);
			// 4 clusters max. : clusters(1-4) X thread points (1-4) X points in line (1-2)

			float[][] thetaArray        = new float[4][2]; // with parameters
			boolean[] scattered = new boolean[peaksXY.get(ii).get(0).size()]; // per location, in amount of peaks detected
			for (int chkPeak=0; chkPeak<peaksXY.get(ii).get(0).size(); chkPeak++) {

				float[] debug = new float[2];

				// 4 points around the peak
				currentPeakPoint[0][0] 	= (int) Math.floor(peaksXY.get(ii).get(0).get(chkPeak)[0]);
				currentPeakPoint[0][1] 	= (int) Math.floor(peaksXY.get(ii).get(0).get(chkPeak)[1]);
				//currentPeakVals[0] 		=       imp.getProcessor().getf(currentPeakPoint[0][0], currentPeakPoint[0][1]);
				// TODO: substitute with line averages/median, maybe one line is enough instead of four!!!
				currentPeakVals[0] 		=  avgAlongLine(currentCenter[0], currentCenter[1], currentPeakPoint[0][0], currentPeakPoint[0][1], imp.getProcessor());
				//logWriter.println("what happened?  peak " + chkPeak + " : " + currentPeakVals[0]);
				//logWriter.println("start (" + currentCenter[0] + "," + currentCenter[1] + ") -> " + currentPeakPoint[0][0] + "," + currentPeakPoint[0][1] + ")");
				//logWriter.println("debug: "+Arrays.toString(debug));
				currentPeakPoint[1][0] = (int) Math.floor(peaksXY.get(ii).get(0).get(chkPeak)[0]);
				currentPeakPoint[1][1] = (int) Math.ceil (peaksXY.get(ii).get(0).get(chkPeak)[1]);
				//currentPeakVals[1] 		=       imp.getProcessor().getf(currentPeakPoint[1][0], currentPeakPoint[1][1]);
				// TODO
				currentPeakVals[1]		=  avgAlongLine(currentCenter[0], currentCenter[1], currentPeakPoint[1][0], currentPeakPoint[1][1], imp.getProcessor());

				currentPeakPoint[2][0] = (int) Math.ceil (peaksXY.get(ii).get(0).get(chkPeak)[0]);
				currentPeakPoint[2][1] = (int) Math.floor(peaksXY.get(ii).get(0).get(chkPeak)[1]);
				//currentPeakVals[2] 		=       imp.getProcessor().getf(currentPeakPoint[2][0], currentPeakPoint[2][1]);
				// TODO
				currentPeakVals[2]		=  avgAlongLine(currentCenter[0], currentCenter[1], currentPeakPoint[2][0], currentPeakPoint[2][1], imp.getProcessor());

				currentPeakPoint[3][0] = (int) Math.ceil (peaksXY.get(ii).get(0).get(chkPeak)[0]);
				currentPeakPoint[3][1] = (int) Math.ceil (peaksXY.get(ii).get(0).get(chkPeak)[1]);
				//currentPeakVals[3] 		=       imp.getProcessor().getf(currentPeakPoint[3][0], currentPeakPoint[3][1]);
				// TODO
				currentPeakVals[3]		=  avgAlongLine(currentCenter[0], currentCenter[1], currentPeakPoint[3][0], currentPeakPoint[3][1], imp.getProcessor());

				// thetaArray[chkPeak][0]
				/****/
				if (useMax)
					thetaArray[chkPeak][0] = maxArray(currentPeakVals);// is max here
				else
					thetaArray[chkPeak][0] = _NthLowest(currentPeakVals, 3); // N lowest out of 4

				//logWriter.println(" max " + thetaArray[chkPeak][0]);

				ArrayList<ArrayList<float[]>> locationCluster = new ArrayList<ArrayList<float[]>>(4); // thread points (1-4) X points in line (1-2)

				for (int kk=0; kk<4; kk++) {  // check every point

					int findInList = Profiler.indexInList(currentPeakPoint[kk]); // check if it is in location list

					if (findInList >= 0) { // was in the mask (got it's peaks)

						ArrayList<float[]> locationThread = new ArrayList<float[]>();

						// point 1
						locationThread.add(0, new float[]{currentPeakPoint[kk][0], currentPeakPoint[kk][1]});
						// point 2
						if (nextPeakLoc(peaksXY.get(findInList).get(0), currentCenter, currentPeakPoint[kk], nextPeakPoint)) { // peaksXY.get(ii) would be for current

							locationThread.add(1, new float[]{nextPeakPoint[0], nextPeakPoint[1]});// nextPeakPoint.clone()

							// thetaArray[chkPeak][1] is max here  // median will be used in case there was more than 2, later re-assigned
							//debug = new float[2];
							thetaArray[chkPeak][1] = Math.max(thetaArray[chkPeak][1],
																	 //Interpolator.interpolateAt(nextPeakPoint[0], nextPeakPoint[1], (FloatProcessor) imp.getProcessor()));
																	 avgAlongLine(
																						 currentPeakPoint[kk][0],
																						 currentPeakPoint[kk][1],
																						 nextPeakPoint[0],
																						 nextPeakPoint[1],
																						 imp.getProcessor()
																						 )
							);

						}

						locationCluster.add(locationThread);

					}

				}

				locationTopology.add(locationCluster);

				// check the topology for this peak, how scattered they are at second ring, use locationCluster
				// check scattering (using positons)
				// put constrain so that the last one has to have more than 1 following point (out of 4)

					if (locationCluster.size()>=2) {

						ArrayList<Float> 	values2  	= new ArrayList<Float>(4);
						ArrayList<float[]> 	locs2  		= new ArrayList<float[]>(4);

						for (int aa=0; aa<locationCluster.size(); aa++) {
							if (locationCluster.get(aa).size()>1) {
								// there is 2nd point, add it
								float[] loc2 = new float[2];
								loc2[0] = locationCluster.get(aa).get(1)[0];
								loc2[1] = locationCluster.get(aa).get(1)[1];
								values2.add(Interpolator.interpolateAt(loc2[0], loc2[1], (FloatProcessor) imp.getProcessor()));
								locs2.add(loc2);
							}
						}

						/****/
						// USE values2 correct  thetaArray[chkPeak][1] - use median instead of previously obtained max
						if (values2.size()>=2 && !useMax) {
							if (values2.size()==2) thetaArray[chkPeak][1] = _NthLowest(values2, 1);
							else if (values2.size()==3) thetaArray[chkPeak][1] = _NthLowest(values2, 2);
							else if (values2.size()==4) thetaArray[chkPeak][1] = _NthLowest(values2, 2);
							else {}
						}

						// assign as scattered if there was only one point at last stage
						if (locs2.size()==1) {
							scattered[chkPeak] = true; // will give a negative vote (exclude it) for this one when it comes to detection
						}

						// USE locs2 check if it is scattered - find inter distance
						if (locs2.size()==2) {
							float interD = (float) (Math.pow(locs2.get(0)[0]-locs2.get(1)[0], 2) + Math.pow(locs2.get(0)[1] - locs2.get(1)[1], 2)); // 0 1
							if (interD>SCATTER_D2) {
								scattered[chkPeak] = true; // will stop fuzzification/detection later
							}
						}
						else if (locs2.size()==3){
							float[] interD = new float[3];
							interD[0] = (float) (Math.pow(locs2.get(0)[0]-locs2.get(1)[0], 2) + Math.pow(locs2.get(0)[1] - locs2.get(1)[1], 2)); // 0 1
							interD[1] = (float) (Math.pow(locs2.get(1)[0]-locs2.get(2)[0], 2) + Math.pow(locs2.get(1)[1] - locs2.get(2)[1], 2)); // 1 2
							interD[2] = (float) (Math.pow(locs2.get(0)[0]-locs2.get(2)[0], 2) + Math.pow(locs2.get(0)[1] - locs2.get(2)[1], 2)); // 0 2
							if (_NthLowest(interD, 1)>SCATTER_D2) {
								scattered[chkPeak] = true;  // will stop fuzzification later
							}

						}
						else if (locs2.size()==4) {
							float[] interD = new float[6];
							interD[0] = (float) (Math.pow(locs2.get(0)[0]-locs2.get(1)[0], 2) + Math.pow(locs2.get(0)[1] - locs2.get(1)[1], 2)); // 0 1
							interD[1] = (float) (Math.pow(locs2.get(0)[0]-locs2.get(2)[0], 2) + Math.pow(locs2.get(0)[1] - locs2.get(2)[1], 2)); // 0 2
							interD[2] = (float) (Math.pow(locs2.get(0)[0]-locs2.get(3)[0], 2) + Math.pow(locs2.get(0)[1] - locs2.get(3)[1], 2)); // 0 3
							interD[3] = (float) (Math.pow(locs2.get(1)[0]-locs2.get(2)[0], 2) + Math.pow(locs2.get(1)[1] - locs2.get(2)[1], 2)); // 1 2
							interD[4] = (float) (Math.pow(locs2.get(1)[0]-locs2.get(3)[0], 2) + Math.pow(locs2.get(1)[1] - locs2.get(3)[1], 2)); // 1 3
							interD[5] = (float) (Math.pow(locs2.get(2)[0]-locs2.get(3)[0], 2) + Math.pow(locs2.get(2)[1] - locs2.get(3)[1], 2)); // 2 3
							if (_NthLowest(interD, 3)>SCATTER_D2) {
								scattered[chkPeak] = true;  // will stop fuzzification later
							}
						}
						else {
							// nothing
						}



					}

			} // thetaArray & locationTopology are done

			logWriter.println("how scattered? "+Arrays.toString(scattered));
			// at least 3 with false scattered[]


			// FUZZY DETECTION (use peaksXY.get(ii).get(0) and thetaArray)
			int nrPeaks = peaksXY.get(ii).get(0).size();
			if (nrPeaks>=3 ) {

				// count amount of "falses" in scattered[nrPeaks]
				int cntNonScattered = 0;
				for (int ee=0; ee<nrPeaks; ee++) if (scattered[ee]==false) cntNonScattered++;

				if (cntNonScattered >=3) { // at least 3 with false in scattered[nrPeaks] is OK to proceed


				// makes sense to carry out fuzzy decision
				//read from peaksXY.get(ii).get(0) and thetaArray ---> form gg[][]
				float[][] gg = new float[3][3]; // final input for Fuzzy 6 values + 3 backgrounds
				for (int i=0; i<nrPeaks; i++) {  // peaksXY.get(ii).get(0).size()

					float iBackgP = Interpolator.interpolateAt(peaksXY.get(ii).get(0).get(i)[0], peaksXY.get(ii).get(0).get(i)[1], Masker.backgr);

					if (i<=2) { // take first three automatically  !!!! what if one is scattered
						gg[i][0] = thetaArray[i][0]; 	// theta 1
						gg[i][1] = thetaArray[i][1]; 	// theta 2
						gg[i][2] = iBackgP; 			// bg
					}
					else {

						boolean overW = false;

						// first overwrite the one that was scattered in first 3  - there can be one such at max.
						int scatteredIdx = -1;
						if (scattered[0]==true) {overW = true; scatteredIdx = 0; }
						if (scattered[1]==true) {overW = true; scatteredIdx = 1; }
						if (scattered[2]==true) {overW = true; scatteredIdx = 2; }

						if (overW) {

							logWriter.println("substituting cluster idx: "+scatteredIdx+" it was scattered");

							gg[scatteredIdx][0] = thetaArray[i][0];
							gg[scatteredIdx][1] = thetaArray[i][1];
							gg[scatteredIdx][2] = iBackgP;

						}
						else {

							// overwrite lowest if necessary
							int lowestIdx = -1;
							if ( Math.min(gg[0][0], gg[0][1]) <= Math.min( Math.min(gg[1][0], gg[1][1]) , Math.min(gg[2][0], gg[2][1]) ) ) lowestIdx = 0;
							if ( Math.min(gg[1][0], gg[1][1]) <= Math.min( Math.min(gg[0][0], gg[0][1]) , Math.min(gg[2][0], gg[2][1]) ) ) lowestIdx = 1;
							if ( Math.min(gg[2][0], gg[2][1]) <= Math.min( Math.min(gg[0][0], gg[0][1]) , Math.min(gg[1][0], gg[1][1]) ) ) lowestIdx = 2;

							logWriter.println("substituting cluster idx: "+lowestIdx+" it was smallest");

							if ( Math.min( thetaArray[i][0] , thetaArray[i][1] ) > Math.min( gg[lowestIdx][0] , gg[lowestIdx][1] ) ) {
								gg[lowestIdx][0] = thetaArray[i][0];
								gg[lowestIdx][1] = thetaArray[i][1];
								gg[lowestIdx][2] = iBackgP;
							}

						}

					}
					//float iBackg = Masker.backgr.getf(centersXY.get(where)[0], centersXY.get(where)[1]);
				}

				float centralValue = imp.getProcessor().getf(currentCenter[0], currentCenter[1]);

				float fuzzyValue =
						fz.bifurcationess(centralValue-iBackgC,
												 gg[0][0] - gg[0][2], gg[0][1] - gg[0][2],
												 gg[1][0] - gg[1][2], gg[1][1] - gg[1][2],
												 gg[2][0] - gg[2][2], gg[2][1] - gg[2][2],
												 false);

				fuzzyScores.setf(centersXY.get(ii)[0], centersXY.get(ii)[1], fuzzyValue);

				// log values
				//logWriter.println("fuzzy input: iDiff is "+iDiff);
				logWriter.println(	"central fuzzy pair -> "+(centralValue - iBackgC)+" : ("+centralValue+","+iBackgC+")"); // :"+centralValue+" , " + iBackgC + "
				for (int ee=0; ee<gg.length; ee++) {
					logWriter.println(	(ee+1)+	". fuzzy pair -> "+(gg[ee][0] - gg[ee][2])+" ("+gg[ee][0]+","+gg[ee][2]+") : "+(gg[ee][1] - gg[ee][2])+" ("+gg[ee][1]+","+gg[ee][2]+")");
					// :"+gg[ee][0]+" , "+gg[ee][1]+" , " + gg[ee][2]+"
				}
				logWriter.println("fuzzy sc: "+fuzzyValue);

				}
			}
			else {
				// no need - peak detection said it is not!
			}

/*			// log locationTopology
			logWriter.println("location topology:");
			for (int ee=0; ee<locationTopology.size(); ee++) {
				logWriter.println("cluster "+ee+" : ");
				for (int ff=0; ff<locationTopology.get(ee).size(); ff++) {
					logWriter.print("thread "+ff+" : ");
					for (int gg=0; gg<locationTopology.get(ee).get(ff).size(); gg++) {
						logWriter.print(locationTopology.get(ee).get(ff).get(gg)[0]+","+locationTopology.get(ee).get(ff).get(gg)[1]+"  ,  ");
					}
					logWriter.print(" | ");
				}
				logWriter.println();
			}*/

			// updates
			topologyXY.add(locationTopology);
			theta.add(thetaArray);

		}
		t2 = System.currentTimeMillis();
		IJ.log("done "+((t2-t1)/1000f)+" sec.");IJ.log("done");
		/***********************************************************/

		logWriter.close();
		IJ.log("done, log saved to "+new File(logFile).getAbsolutePath());

		//ImagePlus fuzzyScoresImagePlus = new ImagePlus("fuzzyScores", fuzzyScores);
		//fuzzyScoresImagePlus.show();

		/***********************************************************/
		t1 = System.currentTimeMillis();
		ByteProcessor score = new ByteProcessor(imp.getWidth(), imp.getHeight());
		for (int ii=0; ii<imp.getWidth()*imp.getHeight(); ii++) if (fuzzyScores.getf(ii) >= MIN_FUZZY_SCORE) score.set(ii, 255);
		Find_Connected_Regions conn_reg = new Find_Connected_Regions(new ImagePlus("", score), true);
		conn_reg.run("");
		//int nr_regions = conn_reg.getNrConnectedRegions();

		//ImagePlus imageLabels = conn_reg.showLabels();
		//imageLabels.show();

		detectionOverlay = formPointOverlay(conn_reg.getConnectedRegions(), MIN_SIZE);
		resTab = formResultsTable(conn_reg.getConnectedRegions(), fuzzyScores, MIN_SIZE);
		resTab.show("BIFURCATIONS");
		int nr_regions = detectionOverlay.size();
		System.out.println("### "+detectionOverlay.size() + "detections, no prunning!");


		if (showAll) {
		// add points used to create regions of connected components
		for (int ii=0; ii<imp.getWidth()*imp.getHeight(); ii++) {
			if (score.get(ii)>=255) detectionOverlay.add(new PointRoi(ii%imp.getWidth()+.5, ii/imp.getWidth()+.5));
		}
		}


		t2 = System.currentTimeMillis();
		IJ.log(nr_regions+" regions extracted.\n" + "elapsed: "+((t2-t1)/1000f)+ " seconds.");
		/***********************************************************/


		canvas.setOverlay(detectionOverlay);
		canvas.addMouseListener(this);
		canvas.addMouseMotionListener(this);
		canvas.addKeyListener(this);

		//ImagePlus outIm1 = new ImagePlus("dets", imp.getProcessor());
//		outIm1.show();

//        ImagePlus outIm = new ImagePlus("viz", imp.getProcessor());
//        outIm.setOverlay(formOverlay(conn_reg.getConnectedRegions()));
//        outIm.show();

        IJ.selectWindow(inimgTitle);
        IJ.setTool("hand");

	}

	private float _NthLowest(float[] array, int N)
	{
		boolean[] checked = new boolean[array.length];
		float minThisLoop = Float.MAX_VALUE;
		int minIdxThisLoop;

		for (int round=0;round<N;round++) {

			// loop through each time
			minThisLoop = Float.MAX_VALUE;
			minIdxThisLoop = -1;

			for (int i=0; i<array.length; i++) {

				if (!checked[i]) {
					if (array[i]<minThisLoop) {
						minThisLoop = array[i];
						minIdxThisLoop = i;
					}
				}

			}
			checked[minIdxThisLoop] = true;

		}

		return minThisLoop;
	}

	private float _NthLowest(ArrayList<Float> array, int N)
	{
		boolean[] checked = new boolean[array.size()];
		float minThisLoop = Float.MAX_VALUE;
		int minIdxThisLoop;

		for (int round=0;round<N;round++) {

			// loop through each time
			minThisLoop = Float.MAX_VALUE;
			minIdxThisLoop = -1;

			for (int i=0; i<array.size(); i++) {

				if (!checked[i]) {
					if (array.get(i)<minThisLoop) {
						minThisLoop = array.get(i);
						minIdxThisLoop = i;
					}
				}

			}
			checked[minIdxThisLoop] = true;

		}

		return minThisLoop;
	}

	private float maxArray(float[] inArray)
	{
		float maxValue = inArray[0];
		for (int i=1; i<inArray.length; i++) {
			if (inArray[i]>maxValue)
				maxValue = inArray[i];
		}
		return maxValue;
	}

	private static boolean nextPeakLoc(ArrayList<float[]> peakList, int[] center, int[] peakOrigin, float[] nextLocExtracted)
	{

		//float[] extLoc = new float[2];
		float minAngDev = MIN_COS_ANG;
		boolean found = false;

		for (int chk=0; chk<peakList.size(); chk++) {    // loops all detected by Profiler-Analyzer

			float[] currentPoint2 = peakList.get(chk);

			//if (Profiler.indexInList(currentPoint2)>=0) { // check if it is not really background

				// if directionality ok - accept it
				// use center here: centersXY, and addLoc
				float angDev = Tools.angularDeviation(center[0], center[1], peakOrigin[0], peakOrigin[1], currentPoint2[0], currentPoint2[1]);

				if (angDev > minAngDev) {
					found = true;
					nextLocExtracted[0] = currentPoint2[0];
					nextLocExtracted[1] = currentPoint2[1];
					minAngDev = angDev;
				}

			//}

		}

		return found;

	}

	public void mouseClicked(MouseEvent e)
	{

		int offscreenX = canvas.offScreenX(e.getX());
		int offscreenY = canvas.offScreenY(e.getY());

		int where = Profiler.indexInList(offscreenX, offscreenY);

		IJ.log(""+offscreenX+", "+offscreenY+" loc idx "+where+"\n");

		if (where>=0) {

			IJ.log(topologyXY.get(where).size()+" extracted clusters here\n");
			IJ.log(peaksXY.get(where).get(0).size()+" peaks here\n");

			if (peaksXY.get(where).get(0).size()<3) {
				IJ.log("not necessary to fuzzify, not enough points...");
				return;
			}


			IJ.log("**********");
			float[][] gg = new float[3][3]; // final input for Fuzzy 6 values + 3 backgrounds
			for (int i=0; i<peaksXY.get(where).get(0).size(); i++) {

				float iBackgP = Interpolator.interpolateAt(peaksXY.get(where).get(0).get(i)[0], peaksXY.get(where).get(0).get(i)[1], Masker.backgr);
				IJ.log("theta: " + theta.get(where)[i][0] + " && " + theta.get(where)[i][1] + " , b = " + iBackgP);

				if (i<=2) {
					gg[i][0] = theta.get(where)[i][0]; 	// theta 1
					gg[i][1] = theta.get(where)[i][1]; 	// theta 2
					gg[i][2] = iBackgP; 				// bg
				}
				else {

					// substitute if necessary
					int lowestIdx = -1;
					if ( Math.min(gg[0][0], gg[0][1]) <= Math.min( Math.min(gg[1][0], gg[1][1]) , Math.min(gg[2][0], gg[2][1]) ) ) lowestIdx = 0;
					if ( Math.min(gg[1][0], gg[1][1]) <= Math.min( Math.min(gg[0][0], gg[0][1]) , Math.min(gg[2][0], gg[2][1]) ) ) lowestIdx = 1;
					if ( Math.min(gg[2][0], gg[2][1]) <= Math.min( Math.min(gg[0][0], gg[0][1]) , Math.min(gg[1][0], gg[1][1]) ) ) lowestIdx = 2;


					IJ.log("lowest idx : "+lowestIdx);

					if ( Math.min( theta.get(where)[i][0] , theta.get(where)[i][1] ) > Math.min( gg[lowestIdx][0] , gg[lowestIdx][1] ) ) {
						gg[lowestIdx][0] = theta.get(where)[i][0];
						gg[lowestIdx][1] = theta.get(where)[i][1];
						gg[lowestIdx][2] = iBackgP;
					}

				}
				//float iBackg = Masker.backgr.getf(centersXY.get(where)[0], centersXY.get(where)[1]);
			}
			IJ.log("**********");
			for (int i=0; i<gg.length; i++) IJ.log("theta -> " + gg[i][0] + " <-> " + gg[i][1] + " ----  " + gg[i][2]);
			IJ.log("**********");
			fz.bifurcationess(       imp.getProcessor().getf(offscreenX, offscreenY) - Masker.backgr.getf(offscreenX, offscreenY) ,
									 gg[0][0] - gg[0][2], gg[0][1] - gg[0][2],
									 gg[1][0] - gg[1][2], gg[1][1] - gg[1][2],
									 gg[2][0] - gg[2][2], gg[2][1] - gg[2][2],
									 false);
			IJ.log("**********");



			// plot topology
			currOvl.clear();
			currOvl = detectionOverlay.duplicate();

			for (int ii=0; ii<topologyXY.get(where).size(); ii++) { // clusters

				for (int jj=0; jj<topologyXY.get(where).get(ii).size(); jj++) { // threads

					for (int kk=0; kk<topologyXY.get(where).get(ii).get(jj).size(); kk++) {

						float x = topologyXY.get(where).get(ii).get(jj).get(kk)[0]  ;
						float y = topologyXY.get(where).get(ii).get(jj).get(kk)[1];
						PointRoi pt = new PointRoi(x+.5f, y+.5f);
						pt.setStrokeColor(getColor(ii));
						currOvl.add(pt);

					}

				}

			}

		}
		else {
			IJ.log("nothing to do");
		}


		canvas.setOverlay(currOvl);

	}

    private Overlay formOverlay(ArrayList<ArrayList<int[]>> regs)
    {

		double VERY_SMALL_POSITIVE = 0.000001;
        Overlay detections = new Overlay();

        for (int i=0; i<regs.size(); i++) {
            if (regs.get(i).size()>0) {
                float[] ellipseParams = Tools.extractEllipse(regs.get(i));

                /*
                PointRoi pt = new PointRoi(ellipseParams[1]+.5, ellipseParams[0]+.5);
                pt.setStrokeColor(Color.BLUE);
                detections.add(pt);
                */


                float A =   (float)Math.sqrt(ellipseParams[3]);
                float B =   (float)Math.sqrt(ellipseParams[2]);

                if (!(B>VERY_SMALL_POSITIVE)) {

                    // set B to 1 and A such that it covers the area
                    B = 1;
                    A = (float) (regs.get(i).size() / Math.PI);

                }
                else {

                    float k =   A/B;
                    B = (float) Math.sqrt(regs.get(i).size()/(k*Math.PI));
                    A = k*B;

                }

                float scalePlot = 2;
                detections.add(Tools.drawEllipse(ellipseParams[1]+0.5f, ellipseParams[0]+0.5f, scalePlot*A, scalePlot*B, ellipseParams[4], Color.YELLOW, 2, 50));

            }
        }

        return detections;

    }

    private Overlay formPointOverlay(ArrayList<ArrayList<int[]>> regs, int minSize)
    {

        Overlay detections = new Overlay();

        for (int i=0; i<regs.size(); i++) {
            if (regs.get(i).size()>minSize) {

                float Cx=0, Cy=0, R= (float) Math.sqrt((float)regs.get(i).size()/Math.PI);
                R = (R<1)? 1 : R ;

                for (int aa=0; aa<regs.get(i).size(); aa++) {
                    Cx += regs.get(i).get(aa)[1];
                    Cy += regs.get(i).get(aa)[0];
                }
                Cx /= regs.get(i).size();
                Cy /= regs.get(i).size();

                OvalRoi ovroi = new OvalRoi(Cx-R+.5, Cy-R+.5, 2*R, 2*R);
                ovroi.setStrokeWidth(2);
				ovroi.setStrokeColor(Color.YELLOW);
                detections.add(ovroi);

            }
        }

        return detections;

    }

	private ResultsTable formResultsTable(ArrayList<ArrayList<int[]>> regs, FloatProcessor fuzzyScores, int minSize)
	{
		ResultsTable rt = new ResultsTable();

		for (int i=0; i<regs.size(); i++) {
			if (regs.get(i).size()>minSize) {

				float Cx=0, Cy=0, avgBifurcationess=0, minBifurcationess=Float.MAX_VALUE, maxBifurcationess=Float.MIN_VALUE;

				for (int aa=0; aa<regs.get(i).size(); aa++) {
					Cx += regs.get(i).get(aa)[1];
					Cy += regs.get(i).get(aa)[0];
					float currBifurcationess = fuzzyScores.getf(regs.get(i).get(aa)[1], regs.get(i).get(aa)[0]);
					avgBifurcationess += currBifurcationess;
					if (currBifurcationess>maxBifurcationess) maxBifurcationess = currBifurcationess;
					if (currBifurcationess<minBifurcationess) minBifurcationess = currBifurcationess;
				}

				Cx /= regs.get(i).size();
				Cy /= regs.get(i).size();
				avgBifurcationess /= regs.get(i).size();

				rt.incrementCounter();
				rt.addValue("X", Cx);
				rt.addValue("Y", Cy);
				rt.addValue("size", regs.get(i).size());
				rt.addValue("avg. bness.", avgBifurcationess);
				rt.addValue("min. bness.", minBifurcationess);
				rt.addValue("max. bness.", maxBifurcationess);

			}
		}

		return rt;
	}

	private static float avgAlongLine(float ox1, float oy1, float ox2, float oy2, ImageProcessor ipSource) // , float[] debug
	{
		float 	dl 	= 0.75f;
		float 	l 	= (float) Math.sqrt(Math.pow(ox2 - ox1, 2) + Math.pow(oy2 - oy1, 2));

		float 	avg = 0;
		int 	cnt = 0;

		int w = ipSource.getWidth();
		int h = ipSource.getHeight();

		if (ox1<0 || oy1<0)
			return 0;

//		debug[0] = ox1;
//		debug[1] = oy1;

		for (float x = ox1, y = oy1, ll = 0; ll<=l; x+=dl*(ox2-ox1)/l, y+=dl*(oy2-oy1)/l, ll+=dl) { // x <= ox2 && y <= oy2

			if (x>=w-1 || y>=h-1){
				break;
			}

			avg += Interpolator.interpolateAt(x, y, (FloatProcessor) ipSource);
			cnt ++;

		}



		return (cnt>0)? avg/cnt : 0;

	}

    private void addProfilerList(ArrayList<ArrayList<float[]>> profileListToUpdate, float[][] profilesToAppend)
    // uses Profiler to update  list of profiles
	// rows correspond to locations, columns to different scale
	// one Profiler can contribute one column at the time
    {

		if (profileListToUpdate.size() == 0) {

			// store it in the list for each loc
			for (int loopLocations=0; loopLocations<profilesToAppend.length; loopLocations++) {

				// profileToAdd from Profiler.profiles
				float[] profileToAdd = new float[profilesToAppend[loopLocations].length];
				for (int k=0; k<profileToAdd.length; k++) {
					profileToAdd[k] = profilesToAppend[loopLocations][k];
				}

				ArrayList<float[]> dummy = new ArrayList<float[]>();
				dummy.add(profileToAdd);
				profileListToUpdate.add(dummy);

			}

		}
		else {

			// it exists already
			// store it in the list for each loc
			for (int loopLocations=0; loopLocations<profilesToAppend.length; loopLocations++) {

				// profileToAdd from Profiler.profiles
				float[] profileToAdd = new float[profilesToAppend[loopLocations].length];
				for (int k=0; k<profileToAdd.length; k++) {
					profileToAdd[k] = profilesToAppend[loopLocations][k];
				}

				//ArrayList<float[]> dummy = new ArrayList<float[]>();
				//dummy.add(profileToAdd);
				profileListToUpdate.get(loopLocations).add(profileToAdd);

			}

		}

    }

    public void mouseMoved(MouseEvent e) {
        //mouseClicked(e);
    }

    public void keyPressed(KeyEvent e) {

        // pressing key 'U' will export data to detection.log

        if (e.getKeyCode() == KeyEvent.VK_U) {

            if (pwP_A!=null) {

                // export to csv when clicked


                // empty the file
                PrintWriter writer = null;
                try {
                    writer = new PrintWriter(profileFile);  writer.print("");   writer.close();
                } catch (FileNotFoundException ex) {}

                try {
                    PrintWriter out;

                    out = new PrintWriter(new BufferedWriter(new FileWriter(profileFile, true)));

                    // detection A
                    out.print("angle[deg]");for (int idx=0; idx<profile_A.length; idx++)    out.print(" " + ang_A[idx]);                out.print("\n");
                    out.print("profile_A"); for (int idx=0; idx<profile_A.length; idx++)    out.print(" " + profile_A[idx]);            out.print("\n");
                    out.print("peak_A_ang");for (int idx=0; idx<peakIdx_A.length; idx++)    out.print(" " + peakAng_A[idx] * Rad2Deg);  out.print("\n");
                    out.print("peak_A_val");for (int idx=0; idx<peakIdx_A.length; idx++)    out.print(" " + Tools.interp1Darray(peakIdx_A[idx], profile_A)); out.print("\n");

                    // detection B
                    out.print("angle[deg]");for (int idx=0; idx<profile_B.length; idx++)    out.print(" " + ang_B[idx]);                out.print("\n");
                    out.print("profile_B"); for (int idx=0; idx<profile_B.length; idx++)    out.print(" " + profile_B[idx]);            out.print("\n");
                    out.print("peak_B_ang");for (int idx=0; idx<peakIdx_B.length; idx++)    out.print(" " + peakAng_B[idx] * Rad2Deg);  out.print("\n");
                    out.print("peak_B_val");for (int idx=0; idx<peakIdx_B.length; idx++)    out.print(" " + Tools.interp1Darray(peakIdx_B[idx], profile_B)); out.print("\n");

                    // intensities A
                    out.print("i_A");           for (int idx=0; idx<profile_A.length; idx++)    out.print(" " + i_A[idx]);    out.print("\n");
                    out.print("iTh_A");         for (int idx=0; idx<profile_A.length; idx++)    out.print(" " + iTh_A[idx]);  out.print("\n");
                    out.print("peak_A_Low");    for (int idx=0; idx<peakIdx_A.length; idx++)    out.print(" " + Tools.interp1Darray(peakIdx_A[idx], iTh_A)); out.print("\n");
                    out.print("peak_A_high");   for (int idx=0; idx<peakIdx_A.length; idx++)    out.print(" " + Tools.interp1Darray(peakIdx_A[idx], i_A));   out.print("\n");

                    // intensities B
                    out.print("i_B");        for (int idx=0; idx<profile_B.length; idx++)    out.print(" " + i_B[idx]);              out.print("\n");
                    out.print("iTh_B");      for (int idx=0; idx<profile_B.length; idx++)    out.print(" " + iTh_B[idx]);            out.print("\n");
                    out.print("peak_B_Low"); for (int idx=0; idx<peakIdx_B.length; idx++)    out.print(" " + Tools.interp1Darray(peakIdx_B[idx], iTh_B));  out.print("\n");
                    out.print("peak_B_high");for (int idx=0; idx<peakIdx_B.length; idx++)    out.print(" " + Tools.interp1Darray(peakIdx_B[idx], i_B)); out.print("\n");

                    out.close();

                } catch (IOException e1) {}

                IJ.log("exported to : \n" + new File(profileFile).getAbsolutePath()+ " \n ");

            }

        }

    }

    @Override
    public void keyReleased(KeyEvent e) {}
    public void mouseDragged(MouseEvent e) {}
    public void keyTyped(KeyEvent e) {}
	public void mousePressed(MouseEvent e) {}
	public void mouseReleased(MouseEvent e) {}
	public void mouseEntered(MouseEvent e) {}
	public void mouseExited(MouseEvent e) {}

}
