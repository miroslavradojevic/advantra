package detection;

import aux.Tools;
import conn.Find_Connected_Regions;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.measure.ResultsTable;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.*;
import java.io.*;
import java.util.ArrayList;
import java.util.Vector;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 9/9/13
 * Time: 11:10 AM
 */
public class Detector {

	/*
	main input is ImagePlus, output is Overlay with detections
	 */

	ImagePlus 			imp;
	FloatProcessor 		fuzzyScores;

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

	Overlay detectionOverlay;	// = new Overlay();
	ResultsTable resTab;		// = new ResultsTable();

	Fuzzy fz;

	int CPU_NR;

	private static float 	Deg2Rad = (float) (Math.PI/180f);
	private static float 	Rad2Deg = (float) (180f/Math.PI);

	public static float 	MIN_COS_ANG = .6f;
	public static float 	MIN_FUZZY_SCORE = .7f;
	public static float 	SCATTER_D2 	= 5;
	public static int		W_STD_RATIO_WRT_TO_D = 3;
	//private static int 	 	MIN_SIZE 	= 2;

	private static float 	LOCAL_NBHOOD = 1.5f;
	private static boolean 	useMax 		= true;
	//private static boolean 	showAll   	= false;

	public Detector(
						   float minCosAng,
						   float minFuzzyScore,
						   float scatterDistance2,
						   int 	 wStdRatioToD
	)
	{

		MIN_COS_ANG = minCosAng;
		MIN_FUZZY_SCORE = minFuzzyScore;
		SCATTER_D2 = scatterDistance2;
		W_STD_RATIO_WRT_TO_D = wStdRatioToD;

		scale = 1.5f;
		resolDeg = Profiler.getResolDeg(scale);



	}

	public ArrayList<ArrayList<int[]>> run(ImagePlus inimg, double neuronDiameter, float iDiff1)
	{
		// load image
		if(inimg==null) return null;
		imp = Tools.convertToFloatImage(inimg);

		D = neuronDiameter;
		iDiff = iDiff1;

		// some important parameters
		CPU_NR = Runtime.getRuntime().availableProcessors();
		r = (float) (scale * D);

		fz = new Fuzzy(iDiff);

		// paralelization parameteres
		int 	totalJobs;

		// MASK
		/***********************************************************/
		int neighbourhoodR = (int) Math.ceil(LOCAL_NBHOOD * scale * D);
		System.out.println("\nextracting background... neigh. radius = "+neighbourhoodR + ", neuron diameter = "+D+" pixels, iDiff = "+iDiff+" ...  ");
		//t1 = System.currentTimeMillis();
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
		// mask is in Masker.mask //
		//ImagePlus inmask = new ImagePlus("inmask", Masker.mask);
		//inmask.show();
		//imp.show();
		//IJ.selectWindow("inimg");
		//IJ.run("Add Image...", "image=inmask x="+0+" y="+0+" opacity=30");  detectionOverlay = imp.getOverlay();
		//inmask.close();

		int nr = Masker.getMaskLocationsTotal(); // will be used later in detection

		float ratio = 100*(float)nr/(totalJobs);
		System.out.println(nr + " locations extracted by Masker (" + ratio + " % total kept) ");
		//GenericDialog cont = new GenericDialog("Continue?");
		//cont.addMessage("Is background eliminated? ("+ratio+" % kept)");
		//cont.enableYesNoCancel();
		//cont.showDialog();
		//if (cont.wasCanceled()||!cont.wasOKed()) {inmask.close(); return; }
		//inmask.close();
		/***********************************************************/

		// PROFILES
		/***********************************************************/
		System.out.print("calculating profiles... ");
		// set the template first
		Profiler.loadTemplate(imp.getProcessor(), Masker.mask, D, scale, W_STD_RATIO_WRT_TO_D, false); // parameters for profile extraction
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
		//Profiler.getLocation2IndexMapping().show(); // don't show index mapping
		System.out.println(" done.");
		/***********************************************************/

		// PEAKS
		/***********************************************************/
		System.out.print("calculating peaks... ");
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
		System.out.println(" done.");
		/***********************************************************/

		// DETECTION (TODO : revise it - input-output, needs more code optimization, can be done in parallel as well)
		/***********************************************************/
		System.out.print("detection... ");
		fuzzyScores = new FloatProcessor(imp.getWidth(), imp.getHeight());
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

		// extract locations and fuzzy scores
		topologyXY	= new ArrayList<ArrayList<ArrayList<ArrayList<float[]>>>>(nr); // spatial topology
		theta		= new ArrayList<float[][]>(nr);   // parameters describing intensity vs. local background difference
		int[][] currentPeakPoint  	= new int[4][2];  // XY locations where profile peaks were detected
		float[]	currentPeakVals		= new float[4];   // value assigned to each peak XY location - it's weight - average along line connecting center & peak XY location
		float[] nextPeakPoint		= new float[2];   // XY locations where NEXT profile peaks were detected

		for (int ii=0; ii<nr; ii++) { // loops all selected locations

			int[] currentCenter = centersXY.get(ii);
			float iBackgC = Masker.backgr.getf(currentCenter[0], currentCenter[1]);
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
				//currentPeakVals[0] 		=  avgAlongLine(currentCenter[0], currentCenter[1], currentPeakPoint[0][0], currentPeakPoint[0][1], imp.getProcessor());
				currentPeakVals[0] 		=  medAlongLine(currentCenter[0], currentCenter[1], currentPeakPoint[0][0], currentPeakPoint[0][1], imp.getProcessor());

				currentPeakPoint[1][0] = (int) Math.floor(peaksXY.get(ii).get(0).get(chkPeak)[0]);
				currentPeakPoint[1][1] = (int) Math.ceil (peaksXY.get(ii).get(0).get(chkPeak)[1]);
				//currentPeakVals[1]		=  avgAlongLine(currentCenter[0], currentCenter[1], currentPeakPoint[1][0], currentPeakPoint[1][1], imp.getProcessor());
				currentPeakVals[1]		=  medAlongLine(currentCenter[0], currentCenter[1], currentPeakPoint[1][0], currentPeakPoint[1][1], imp.getProcessor());


				currentPeakPoint[2][0] = (int) Math.ceil (peaksXY.get(ii).get(0).get(chkPeak)[0]);
				currentPeakPoint[2][1] = (int) Math.floor(peaksXY.get(ii).get(0).get(chkPeak)[1]);
				//currentPeakVals[2]		=  avgAlongLine(currentCenter[0], currentCenter[1], currentPeakPoint[2][0], currentPeakPoint[2][1], imp.getProcessor());
				currentPeakVals[2]		=  medAlongLine(currentCenter[0], currentCenter[1], currentPeakPoint[2][0], currentPeakPoint[2][1], imp.getProcessor());

				currentPeakPoint[3][0] = (int) Math.ceil (peaksXY.get(ii).get(0).get(chkPeak)[0]);
				currentPeakPoint[3][1] = (int) Math.ceil (peaksXY.get(ii).get(0).get(chkPeak)[1]);
				//currentPeakVals[3]		=  avgAlongLine(currentCenter[0], currentCenter[1], currentPeakPoint[3][0], currentPeakPoint[3][1], imp.getProcessor());
				currentPeakVals[3]		=  medAlongLine(currentCenter[0], currentCenter[1], currentPeakPoint[3][0], currentPeakPoint[3][1], imp.getProcessor());

				// thetaArray[chkPeak][0] calculation
				thetaArray[chkPeak][0] = maxArray(currentPeakVals);			// is max here
				//	thetaArray[chkPeak][0] = _NthLowest(currentPeakVals, 3); // N lowest out of 4

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

							// thetaArray[chkPeak][1] is calculated here
							thetaArray[chkPeak][1] = Math.max(thetaArray[chkPeak][1],
							//avgAlongLine(currentPeakPoint[kk][0],currentPeakPoint[kk][1],nextPeakPoint[0],nextPeakPoint[1],imp.getProcessor())
							medAlongLine(currentPeakPoint[kk][0],currentPeakPoint[kk][1],nextPeakPoint[0],nextPeakPoint[1],imp.getProcessor())
							);

						}

						locationCluster.add(locationThread);

					}

				}

				locationTopology.add(locationCluster);

				// check the topology for this peak, how scattered they are at second ring, use locationCluster
				// in case locationCluster has >=3 samples => 6 or 3 combinations inside

				if (locationCluster.size()>=2) { // if what was added to topology has at least 2 points == there was extension

					//ArrayList<Float> 	values2  	= new ArrayList<Float>(4);
					ArrayList<float[]> 	locs2  		= new ArrayList<float[]>(4);

					for (int aa=0; aa<locationCluster.size(); aa++) {
						if (locationCluster.get(aa).size()>1) {
							// there is 2nd point, add it
							float[] loc2 = new float[2];
							loc2[0] = locationCluster.get(aa).get(1)[0];
							loc2[1] = locationCluster.get(aa).get(1)[1];
							//values2.add(Interpolator.interpolateAt(loc2[0], loc2[1], (FloatProcessor) imp.getProcessor()));
							locs2.add(loc2);
						}
					}
					/****/
//					// USE values2 correct  thetaArray[chkPeak][1] - use median instead of previously obtained max
//					if (values2.size()>=2 && !useMax) {
//						if (values2.size()==2) thetaArray[chkPeak][1] = _NthLowest(values2, 1);
//						else if (values2.size()==3) thetaArray[chkPeak][1] = _NthLowest(values2, 2);
//						else if (values2.size()==4) thetaArray[chkPeak][1] = _NthLowest(values2, 2);
//						else {}
//					}

					// assign as scattered if there was only one point at last stage
					if (locs2.size()==1) {
						scattered[chkPeak] = true; // will give a negative vote (exclude it) for this one when it comes to detection
					}

					// USE locs2 check if it is scattered - find inter distance
					if (locs2.size()==2) {
						float interD = (float) (Math.pow(locs2.get(0)[0]-locs2.get(1)[0], 2) + Math.pow(locs2.get(0)[1] - locs2.get(1)[1], 2)); // 0 1
						if (interD>SCATTER_D2) {
							scattered[chkPeak] = true; // might stop fuzzification later
						}
					}
					else if (locs2.size()==3){
						float[] interD = new float[3];
						interD[0] = (float) (Math.pow(locs2.get(0)[0]-locs2.get(1)[0], 2) + Math.pow(locs2.get(0)[1] - locs2.get(1)[1], 2)); // 0 1
						interD[1] = (float) (Math.pow(locs2.get(1)[0]-locs2.get(2)[0], 2) + Math.pow(locs2.get(1)[1] - locs2.get(2)[1], 2)); // 1 2
						interD[2] = (float) (Math.pow(locs2.get(0)[0]-locs2.get(2)[0], 2) + Math.pow(locs2.get(0)[1] - locs2.get(2)[1], 2)); // 0 2
						if (_NthLowest(interD, 1)>SCATTER_D2) {
							scattered[chkPeak] = true;  // might stop fuzzification later
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
							scattered[chkPeak] = true;  // might stop fuzzification later
						}
					}
//					else {
//						// nothing
//					}



				}

			} // thetaArray & locationTopology are done

			//logWriter.println("how scattered? "+Arrays.toString(scattered));
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

//								logWriter.println("substituting cluster idx: "+scatteredIdx+" it was scattered");

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

//								logWriter.println("substituting cluster idx: "+lowestIdx+" it was smallest");

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
//					logWriter.println(		0+		". fuzzy pair -> "+(centralValue - iBackgC)+" : "); // :"+centralValue+" , " + iBackgC + "
//					for (int ee=0; ee<gg.length; ee++) {
//						logWriter.println(	(ee+1)+	". fuzzy pair -> "+(gg[ee][0] - gg[ee][2])+" : "+(gg[ee][1] - gg[ee][2])); // :"+gg[ee][0]+" , "+gg[ee][1]+" , " + gg[ee][2]+"
//					}
//					logWriter.println("fuzzy sc: "+fuzzyValue);

				}
			}
//			else {
//				// peak detection said it is not bifurcation
//			}

			/*
			// log locationTopology
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
			topologyXY.add(locationTopology); // just to show & debug, thetas are used to make the decision
			theta.add(thetaArray);  // no need to save unless logging

		}
		System.out.println(" done.");
		/***********************************************************/

		//new ImagePlus("score", fuzzyScores).show();

		//  EXTRACT CONNECTED REGIONS WITH HIGH RESPONSE
		/***********************************************************/
		System.out.print("extract regions (connected components grouping)... ");
		ByteProcessor score = new ByteProcessor(imp.getWidth(), imp.getHeight());
		for (int ii=0; ii<imp.getWidth()*imp.getHeight(); ii++) if (fuzzyScores.getf(ii) >= MIN_FUZZY_SCORE) score.set(ii, 255);
		Find_Connected_Regions conn_reg = new Find_Connected_Regions(new ImagePlus("", score), true);  // true means save locations
		conn_reg.run("");
		System.out.println(" done.");
		return conn_reg.getConnectedRegions();

	}

	public ArrayList<ArrayList<int[]>> run_old(ImagePlus inimg, double neuronDiameter, float iDiff1)
	{

		// working version of the JunctionDet_v11 plugin before it was modified (wider neighbourhood for local background calculation & averaging along the line)

		if(inimg==null) return null;
		imp = Tools.convertToFloatImage(inimg);

		//scale 	= 1.5f;
		D		= neuronDiameter;
		iDiff 	= iDiff1;

		CPU_NR 		= Runtime.getRuntime().availableProcessors();
		r 			= (float) (scale * D);
		//resolDeg 	= Profiler.getResolDeg(scale);
		fz 			= new Fuzzy(iDiff);

		int totalJobs;

		/***********************************************************/
		int neighbourhoodR = (int) Math.ceil(scale*D);
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
		//ImagePlus inmask = new ImagePlus("inmask", Masker.mask);
		//inmask.show();
		//imp.show();
		//IJ.selectWindow("inimg");
		//inmask.close();
		int nr = Masker.getMaskLocationsTotal();
		float ratio = 100*(float)nr/(totalJobs);
		System.out.println("total "+nr+" given by Masker ("+ratio+" % kept) "+Profiler.getProfilerLocationsTotal()); // done "+((t2-t1)/1000f)+" sec.
		//GenericDialog cont = new GenericDialog("Continue?");
		//cont.addMessage("Is background eliminated? ("+ratio+" % kept)");
		//cont.enableYesNoCancel();
		//cont.showDialog();
		//if (cont.wasCanceled()||!cont.wasOKed()) {inmask.close(); return; }
		//inmask.close();
		/***********************************************************/



		/***********************************************************/
		Profiler.loadTemplate(imp.getProcessor(), Masker.mask, D, scale, W_STD_RATIO_WRT_TO_D, false); // set the template first
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
		//Profiler.getLocation2IndexMapping().show(); // don't show index mapping
		/***********************************************************/


		/***********************************************************/
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
		/***********************************************************/


		/***********************************************************/
		fuzzyScores = new FloatProcessor(imp.getWidth(), imp.getHeight());

		centersXY 	= new ArrayList<int[]>(nr); // from Profiler
		for (int qq=0; qq<nr; qq++) centersXY.add(qq, new int[]{Profiler.locations[qq][0], Profiler.locations[qq][1]});

		anglesRad = Analyzer.exportPeakIdx();
		for (int ii=0; ii<anglesRad.size(); ii++) {
			for (int jj=0; jj<anglesRad.get(ii).size(); jj++) {
				for (int kk=0; kk<anglesRad.get(ii).get(jj).size(); kk++) {
					float newValue = anglesRad.get(ii).get(jj).get(kk) * resolDeg * Deg2Rad;
					anglesRad.get(ii).get(jj).set(kk, newValue);
				}
			}
		}

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
					peaksXY2.add(kk, pkXY);
				}
				peaksXY1.add(jj, peaksXY2);
			}
			peaksXY.add(ii, peaksXY1);
		}

		// extract locations and fuzzy scores
		topologyXY	= new ArrayList<ArrayList<ArrayList<ArrayList<float[]>>>>(nr);
		theta		= new ArrayList<float[][]>(nr);   							// parameters describing intensity
		int[][] currentPeakPoint  	= new int[4][2];
		float[]	currentPeakVals		= new float[4];
		float[] nextPeakPoint		= new float[2];

		for (int ii=0; ii<nr; ii++) { // loop through all locations selected

			int[] currentCenter = centersXY.get(ii);
			float iBackgC = Masker.backgr.getf(currentCenter[0], currentCenter[1]);

			// check how many clusters there could be
			ArrayList<ArrayList<ArrayList<float[]>>> locationTopology = new ArrayList<ArrayList<ArrayList<float[]>>>(4);
			// 4 clusters max. : clusters(1-4) X thread points (1-4) X points in line (1-2)

			float[][] thetaArray        = new float[4][2]; // with parameters
			boolean[] scattered = new boolean[peaksXY.get(ii).get(0).size()]; // per location, in amount of peaks detected
			for (int chkPeak=0; chkPeak<peaksXY.get(ii).get(0).size(); chkPeak++) {

				// 4 points around the peak
				currentPeakPoint[0][0] 	= (int) Math.floor(peaksXY.get(ii).get(0).get(chkPeak)[0]);
				currentPeakPoint[0][1] 	= (int) Math.floor(peaksXY.get(ii).get(0).get(chkPeak)[1]);
				currentPeakVals[0] 		=       imp.getProcessor().getf(currentPeakPoint[0][0], currentPeakPoint[0][1]);

				currentPeakPoint[1][0] = (int) Math.floor(peaksXY.get(ii).get(0).get(chkPeak)[0]);
				currentPeakPoint[1][1] = (int) Math.ceil (peaksXY.get(ii).get(0).get(chkPeak)[1]);
				currentPeakVals[1] 		=       imp.getProcessor().getf(currentPeakPoint[1][0], currentPeakPoint[1][1]);

				currentPeakPoint[2][0] = (int) Math.ceil (peaksXY.get(ii).get(0).get(chkPeak)[0]);
				currentPeakPoint[2][1] = (int) Math.floor(peaksXY.get(ii).get(0).get(chkPeak)[1]);
				currentPeakVals[2] 		=       imp.getProcessor().getf(currentPeakPoint[2][0], currentPeakPoint[2][1]);

				currentPeakPoint[3][0] = (int) Math.ceil (peaksXY.get(ii).get(0).get(chkPeak)[0]);
				currentPeakPoint[3][1] = (int) Math.ceil (peaksXY.get(ii).get(0).get(chkPeak)[1]);
				currentPeakVals[3] 		=       imp.getProcessor().getf(currentPeakPoint[3][0], currentPeakPoint[3][1]);

				// thetaArray[chkPeak][0] is max here
				/****/
				if (useMax)
					thetaArray[chkPeak][0] = maxArray(currentPeakVals);
				else
					thetaArray[chkPeak][0] = _NthLowest(currentPeakVals, 3); // N lowest out of 4

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
							thetaArray[chkPeak][1] = Math.max(thetaArray[chkPeak][1],
																	 Interpolator.interpolateAt(nextPeakPoint[0], nextPeakPoint[1], (FloatProcessor) imp.getProcessor()));

						}
						locationCluster.add(locationThread);

					}

				}

				locationTopology.add(locationCluster);

				// check the topology for this peak, how scattered they are at second ring, use locationCluster
				// in case locationCluster has >=3 samples => 6 or 3 combinations inside

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

					// USE locs2 check if it is scattered - find inter distance
					if (locs2.size()==2) {
						float interD = (float) (Math.pow(locs2.get(0)[0]-locs2.get(1)[0], 2) + Math.pow(locs2.get(0)[1] - locs2.get(1)[1], 2)); // 0 1
						if (interD>SCATTER_D2) {
							scattered[chkPeak] = true; // will stop fuzzification later
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

								//logWriter.println("substituting cluster idx: "+scatteredIdx+" it was scattered");

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

								//logWriter.println("substituting cluster idx: "+lowestIdx+" it was smallest");

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

				}
			}
			else {
				// no need - peak detection said it is not!
			}

			// updates
			topologyXY.add(locationTopology);
			theta.add(thetaArray);

		}
		/***********************************************************/


		/***********************************************************/
		ByteProcessor score = new ByteProcessor(imp.getWidth(), imp.getHeight());
		for (int ii=0; ii<imp.getWidth()*imp.getHeight(); ii++) if (fuzzyScores.getf(ii) >= MIN_FUZZY_SCORE) score.set(ii, 255);
		Find_Connected_Regions conn_reg = new Find_Connected_Regions(new ImagePlus("", score), true);
		conn_reg.run("");
		//int nr_regions = conn_reg.getNrConnectedRegions();

		//ImagePlus imageLabels = conn_reg.showLabels();
		//imageLabels.show();


		return conn_reg.getConnectedRegions();

		//detectionOverlay = formPointOverlay(conn_reg.getConnectedRegions(), MIN_SIZE);
		//resTab = formResultsTable(conn_reg.getConnectedRegions(), fuzzyScores, MIN_SIZE);
		//resTab.show("BIFURCATIONS");
		//int nr_regions = detectionOverlay.size();

//		if (showAll) {
//			// add points used to create regions of connected components
//			for (int ii=0; ii<imp.getWidth()*imp.getHeight(); ii++) {
//				if (score.get(ii)>=255) detectionOverlay.add(new PointRoi(ii%imp.getWidth()+.5, ii/imp.getWidth()+.5));
//			}
//		}

		//System.out.println("detection overlay formed "+detectionOverlay.size()+" detections.");
		//return detectionOverlay;

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

	private static float medAlongLine(float ox1, float oy1, float ox2, float oy2, ImageProcessor ipSource) // , float[] debug
	{
		float 	dl 	= 0.50f;
		float 	l 	= (float) Math.sqrt(Math.pow(ox2 - ox1, 2) + Math.pow(oy2 - oy1, 2));

		float 	avg = 0;
		int 	cnt = 0;

		int w = ipSource.getWidth();
		int h = ipSource.getHeight();

		if (ox1<0 || oy1<0)
			return 0;

//		debug[0] = ox1;
//		debug[1] = oy1;

		// count them first
		for (float x = ox1, y = oy1, ll = 0; ll<=l; x+=dl*(ox2-ox1)/l, y+=dl*(oy2-oy1)/l, ll+=dl) {
			if (x>=w-1 || y>=h-1){
				break;
			}
			cnt++;
		}

		// alloc array
		float[] valsAlongLine = new float[cnt];

		cnt=0;
		for (float x = ox1, y = oy1, ll = 0; ll<=l; x+=dl*(ox2-ox1)/l, y+=dl*(oy2-oy1)/l, ll+=dl) { // x <= ox2 && y <= oy2

			if (x>=w-1 || y>=h-1){
				break;
			}
			valsAlongLine[cnt] = Interpolator.interpolateAt(x, y, (FloatProcessor) ipSource);
//			avg += Interpolator.interpolateAt(x, y, (FloatProcessor) ipSource);
			cnt ++;

		}

		return (float) Tools.median_Wirth(valsAlongLine);

		//return (cnt>0)? avg/cnt : 0;

	}

	private static float maxArray(float[] inArray)
	{
		float maxValue = inArray[0];
		for (int i=1; i<inArray.length; i++) {
			if (inArray[i]>maxValue)
				maxValue = inArray[i];
		}
		return maxValue;
	}

	private static float _NthLowest(float[] array, int N)
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

	private static float _NthLowest(ArrayList<Float> array, int N)
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

	public Overlay formPointOverlay(Vector<float[]> detectionList, Color col)
	{

		Overlay detections = new Overlay();

		for (int a=0; a<detectionList.size(); a++) {

			float Cx = detectionList.get(a)[0];
			float Cy = detectionList.get(a)[1];
			float R  = detectionList.get(a)[2];
			OvalRoi ovroi = new OvalRoi(Cx-R+.5, Cy-R+.5, 2*R, 2*R);
			ovroi.setStrokeWidth(3);
			ovroi.setStrokeColor(col);
			detections.add(ovroi);

		}

		return detections;

	}

	public Vector<float[]> formDetectionList(ArrayList<ArrayList<int[]>> regs, int minSize) // float[] -> [x, y, r]
	{

		Vector<float[]> out = new Vector<float[]>();

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

				out.add(new float[]{Cx, Cy, R});

			}
		}

		// prune those smallest that are around the big ones - 2d clustering
		//return clustering(out);

		return out;

	}

	public ResultsTable formResultsTable(ArrayList<ArrayList<int[]>> regs, FloatProcessor fuzzyScores, int minSize)
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

    public void formChimeraScript(ArrayList<ArrayList<int[]>> regs, FloatProcessor fuzzyScores, int minSize){
    // will output UCSF Chimera .com file that can be loaded with image for visuelization

        PrintWriter outWriter   = null;
        String      outName     = "chimeraList"+".com"; // imp.getTitle()



        // initialize detection log file
        try {
            outWriter = new PrintWriter(new BufferedWriter(new FileWriter(outName, true)));
            // // dateFormat.format(date)
        } catch (IOException e) {}


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
                // (int)Math.round()
                outWriter.println("shape sphere radius "+Math.sqrt(regs.get(i).size() / Math.PI)+" center "+IJ.d2s(Cx, 3, 5)+","+IJ.d2s(Cy, 3, 5)+",0 color red");

//                rt.incrementCounter();
//                rt.addValue("X", Cx);
//                rt.addValue("Y", Cy);
//                rt.addValue("size", regs.get(i).size());
//                rt.addValue("avg. bness.", avgBifurcationess);
//                rt.addValue("min. bness.", minBifurcationess);
//                rt.addValue("max. bness.", maxBifurcationess);

            }
        }


        System.out.println(new File(outName).getAbsolutePath()+" is exported");
        outWriter.close();

    }

	public static int[] clustering(Vector<float[]> disks, float tolerance) //list with [x, y, r]
	{

		// this method will cluster all the disks that overlap or are within some tolerance in the same cluster (give them the same label)

		int[] labels = new int[disks.size()];
		for (int i = 0; i < labels.length; i++) labels[i] = i;

//		System.out.println("INIT. LABELS:");
//		for (int i = 0; i < labels.length; i++)
//			System.out.print(labels[i]+" ");
//		System.out.println();

		for (int i = 0; i < disks.size(); i++) {

			// one versus the rest
			for (int j = 0; j < disks.size(); j++) {

				if (i!=j) {

					double dst2 	= Math.pow(disks.get(i)[0]-disks.get(j)[0], 2) + Math.pow(disks.get(i)[1]-disks.get(j)[1], 2);
					double rd2 		= Math.pow(disks.get(i)[2]+disks.get(j)[2], 2);
					double tol2		= Math.pow(tolerance, 2);

					if (dst2<=rd2 || dst2<=tol2) {  // they are neighbours

						if (labels[j]!=labels[i]) {

							int currLabel = labels[j];
							int newLabel  = labels[i];

							labels[j] = newLabel;

							//set all that also were currLabel to newLabel
							for (int k = 0; k < labels.length; k++)
								if (labels[k]==currLabel)
									labels[k] = newLabel;

						}

					}

				}

			}

		}

//		System.out.println("OUT LABELS:");
//		for (int ii = 0; ii < labels.length; ii++)
//			System.out.print(labels[ii]+" ");
//		System.out.println();

		return labels; // cluster labels for each disc

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

			// it exists already, store it in the list for each loc
			for (int loopLocations=0; loopLocations<profilesToAppend.length; loopLocations++) {

				// profileToAdd from Profiler.profiles
				float[] profileToAdd = new float[profilesToAppend[loopLocations].length];
				for (int k=0; k<profileToAdd.length; k++) {
					profileToAdd[k] = profilesToAppend[loopLocations][k];
				}

				profileListToUpdate.get(loopLocations).add(profileToAdd);

			}

		}

	}

	public static Vector<float[]> groupClusters(int[] labels, Vector<float[]> detections){

		Vector<float[]> afterClustering = new Vector<float[]>();

		boolean[] checked = new boolean[labels.length];

		// loop labels to recreate list
		for (int l = 0; l < labels.length; l++) {

			if (!checked[l]) {

				checked[l] = true;
				afterClustering.add(detections.get(l).clone());

				int currLabel = labels[l];

				for (int j = l+1; j < labels.length; j++) {

					if (labels[j]==currLabel) {

						checked[j] = true;

						if (detections.get(l)[2]>detections.get(j)[2]) {
							int currentListIdx = afterClustering.size()-1;
							afterClustering.set(currentListIdx, detections.get(l).clone());
						}

					}

				}

			}

		}

		return afterClustering;

	}

//	public static float[] classificationPerformance(Vector<float[]> dets, float[][] gndTth, float tolerance) // [x,y,r],[x,y], tolerance to location acc.
//	{
//
//		float tp=0, fp=0, fn=0; // tp + fp == dets.size()   fn -> those that were missed
//		boolean[] isHit = new boolean[gndTth.length];
//
//		for (int detIdx = 0; detIdx < dets.size(); detIdx++) {
//
//			boolean isFoundInGndTth = false;
//
//			// compare with each gndTth element whether it is covered & check if yes
//			for (int gndTthIdx = 0; gndTthIdx < gndTth.length; gndTthIdx++) {
//
//
//				float squaredDistFromGndTthElement = dist2dSquared(gndTth[gndTthIdx][0], gndTth[gndTthIdx][1],
//																		  dets.get(detIdx)[0], dets.get(detIdx)[1]);
//
//				if (squaredDistFromGndTthElement<=Math.pow(tolerance,2)) { // dets.get(detIdx)[2]
//
//					isFoundInGndTth 	= true;
//					isHit[gndTthIdx] 	= true;
//					tp++;
//					break; // stop checking the rest
//
//				}
//
//			}
//
//			if (!isFoundInGndTth) fp++;
//
//		} // tp & fp are calculated
//
//		// check for those that were missed
//		for (int ii=0; ii<isHit.length; ii++)
//			if (!isHit[ii]) fn++;
//
//		float[] out = new float[]{tp, fp, fn};
//		return out;
//
//	}

	public static Overlay classificationPerformanceOverlay(Vector<float[]> dets, float[][] gndTth, float tolerance, float[] eval) // [x,y,r],[x,y], tolerance to location acc.
	{

		Overlay ov = new Overlay();

		// plot annotations with tolerance
		for (int gndTthIdx = 0; gndTthIdx < gndTth.length; gndTthIdx++) {
			OvalRoi o = new OvalRoi(
										   gndTth[gndTthIdx][0]-tolerance+.5,
										   gndTth[gndTthIdx][1]-tolerance+.5,
										   2*tolerance,
										   2*tolerance);
			o.setStrokeColor(Color.BLUE);
			ov.add(o);
		}

		float tp=0, fp=0, fn=0; // tp + fp == dets.size()   fn -> those that were missed
		boolean[] isHit = new boolean[gndTth.length];

		for (int detIdx = 0; detIdx < dets.size(); detIdx++) {

			boolean isFoundInGndTth = false;

			// compare with each gndTth element whether it is covered & check if yes
			for (int gndTthIdx = 0; gndTthIdx < gndTth.length; gndTthIdx++) {


				float squaredDistFromGndTthElement = dist2dSquared(gndTth[gndTthIdx][0], gndTth[gndTthIdx][1],
																		  dets.get(detIdx)[0], dets.get(detIdx)[1]);

				if (squaredDistFromGndTthElement<=Math.pow(tolerance, 2)) { // dets.get(detIdx)[2]

					isFoundInGndTth 	= true;
					isHit[gndTthIdx] 	= true;
					tp++;

					// ADD TP
					PointRoi p = new PointRoi(dets.get(detIdx)[0]+.5, dets.get(detIdx)[1]+.5);
					p.setStrokeColor(Color.RED);
					ov.add(p);

					break; // stop checking the rest

				}

			}

			if (!isFoundInGndTth) {
				fp++;

				// ADD FP
				PointRoi p = new PointRoi(dets.get(detIdx)[0]+.5, dets.get(detIdx)[1]+.5);
				p.setStrokeColor(Color.YELLOW);
//				OvalRoi o = new OvalRoi(
//											   dets.get(detIdx)[0]-dets.get(detIdx)[2]+.5,
//											   dets.get(detIdx)[1]-dets.get(detIdx)[2]+.5,
//											   2*dets.get(detIdx)[2],
//											   2*dets.get(detIdx)[2]);
//				o.setStrokeColor(Color.YELLOW);
//				ov.add(p);
				ov.add(p);
			}

		} // tp & fp are calculated

		// check for those that were missed
		for (int ii=0; ii<isHit.length; ii++)
			if (!isHit[ii]) {
				fn++;

				// ADD FN
				PointRoi p = new PointRoi(gndTth[ii][0]+.5, gndTth[ii][1]+.5);
				p.setStrokeColor(Color.BLUE);
				ov.add(p);
			}

		//float[] out = new float[]{tp, fp, fn};
		eval[0] = tp;
		eval[1] = fp;
		eval[2] = fn;

		return ov;

	}

	private static float dist2dSquared(float x1, float y1, float x2, float y2)
	{
		return (float) (Math.pow(x2-x1, 2)+Math.pow(y2-y1, 2));
	}

}