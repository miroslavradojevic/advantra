package detection;

import aux.Tools;
import conn.Find_Connected_Regions;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.*;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 7/7/13
 * Time: 8:13 PM
 */
public class JunctionDet implements PlugInFilter, MouseListener {

	ImagePlus 	    inimg;
    ByteProcessor   scoreimg;

    // visualization  (for selected locations)
    ImagePlus       vizProfileImage;
    ImageStack      vizProfileStack;

    ArrayList<ArrayList<float[]>>   profiles;
    ArrayList<ArrayList<Integer>>   angularRes;
    ArrayList<ArrayList<String>>    profilesName;

    int                             totalCfgs;
//    double                          neuronD;
    double[]                        scales;     // per cfg

	// ms output
    ArrayList<ArrayList<float[]>>   angles;   	// 3 angles  (necessary for matching)
    ArrayList<ArrayList<float[]>>   idxs;   	// 3 detection indexes
    ArrayList<ArrayList<float[]>>   peaks;      // each angle one peak score

    ArrayList<Float>   				bacgr;      // per loc
    ArrayList<Float>   	            centrAvg;   // per loc

	Overlay detections;

    // values used to form the configuration
    // mean-shift will be supplying with these
    double[][][]    P; // point coordinates
    double[][]      G; // peak values of the detection
    double[]        A; // angles (deg)
    double[]        R; // will hold radiuses for each configuration

    /*
    parameters
     */
    double neuronDiamMax, D;//, eps, dMin; // scaleMin, scaleMax, neuronDiamMin
    float iDiff;
	int CPU_NR;//H, M, MaxIter;
	int weightStdRatioToD;

	public int setup(String s, ImagePlus imagePlus) {

		if(imagePlus==null) return DONE;
		inimg = Tools.convertToFloatImage(imagePlus);
		inimg.setTitle("inimg");
		return DOES_8G+DOES_32+NO_CHANGES;

	}

	public void run(ImageProcessor imageProcessor) {

        int totalJobs; // used in paralelization lines
        long t1, t2;

		// fix scales
		scales = new double[3];
		scales[0] = 1.5;
		scales[1] = 2.0;
		scales[2] = 2.5;

		// profile shape some middle bell
		weightStdRatioToD = 4;

        totalCfgs = scales.length;

		// these will be modified each time extractScore is called
        P = new double[3][scales.length][2]; 	// point coordinates per cluster per point 2d
        G = new double[3][scales.length]; 		// values of the peaks per cluster per point
        A = new double[3+1]; 					// sums votes per cluster elements plus the central

		R = new double[scales.length];
        /*
        generic dialog
        */
        //neuronDiamMin       =  Prefs.get("advantra.critpoint.neuronDiamMin", 3);
        neuronDiamMax       =  Prefs.get("advantra.critpoint.neuronDiamMax", 3);
        //scaleMin            =  Prefs.get("advantra.critpoint.scaleMin", 1.5);
        //scaleMax            =  Prefs.get("advantra.critpoint.scaleMax", 2.5);
        D                   =  Prefs.get("advantra.critpoint.D", 10);
//        H       = (int) Prefs.get("advantra.critpoint.ms2d.H", 4);
//        MaxIter = (int) Prefs.get("advantra.critpoint.ms2d.MaxIter", 200);
//        eps     =       Prefs.get("advantra.critpoint.ms2d.eps", 0.0001);
//        dMin       =       Prefs.get("advantra.critpoint.ms2d.d", 0.5);
//        M       = (int) Prefs.get("advantra.critpoint.ms2d.M", 5);

		iDiff = 5;

        CPU_NR              = (int) Prefs.get("advantra.critpoint.CPU_NR", 4);

        GenericDialog gd = new GenericDialog("JUNCTION DET.");
        gd.addNumericField("neuronD ",  neuronDiamMax, 0, 10, " MAX");
        gd.addNumericField("D",         D, 1, 10, "score param.");

		/*gd.addMessage("--- LOCAL PEAKS (MS on 2D image) ---");
        gd.addNumericField("H ",        H,      0, 10, "spatial neighbourhood");
        gd.addNumericField("Max Iter",  MaxIter,0, 10, "max #iterations");
        gd.addNumericField("Eps",       eps,    6, 10, "convergence threshold");
        gd.addMessage("---");
        gd.addNumericField("d",         dMin,      1, 10, "min intra-cluster distance");
        gd.addNumericField("M",         M,      0, 10, "min # points in cluster");
        */

        gd.addNumericField("CPU_NR ",   CPU_NR, 0, 10, "");

        gd.showDialog();
        if (gd.wasCanceled()) return;

        neuronDiamMax       =  gd.getNextNumber();
        Prefs.set("advantra.critpoint.neuronDiamMax",   neuronDiamMax);
        D                   =  gd.getNextNumber();
        Prefs.set("advantra.critpoint.D",               D);
		CPU_NR              =  (int) gd.getNextNumber();
        Prefs.set("advantra.critpoint.CPU_NR", 	        CPU_NR);

        /*
        ***********************************************************
         */
        IJ.log("extracting background...");
        t1 = System.currentTimeMillis();
        int neighbourhoodR = (int) Math.ceil(4*neuronDiamMax);
        Masker.loadTemplate(inimg.getProcessor(), neighbourhoodR, iDiff);
        totalJobs = inimg.getHeight()*inimg.getWidth();
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
        IJ.log("done "+((t2-t1)/1000f)+" sec.");

        inimg.show();

		// overlay mask
		ImagePlus maskViz = new ImagePlus("inmask", Masker.mask);
        maskViz.show();
        IJ.selectWindow("inimg");
        IJ.run("Add Image...", "image=inmask x="+0+" y="+0+" opacity=30");
        maskViz.close();

        /*
        score image initialize
         */
        scoreimg = new ByteProcessor(inimg.getWidth(), inimg.getHeight());
		//Profiler.loadTemplate(inimg.getProcessor(), Masker.mask);
        // can see the locations now

        int totalLocations = Masker.getMaskLocationsTotal();

        IJ.log("total "+totalLocations+" locations given by mask ("+(100*(float)totalLocations/(inimg.getWidth()*inimg.getHeight()))+" percent kept)");

        profiles     	= new ArrayList<ArrayList<float[]>>(totalLocations);
		angularRes    	= new ArrayList<ArrayList<Integer>>(totalLocations);
        profilesName  	= new ArrayList<ArrayList<String>>(totalLocations);

        /*
        ***********************************************************
         */

		IJ.log("calculating profiles... ");
		t1 = System.currentTimeMillis();

		//int totalConfs = 0;
        //for (double neuronDiam = neuronDiamMax; neuronDiam<=neuronDiamMax; neuronDiam++) {

            for (int scaleIdx=0; scaleIdx<totalCfgs; scaleIdx++) {

                R[scaleIdx] = neuronDiamMax*scales[scaleIdx];
                Profiler.loadTemplate(inimg.getProcessor(), Masker.mask, neuronDiamMax, scales[scaleIdx], weightStdRatioToD, false);

                totalJobs = Profiler.offsets.size();
                Profiler profiler_jobs[] = new Profiler[CPU_NR];

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

                updateList(scales[scaleIdx]);

            }

		t2 = System.currentTimeMillis();
		IJ.log("done "+((t2-t1)/1000f)+" sec.");

        /*
        ***********************************************************
         */

        IJ.log("calculating peaks... ");
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

        /*
        ***********************************************************
         */

        IJ.log("clustering corresponding peaks...");

        angles    	= new ArrayList<ArrayList<float[]>>(totalLocations);
        idxs    	= new ArrayList<ArrayList<float[]>>(totalLocations);
        peaks       = new ArrayList<ArrayList<float[]>>(totalLocations);
		centrAvg   	= new ArrayList<Float>(totalLocations);
		bacgr       = new ArrayList<Float>(totalLocations);

        for (int locIdx=0; locIdx<totalLocations; locIdx++) {

            int atX = Profiler.locations[locIdx][0];
            int atY = Profiler.locations[locIdx][1];

            ArrayList<float[]>  anglesPerCfg = new ArrayList<float[]>(totalCfgs);
            ArrayList<float[]>  idxsPerCfg = new ArrayList<float[]>(totalCfgs);
            ArrayList<float[]>  peaksPerCfg = new ArrayList<float[]>(totalCfgs);

            boolean isFirst = true;
            float[] refAnglesDeg    = new float[3];

            for (int cfgIdx=0; cfgIdx<totalCfgs; cfgIdx++) {

                /*
                 find central bloc score for this location
                  */

//                centralAvgPerCfg.add();

                // get indexes form Analyzer for this conf & loc
                if (Analyzer.peakIdx.get(locIdx).get(locIdx).size() >=3) {   // there was a change in Analyzer

                    /*
                    converging angles and scores
                     */

                    float[] anglesDeg       = new float[3];
                    float[] indexes       	= new float[3];
                    float[] peaks          	= new float[3];

                    indexes[0] = Analyzer.peakIdx.get(locIdx).get(cfgIdx).get(0);
                    indexes[1] = Analyzer.peakIdx.get(locIdx).get(cfgIdx).get(1);
                    indexes[2] = Analyzer.peakIdx.get(locIdx).get(cfgIdx).get(2);

                    anglesDeg[0] = Analyzer.peakIdx.get(locIdx).get(cfgIdx).get(0) * angularRes.get(locIdx).get(cfgIdx);
                    anglesDeg[1] = Analyzer.peakIdx.get(locIdx).get(cfgIdx).get(1) * angularRes.get(locIdx).get(cfgIdx);
                    anglesDeg[2] = Analyzer.peakIdx.get(locIdx).get(cfgIdx).get(2) * angularRes.get(locIdx).get(cfgIdx);

                    peaks[0] = (float) Tools.interp1Darray(indexes[0], profiles.get(locIdx).get(cfgIdx));
                    peaks[1] = (float) Tools.interp1Darray(indexes[1], profiles.get(locIdx).get(cfgIdx));
                    peaks[2] = (float) Tools.interp1Darray(indexes[2], profiles.get(locIdx).get(cfgIdx));

                    if (isFirst) {

                        refAnglesDeg[0] = anglesDeg[0];// = indexes[0] * angularRes.get(locIdx).get(cfgIdx);
                        refAnglesDeg[1] = anglesDeg[1];// = indexes[1] * angularRes.get(locIdx).get(cfgIdx);
                        refAnglesDeg[2] = anglesDeg[2];// = indexes[2] * angularRes.get(locIdx).get(cfgIdx);

                        isFirst = false;

                    }
                    else {
                        // not the first time, check angular matches
                        // read them and store wrt first

//                        IJ.log("** not first,  **"+atX+" , "+atY+" , cfg="+ cfgIdx);
//                        IJ.log(Arrays.toString(anglesDeg));
//                        IJ.log(Arrays.toString(refAnglesDeg));
//                        IJ.log("---");

                        // rearrange wrt. to the first
                        // match using Hungarian algorithm
                        boolean[][] chkd = new boolean[3][3];
                        double[][] 	dst2 = new double[3][3];

                        for (int i=0; i<3; i++) {
                            for (int j=0; j<3; j++) {
                                dst2[i][j] = Math.abs(Tools.wrap_180(anglesDeg[i] - refAnglesDeg[j]));
                            }
                        }

                        int[] mapping = new int[3];

                        for (int check=0; check<3; check++) {

                            double dst2Min = Double.MAX_VALUE;
                            int imin = -1;
                            int jmin = -1;

                            for (int i=0; i<3; i++) {

                                for (int j=0; j<3; j++) {
                                    if (!chkd[i][j] && dst2[i][j]<dst2Min) {
                                        dst2Min = dst2[i][j];
                                        imin = i;
                                        jmin = j;
                                    }
                                }

                            }

                            // row imin in chkd to true
                            for (int w=0; w<3; w++) chkd[imin][w] = true;
                            // col jmin in chkd to true
                            for (int w=0; w<3; w++) chkd[w][jmin] = true;

                            mapping[imin]=jmin;
//                            IJ.log("map "+imin+" goes to "+mapping[imin]+" ang. dst was "+dst2[imin][jmin]);

                        }

                        float[] anglesDegTemp   = new float[3]; // ugly solution to swap 3
                        float[] indexesTemp     = new float[3]; // ugly solution to swap 3
                        float[] peaksTemp      = new float[3];

                        for (int t=0; t<mapping.length; t++) {
                            anglesDegTemp[mapping[t]] = anglesDeg[t];
                            indexesTemp[mapping[t]] = indexes[t];
                            peaksTemp[mapping[t]]   = peaks[t];
                        }

//                        IJ.log("...");
//                        IJ.log(Arrays.toString(anglesDegTemp));
//                        IJ.log("...");

                        for (int t=0; t<mapping.length; t++) {
                            indexes[t]      = indexesTemp[t];
                            anglesDeg[t]    = anglesDegTemp[t];
                            peaks[t]        = peaksTemp[t];
                        }

                    }

//                    IJ.log("***");
//                    IJ.log(Arrays.toString(anglesDeg));
//                    IJ.log(Arrays.toString(refAnglesDeg));
//                    IJ.log("---");

                    idxsPerCfg.add(indexes);
                    anglesPerCfg.add(anglesDeg);
                    peaksPerCfg.add(peaks);

                }
                else {

                    anglesPerCfg.add(null);
                    peaksPerCfg.add(null);
                    idxsPerCfg.add(null);

                }

            }

            angles.add(anglesPerCfg);
            peaks.add(peaksPerCfg);
            idxs.add(idxsPerCfg);

            centrAvg.add(centralAvg(atX, atY, neuronDiamMax, (FloatProcessor) inimg.getProcessor()));
            bacgr.add(Masker.backgr.getf(atX, atY)); // take the value from Masker

        }

        IJ.log("done...");

        // loop once again and extract scores

        for (int locIdx=0; locIdx<totalLocations; locIdx++) {

			double score =  extractScore(locIdx, false);

            scoreimg.set(Profiler.locations[locIdx][0], Profiler.locations[locIdx][1], ((int) score));

        }

        ImagePlus scoreImagePlus = new ImagePlus("score", scoreimg);
		scoreImagePlus.show();


        /*
        -----------------------------------------------
         */

        t1 = System.currentTimeMillis();
		Find_Connected_Regions conn_reg = new Find_Connected_Regions(scoreImagePlus, true);
		conn_reg.run("");
		t2 = System.currentTimeMillis();
		int nr_regions = conn_reg.getNrConnectedRegions();

		IJ.log(nr_regions+" connected regions extracted.\n" +
					   "elapsed: "+((t2-t1)/1000f)+ " seconds.");

        ArrayList<ArrayList<int[]>> regs = conn_reg.getConnectedRegions();

		ImagePlus imageLabels = conn_reg.showLabels();
        imageLabels.show();

		detections = new Overlay();
		double VERY_SMALL_POSITIVE = 0.000001;

        for (int i=0; i<regs.size(); i++) {
            if (regs.get(i).size()>1) {
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

//                IJ.log("A: "+A+" , B: "+B+"");
//                if(regs.get(i).get(0)[1]==584 && regs.get(i).get(0)[0]==254) {
//                    IJ.log("debug");
//                    IJ.log("x="+regs.get(i).get(0)[1]+", y="+regs.get(i).get(0)[0]+", A="+A+", B="+B+", angle="+ellipseParams[4]);
//                }

                float scalePlot = 2;
                detections.add(Tools.drawEllipse(ellipseParams[1]+0.5f, ellipseParams[0]+0.5f, scalePlot*A, scalePlot*B, ellipseParams[4], Color.YELLOW, 2, 50));

            }
        }

        // detected spots  from the score
//		for (int i1=0; i1<scoreimg.getWidth(); i1++) {
//			for (int i2=0; i2<scoreimg.getHeight(); i2++) {
//				if (scoreimg.get(i1, i2)==255) {
//					PointRoi pt = new PointRoi(i1+0.5, i2+0.5);
//					detections.add(pt);
//				}
//
//			}
//		}

        inimg.setOverlay(detections);
        //scoreImagePlus.setOverlay(detections);

		inimg.getCanvas().addMouseListener(this); // prepare for mouse interaction

        vizProfileImage = new ImagePlus();

//        // overlay mask
//        ImagePlus inmaskImg = new ImagePlus("inmask", Masker.maskip);
//        inmaskImg.show();
//        IJ.selectWindow("inimg");
//        IJ.run("Add Image...", "image=inmask x="+0+" y="+0+" opacity=20");
//        inmaskImg.close();

//        // overlay score
//        ImagePlus scoreImg = new ImagePlus("score", scoreimg);
//        scoreImg.show();
//        scoreImg.setOverlay(ov);
//        IJ.selectWindow("inimg");
//        IJ.run("Add Image...", "image=score x="+0+" y="+0+" opacity=70");
//        scoreImg.close();

        IJ.selectWindow("inimg");
        IJ.setTool("hand");
        inimg.getCanvas().zoomIn(0, 0);

	}

	private double extractScore(int locIdx, boolean printVals) {

        // reset arrys to zero each time
        Arrays.fill(A, 0);
        for (int q=0; q<G.length; q++) {
            Arrays.fill(G[q], 0);
        }
        for (int q=0; q<P.length; q++) {
            for (int q1=0; q1<P[0].length; q1++) {
                Arrays.fill(P[q][q1], 0);
            }
        }

        A[0] = centrAvg.get(locIdx) - bacgr.get(locIdx);
        A[0] = (A[0]>D)? 1 : 0; // *bacgr.get(locIdx)
		if (printVals) IJ.log("A0 ("+centrAvg.get(locIdx)+") = "+A[0]); //A0 = 1-Math.exp(-A0/D);

		int cntCfgs = 0;

		for (int cfgIdx=0; cfgIdx<totalCfgs; cfgIdx++) {

			if (angles.get(locIdx).get(cfgIdx)!=null) { // if ms gave 3 values

				// sorted peak indexes
				float[] sortedIdxs = new float[3];
				for (int i1=0; i1<sortedIdxs.length; i1++) sortedIdxs[i1] = idxs.get(locIdx).get(cfgIdx)[i1];
				Arrays.sort(sortedIdxs);

				//
				if (cfgIdx>0 && !isProfileValid(sortedIdxs, profiles.get(locIdx).get(cfgIdx), bacgr.get(locIdx)+1)) { // +1 in case background was zero
					if (printVals) IJ.log("detection not valid");
					break;
				}

                    for (int clusterIdx = 0; clusterIdx<3; clusterIdx++) {

                        P[clusterIdx][cfgIdx][0] =  R[cfgIdx] * Math.cos(angles.get(locIdx).get(cfgIdx)[clusterIdx]*(Math.PI/180));
                        P[clusterIdx][cfgIdx][1] = -R[cfgIdx] * Math.sin(angles.get(locIdx).get(cfgIdx)[clusterIdx]*(Math.PI/180));

                        double temp = peaks.get(locIdx).get(cfgIdx)[clusterIdx]- bacgr.get(locIdx);
                        temp = (temp>D)? 1 : 0 ; // 1-Math.exp(-A1t/D) // *bacgr.get(locIdx)
                        if (printVals) IJ.log("\ncluster "+clusterIdx+" pk: "+peaks.get(locIdx).get(cfgIdx)[clusterIdx]+", bgr: "+bacgr.get(locIdx)+" , "+temp);
                        A[clusterIdx+1] +=  temp; // because first index was for central value

                    }

				cntCfgs++;

			}
			else {
				if (printVals) IJ.log("not enough points");
				break;
			}

		}

		double score = 0;

        // check  directionality
		if (cntCfgs==3) {

            boolean alined = true;

            for (int clusterIdx=0; clusterIdx<3; clusterIdx++) {

				// global aline
                for (int cfgIdx=1; cfgIdx<totalCfgs; cfgIdx++) {//first one cannot refer to anything before itself

                    if (!isDirectional(P[clusterIdx][cfgIdx], P[clusterIdx][cfgIdx-1], R[cfgIdx], R[cfgIdx-1])) {

                        alined = false;
                        if (printVals) IJ.log("not alined globally in cluster "+clusterIdx);
                        break; // and get out with score 0

                    }

                }
				// local aline
				for (int cfgIdx=2; cfgIdx<totalCfgs; cfgIdx++) {
					if ( !isLocalDirectional(P[clusterIdx][cfgIdx], P[clusterIdx][cfgIdx-1], P[clusterIdx][cfgIdx-2], 45) ) {
						alined = false;
						if (printVals) IJ.log("not alined locally in cluster "+clusterIdx);
						break;
					}
				}

            }

            if (alined) {

                score = (A[0]==1 && (A[1] + A[2] + A[3])>=8)? 255: 0;
                if (printVals) IJ.log("VOTING SCORE = "+score);

            }

		}
		else {
			score = 0; // there was no 3 in a row
		}

		return score;

//            if (A0>0 && A1>0 && A2>0 && A3>0) {
////                score += Math.log(A0);
////                score += Math.log(A1);
////                score += Math.log(A2);
////                score += Math.log(A3);
////                score = Math.exp(score/4);
//            }


	}

    private boolean isDirectional(double[] pNew, double[] pPrev, double rNew, double rPrev)
    {

        return ((pNew[0]*pPrev[0] + pNew[1]*pPrev[1]) / (rNew*rPrev))>rPrev/rNew;

    }

	private boolean isLocalDirectional(double[] pNew, double[] pPrev, double[] pPrevPrev, double angDiffDeg)
	{

			double v1x, v1y, iv1, v2x, v2y, iv2;
			v1x = pNew[0]-pPrev[0];
			v1y = pNew[1]-pPrev[1];
			v2x = pPrev[0]-pPrevPrev[0];
			v2y = pPrev[1]-pPrevPrev[1];
			iv2 = Math.sqrt(v2x*v2x+v2y*v2y);
			iv1 = Math.sqrt(v1x*v1x+v1y*v1y);
			return Math.acos((v1x*v2x+v1y*v2y) / (iv1 * iv2))*(180/Math.PI)<angDiffDeg;

	}

    private boolean isDirectional(double[] pNew, double[] pPrev)  // slower version
    {
        double pNewNorm     = Math.sqrt(pNew[0]*pNew[0]+pNew[1]*pNew[1]);
        double pPrevNorm    = Math.sqrt(pPrev[0]*pPrev[0]+pPrev[1]*pPrev[1]);
        return (pNew[0]*pPrev[0] + pNew[1]*pPrev[1]) / (pNewNorm*pPrevNorm) > pPrevNorm/pNewNorm;
    }

    private void updateList(double currScale1) // uses Profiler to update  profiles, profilesName, angularRes
    {
        // store it in the list for each loc
        for (int loopLocations=0; loopLocations<Profiler.locations.length; loopLocations++) {

            // profileToAdd from Profiler.profiles
            float[] profileToAdd = new float[Profiler.profiles[loopLocations].length];
            for (int k=0; k<profileToAdd.length; k++) {
                profileToAdd[k] = Profiler.profiles[loopLocations][k];
            }

            // stringToAdd from Profiler.neuronDiam, Profiler.scale, Profiler.resolDeg
            String stringToAdd = "nD_"+Profiler.neuronDiam+"_s_"+Profiler.scale+"_a_"+Profiler.getResolDeg(currScale1);

            // intToAdd from Profiler.resolDeg
            int intToAdd = Profiler.getResolDeg(currScale1);

            if (profiles.size()<Profiler.locations.length) {

                ArrayList<float[]> A = new ArrayList<float[]>();
                A.add(profileToAdd);
                profiles.add(A);

                ArrayList<String> B = new ArrayList<String>();
                B.add(stringToAdd);
                profilesName.add(B);

                ArrayList<Integer> C = new ArrayList<Integer>();
                C.add(intToAdd);
				angularRes.add(C);

            }
            else {

				profiles.get(loopLocations).add(profileToAdd);
				profilesName.get(loopLocations).add(stringToAdd);
				angularRes.get(loopLocations).add(intToAdd);

			}

        }
    }

	private boolean isProfileValid(float[] peakIndexes, float[] profileToCheck, float threshold)
	{

		//Arrays.sort(peakIndexes);

		boolean profileValid = true;

		for (int i=0; i<peakIndexes.length; i++) {

			boolean intervalOK = false;

			int startInterval = Math.round(peakIndexes[i]);
			startInterval = (startInterval<0)? 0 : startInterval;
			startInterval = (startInterval>=profileToCheck.length)? profileToCheck.length-1 : startInterval;

			int endInterval = Math.round( (i==peakIndexes.length-1)? peakIndexes[0] : peakIndexes[i+1] );
			endInterval = (endInterval<0)? 0 : endInterval;
			endInterval = (endInterval>=profileToCheck.length)? profileToCheck.length-1 : endInterval ;

//			IJ.log("check between "+startInterval+" and "+endInterval);

			if (endInterval>=startInterval) {
				for (int k=(int)Math.floor(startInterval); k<=(int)Math.ceil(endInterval); k++) {
					if (profileToCheck[k]<threshold) {
						intervalOK = true;
						break;
					}
				}
			}
			else {
				endInterval += profileToCheck.length;
				for (int k=(int)Math.floor(startInterval); k<=(int)Math.ceil(endInterval); k++) {
					int idxToTake = (k>=profileToCheck.length)? (k-profileToCheck.length) : k ;
					if (profileToCheck[idxToTake]<threshold) {
						intervalOK = true;
						break;
					}
				}
			}

			if (intervalOK==false) {
				profileValid = false;
				break;
			}

		}

		return profileValid;

	}

    public float centralAvg(int atX, int atY, double diam, FloatProcessor inip)
    {
        float out = 0;
        int cnt = 0;
        int margin = (int) Math.ceil(diam/2);

        if(atX>=margin && atY>=margin && atX<inip.getWidth()-margin && atY<inip.getHeight()-margin) {

            for (double x = atX-margin; x <= atX+margin; x+=0.5) {
                for (double y = atY-margin; y <= atY+margin; y+=0.5) {

                    if ((x-atX)*(x-atX)+(y-atY)*(y-atY)<= margin*margin) {
                        out += Interpolator.interpolateAt(x, y, inip);
                        cnt++;
                    }

                }
            }

        }

        return (cnt>0)? (out/cnt) : 0;    // /minSumsPerOrt

    }

    public void mouseClicked(MouseEvent e) {

        ImageCanvas srcCanv = (ImageCanvas) e.getSource();
        int atX = 	srcCanv.offScreenX(e.getX());
        int atY = 	srcCanv.offScreenY(e.getY());

        vizProfileStack = new ImageStack(600, 300);

        String fileName = "plotProfiles.r";
        // empty the file
        PrintWriter writer = null;
        try {
            writer = new PrintWriter(fileName);
        } catch (FileNotFoundException ex) {
            //ex.printStackTrace();
        }
        writer.print("");
        writer.close();

        for (int i=0; i<Profiler.locations.length; i++) {  // i will loop locations
            if (Profiler.locations[i][0]==atX && Profiler.locations[i][1]==atY) {

                extractScore(i, true); // with IJ.log

                /*
				overlay for the main image
			    */
                Overlay ov = detections.duplicate();
                PointRoi pt = new PointRoi(atX+.5, atY+.5);
                ov.add(pt);

                for (int g=0; g<totalCfgs; g++) {  // g will loop configurations

                    int currProfileLen = profiles.get(i).get(g).length;

                    /*
                    plot ij  template
                     */

                    int plotPts = 500;

                    float[] plotX = new float[plotPts];
                    for (int q=0; q<plotPts; q++) plotX[q] = q*((float)currProfileLen/plotPts);// * angularRes.get(i).get(g);

                    // interpolate detection to be plotted
                    float[] plotY = new float[plotPts];
                    for (int q=0; q<plotPts; q++) plotY[q] = (float) (Tools.interp1Darray(plotX[q], profiles.get(i).get(g)));

                    // turn plotX to real angles for real plotting
                    for (int q=0; q<plotPts; q++) plotX[q] *= angularRes.get(i).get(g);

                    // detection
                    Plot p = new Plot("", "orient.[deg]", profilesName.get(i).get(g), plotX, plotY);
					p.setSize(600, 300);
					p.draw();

                    /*

                    // add convergence points
                    float[] plotConvX = new float[Analyzer.nrPoints];
                    float[] plotConvY = new float[Analyzer.nrPoints];

                    for (int i1=0; i1<plotConvX.length; i1++)
                        plotConvX[i1] = Analyzer.convIdx.get(i).get(g)[i1] * angularRes.get(i).get(g);

                    for (int i1=0; i1<plotConvY.length; i1++)
                        plotConvY[i1] = (float) Tools.interp1Darray(Analyzer.convIdx.get(i).get(g)[i1], profiles.get(i).get(g));

                    p.addPoints(plotConvX, plotConvY, Plot.X);
*/
                    // backgr
                    float[] currBckgr = new float[plotPts];
                    for (int i1=0; i1<plotPts; i1++) currBckgr[i1] = bacgr.get(i);
                    p.setColor(Color.BLUE);
                    p.addPoints(plotX, currBckgr, Plot.LINE);

                    // plot detections if there were any
					if (idxs.get(i).get(g)!=null) {

                         /*
                            i - location, g - configuration
                          */

                            float pX, pY;
                            pX = angles.get(i).get(g)[0];
                            pY = peaks.get(i).get(g)[0];
                            p.setColor(Color.RED);
                            p.addPoints(new float[]{pX}, new float[]{pY}, Plot.BOX);
                            p.draw();

                            pX = angles.get(i).get(g)[1];
                            pY = peaks.get(i).get(g)[1];
                            p.setColor(Color.GREEN);
                            p.addPoints(new float[]{pX}, new float[]{pY}, Plot.BOX);
                            p.draw();

                            pX = angles.get(i).get(g)[2];
                            pY = peaks.get(i).get(g)[2];
                            p.setColor(Color.BLUE);
                            p.addPoints(new float[]{pX}, new float[]{pY}, Plot.BOX);
                            p.draw();

                        double rd = neuronDiamMax*scales[g]; // plot radius

                        float pointX, pointY;
                        PointRoi point;
                        // cluster 0 - RED
                        // i - location, g-configuration, 0-cluster
                        pointX = (float) (atX+rd*Math.cos(angles.get(i).get(g)[0]*(2*Math.PI/360))+.5);
                        pointY = (float) (atY-rd*Math.sin(angles.get(i).get(g)[0]*(2*Math.PI/360))+.5);
                        point = new PointRoi(pointX, pointY);
                        point.setStrokeColor(Color.RED);
                        ov.add(point);

                        // cluster 1 - GREEN
                        // i - location, g-configuration
                        pointX = (float) (atX+rd*Math.cos(angles.get(i).get(g)[1]*(2*Math.PI/360))+.5);
                        pointY = (float) (atY-rd*Math.sin(angles.get(i).get(g)[1]*(2*Math.PI/360))+.5);
                        point = new PointRoi(pointX, pointY);
                        point.setStrokeColor(Color.GREEN);
                        ov.add(point);

                        // cluster 2 - BLUE
                        // i - location, g-configuration,
                        pointX = (float) (atX+rd*Math.cos(angles.get(i).get(g)[2]*(2*Math.PI/360))+.5);
                        pointY = (float) (atY-rd*Math.sin(angles.get(i).get(g)[2]*(2*Math.PI/360))+.5);
                        point = new PointRoi(pointX, pointY);
                        point.setStrokeColor(Color.BLUE);
                        ov.add(point);

                        ov.add(new OvalRoi(atX-rd+.5, atY-rd+.5, 2*rd, 2*rd));

                    }

                    inimg.setOverlay(ov);
                    vizProfileStack.addSlice("", p.getProcessor());

                    /*
                    export R format
                     */

                    String printProfile = "";
                    String printAngle = "";
                    String xName = "ang_"+profilesName.get(i).get(g);
                    String yName = profilesName.get(i).get(g);

                    printProfile    += yName+" <- c(";
                    printAngle      += xName+" <- c(";

                    for (int i1=0; i1<plotPts; i1++) {

                        printProfile+=plotY[i1] +"";// profiles.get(i).get(g)
                        printAngle  +=plotX[i1] +"";

                        if(i1<plotPts-1) {
                            printProfile    +=", ";
                            printAngle      +=", ";
                        }
                    }

                    printProfile+=")\n";
                    printAngle+=")\n";

                    try {
                        PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(fileName, true)));
                        out.println(printAngle);
                        out.println(printProfile);
                        if (g==0)
                            out.println("plot("+xName+","+yName+",type=\"l\")\n grid()");
                        else
                            out.println("lines("+xName+","+yName+",type=\"l\")");
                        out.close();

                    } catch (IOException e1) {}

                }

            }
        }

        if (vizProfileStack.getSize()>0) {
            vizProfileImage.setStack(vizProfileStack);
            vizProfileImage.updateAndDraw();
            vizProfileImage.setTitle("profiles");
            vizProfileImage.show();
        }

    }

    public void mousePressed(MouseEvent e) {}

    public void mouseReleased(MouseEvent e) {}

    public void mouseEntered(MouseEvent e) {}

    public void mouseExited(MouseEvent e) {}
}
