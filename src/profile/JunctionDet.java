package profile;

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

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 7/7/13
 * Time: 8:13 PM
 */
public class JunctionDet implements PlugInFilter, MouseListener {

	ImagePlus 	    inimg;
    ByteProcessor  scoreimg;

    // visualization  (for selected locations)
    ImagePlus       vizProfileImage;
    ImageStack      vizProfileStack;

    ArrayList<ArrayList<float[]>>   profiles;
    ArrayList<ArrayList<Integer>>   angularRes;
    ArrayList<ArrayList<String>>    profilesName;

	//ArrayList<ArrayList<float[]>>   hdomes;

    ArrayList<Double>               neuronD;    // per cfg
    ArrayList<Double>               scale;      // per cfg

	// ms output
    ArrayList<ArrayList<float[]>>   angles;   	// 3 angles  (necessary for matching)
    ArrayList<ArrayList<float[]>>   idxs;   	// 3 profile indexes
    ArrayList<ArrayList<float[]>>   peaks;      // each angle one peak score

    ArrayList<Float>   				bacgr;      // per loc

    ArrayList<ArrayList<Float>>   	centrAvg;

	Overlay detections;

    /*
    parameters
     */
    double neuronDiamMax, D, eps, dMin; // scaleMin, scaleMax, neuronDiamMin
    int H, M, MaxIter, CPU_NR;

	public int setup(String s, ImagePlus imagePlus) {

		if(imagePlus==null) return DONE;
		inimg = Tools.convertToFloatImage(imagePlus);
		inimg.setTitle("inimg");
		return DOES_8G+DOES_32+NO_CHANGES;

	}

	public void run(ImageProcessor imageProcessor) {

        int totalJobs; // used in paralelization
        long t1, t2;

		// fix scales
		double[] scales = new double[3];
		scales[0] = 1.5;
		scales[1] = 2.0;
		scales[2] = 2.5;


        /*
        generic dialog
        */
        //neuronDiamMin       =  Prefs.get("advantra.critpoint.neuronDiamMin", 3);
        neuronDiamMax       =  Prefs.get("advantra.critpoint.neuronDiamMax", 3);
        //scaleMin            =  Prefs.get("advantra.critpoint.scaleMin", 1.5);
        //scaleMax            =  Prefs.get("advantra.critpoint.scaleMax", 2.5);
        D                   =  Prefs.get("advantra.critpoint.D", 10);

        H       = (int) Prefs.get("advantra.critpoint.ms2d.H", 4);
        MaxIter = (int) Prefs.get("advantra.critpoint.ms2d.MaxIter", 200);
        eps     =       Prefs.get("advantra.critpoint.ms2d.eps", 0.0001);
        dMin       =       Prefs.get("advantra.critpoint.ms2d.d", 0.5);
        M       = (int) Prefs.get("advantra.critpoint.ms2d.M", 5);


        CPU_NR              = (int) Prefs.get("advantra.critpoint.CPU_NR", 4);

        GenericDialog gd = new GenericDialog("JUNCTION DET.");
        //gd.addNumericField("neuronD ",  neuronDiamMin, 0, 10, " MIN");
        gd.addNumericField("neuronD ",  neuronDiamMax, 0, 10, " MAX");
        //gd.addNumericField("scale",     scaleMin, 1, 10, "MIN");
        //gd.addNumericField("scale",     scaleMax, 1, 10, "MAX");
        gd.addNumericField("D",         D, 1, 10, "score param.");

		/*gd.addMessage("--- LOCAL PEAKS (MS on 2D image) ---");
        gd.addNumericField("H ",        H,      0, 10, "spatial neighbourhood");
        gd.addNumericField("Max Iter",  MaxIter,0, 10, "max #iterations");
        gd.addNumericField("Eps",       eps,    6, 10, "convergence threshold");
        gd.addMessage("---");
        gd.addNumericField("d",         dMin,      1, 10, "min intra-cluster distance");
        gd.addNumericField("M",         M,      0, 10, "min # points in cluster");
        */

		gd.addMessage("---");
        gd.addNumericField("CPU_NR ",   CPU_NR, 0, 10, "");

        gd.showDialog();
        if (gd.wasCanceled()) return;

//        neuronDiamMin       =  gd.getNextNumber();
//        Prefs.set("advantra.critpoint.neuronDiamMin", 	neuronDiamMin);
        neuronDiamMax       =  gd.getNextNumber();
        Prefs.set("advantra.critpoint.neuronDiamMax",   neuronDiamMax);
//        scaleMin            =  gd.getNextNumber();
//        Prefs.set("advantra.critpoint.scaleMin", 	    scaleMin);
//        scaleMax            =  gd.getNextNumber();
//        Prefs.set("advantra.critpoint.scaleMax",        scaleMax);
        D                   =  gd.getNextNumber();
        Prefs.set("advantra.critpoint.D",               D);

		/*
		H       =  (int) gd.getNextNumber();
        Prefs.set("advantra.critpoint.ms2d.H", 	H);
		MaxIter =  (int) gd.getNextNumber();
        Prefs.set("advantra.critpoint.ms2d.MaxIter", MaxIter);
        eps     =  gd.getNextNumber();
        Prefs.set("advantra.critpoint.ms2d.eps", eps);
        dMin       =  gd.getNextNumber();
        Prefs.set("advantra.critpoint.ms2d.d", dMin);
        M       =  (int) gd.getNextNumber();
        Prefs.set("advantra.critpoint.ms2d.M", M);
        */


		CPU_NR              =  (int) gd.getNextNumber();
        Prefs.set("advantra.critpoint.CPU_NR", 	        CPU_NR);

        /*
        calculate the mask  (MaskerDemo example), should result in
        ByteProcessor inmask  == Masker.maskip
        */
        IJ.log("extracting background...");
        t1 = System.currentTimeMillis();
        int neighbourhoodR = (int) Math.ceil(4*neuronDiamMax);
        Masker.loadTemplate(inimg.getProcessor(), neighbourhoodR, (float) D);
        totalJobs = Masker.image_height*Masker.image_width;
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

		/*
		// overlay mask
		ImagePlus maskViz = new ImagePlus("inmask", Masker.maskip);
        maskViz.show();
        IJ.selectWindow("inimg");
        IJ.run("Add Image...", "image=inmask x="+0+" y="+0+" opacity=30");
        maskViz.close();
		*/


        /*
        score image initialize
         */
        scoreimg = new ByteProcessor(inimg.getWidth(), inimg.getHeight());
		Profiler.loadTemplate(inimg.getProcessor(), Masker.maskip);
        // can see the locations now
        int totalLocations = Profiler.locations.length;

        IJ.log("total "+totalLocations+" locations given by mask ("+(100*(float)totalLocations/(inimg.getWidth()*inimg.getHeight()))+" percent kept)");

        profiles     	= new ArrayList<ArrayList<float[]>>(totalLocations);
		angularRes    	= new ArrayList<ArrayList<Integer>>(totalLocations);
        profilesName  	= new ArrayList<ArrayList<String>>(totalLocations);

        neuronD = new ArrayList<Double>(); // this is actually per configuration
        scale 	= new ArrayList<Double>();

		IJ.log("calculating profiles... ");
		t1 = System.currentTimeMillis();

		int totalConfs = 0;
        //for (double neuronDiam = neuronDiamMax; neuronDiam<=neuronDiamMax; neuronDiam++) {

			double neuronDiam = neuronDiamMax;

            for (int scaleIdx=0; scaleIdx<scales.length; scaleIdx++) {

				double sc = scales[scaleIdx];
				totalConfs++;

                neuronD.add(neuronDiam);
                scale.add(sc);

                Profiler.loadParams(neuronDiam, sc);

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

                updateList();

            }
        //}
		t2 = System.currentTimeMillis();
		IJ.log("done "+((t2-t1)/1000f)+" sec.");

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

        IJ.log("finding matches and corresponding scores for each location...");

        angles    	= new ArrayList<ArrayList<float[]>>(totalLocations);
        idxs    	= new ArrayList<ArrayList<float[]>>(totalLocations);
        peaks       = new ArrayList<ArrayList<float[]>>(totalLocations);
		centrAvg   	= new ArrayList<ArrayList<Float>>(totalLocations);
		bacgr       = new ArrayList<Float>(totalLocations);

        for (int locIdx=0; locIdx<totalLocations; locIdx++) {

            int atX = Profiler.locations[locIdx][0];
            int atY = Profiler.locations[locIdx][1];

            ArrayList<float[]>  anglesPerCfg = new ArrayList<float[]>(totalConfs);
            ArrayList<float[]>  idxsPerCfg = new ArrayList<float[]>(totalConfs);
            ArrayList<float[]>  peaksPerCfg = new ArrayList<float[]>(totalConfs);
            ArrayList<Float>   centralAvgPerCfg = new ArrayList<Float>(totalConfs);

            boolean isFirst = true;
            float[] refAnglesDeg    = new float[3];

            for (int cfgIdx=0; cfgIdx<totalConfs; cfgIdx++) {

                /*
                 find central bloc score for this location
                  */

                centralAvgPerCfg.add(centralAvg(atX, atY, neuronD.get(cfgIdx), (FloatProcessor) inimg.getProcessor()));

                // get indexes form Analyzer for this conf & loc
                if (Analyzer.peakIdx[locIdx][cfgIdx]!=null) {

                    /*
                    converging angles and scores
                     */

                    float[] anglesDeg       = new float[3];
                    float[] indexes       	= new float[3];
                    float[] peaks          	= new float[3];

                    indexes[0] = Analyzer.peakIdx[locIdx][cfgIdx][0];
                    indexes[1] = Analyzer.peakIdx[locIdx][cfgIdx][1];
                    indexes[2] = Analyzer.peakIdx[locIdx][cfgIdx][2];

                    anglesDeg[0] = Analyzer.peakIdx[locIdx][cfgIdx][0] * angularRes.get(locIdx).get(cfgIdx);
                    anglesDeg[1] = Analyzer.peakIdx[locIdx][cfgIdx][1] * angularRes.get(locIdx).get(cfgIdx);
                    anglesDeg[2] = Analyzer.peakIdx[locIdx][cfgIdx][2] * angularRes.get(locIdx).get(cfgIdx);

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

            centrAvg.add(centralAvgPerCfg);

            bacgr.add(Masker.back.getf(atX, atY)); // take the value from Masker

        }

        IJ.log("done, extracting the scores...");

        // loop once again and extract scores

        for (int locIdx=0; locIdx<totalLocations; locIdx++) {

			double score =  extractScore(locIdx, false);

            scoreimg.set(Profiler.locations[locIdx][0], Profiler.locations[locIdx][1], ((int) score));

        }


        ImagePlus scoreImagePlus = new ImagePlus("score", scoreimg);
		scoreImagePlus.show();
		//IJ.run(scoreImagePlus, "Smooth", "");

/*		IJ.log("detecting local peaks in score...");
        t1= System.currentTimeMillis();
        MS2D.loadTemplate(scoreimg, Masker.maskip, H, MaxIter, eps);
        totalJobs = MS2D.toProcess;
        MS2D ms_jobs1[] = new MS2D[CPU_NR];
        for (int i = 0; i < ms_jobs1.length; i++) {
            ms_jobs1[i] = new MS2D(i*totalJobs/CPU_NR,  (i+1)*totalJobs/CPU_NR);
            ms_jobs1[i].start();
        }
        for (int i = 0; i < ms_jobs1.length; i++) {
            try {
                ms_jobs1[i].join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
        t2 = System.currentTimeMillis();
        IJ.log("done extracting peaks in score "+((t2-t1)/1000f)+" sec.");


        IJ.log("clustering...");
        t1=System.currentTimeMillis();
        ArrayList<ArrayList<double[]>> clusters = MS2D.extractClusters(dMin);
        t2=System.currentTimeMillis();
        IJ.log("done, "+clusters.size()+" clusters " + ((t2 - t1) / 1000f) + " s.");

        IJ.log("cluster centroids...");
        t1=System.currentTimeMillis();
        ArrayList<double[]> centroids = MS2D.extractClusterCentroids(clusters, M);
        t2=System.currentTimeMillis();
        IJ.log("done " + ((t2 - t1) / 1000f) + " sec.");*/

        detections = new Overlay();

		for (int i1=0; i1<scoreimg.getWidth(); i1++) {
			for (int i2=0; i2<scoreimg.getHeight(); i2++) {
				if (scoreimg.get(i1, i2)==255) {
					PointRoi pt = new PointRoi(i1+0.5, i2+0.5);
					detections.add(pt);
				}

			}
		}

/*        for (int i1=0; i1<centroids.size(); i1++) {

            double atX = centroids.get(i1)[0];
            double atY = centroids.get(i1)[1];

            //double rds = Math.sqrt(centroids.get(i1)[2]/Math.PI);
            double rds = 1;//centroids.get(i1)[3]/0.01f;

            OvalRoi oval = new OvalRoi(atX+0.5-rds, atY+0.5-rds, 2*rds, 2*rds);
            oval.setStrokeColor(Color.RED);
			detections.add(oval);

        }*/

//        for (int i1=0; i1<Tnew.length; i1++) {
//            PointRoi pt = new PointRoi(Tnew[i1][1]+0.5, Tnew[i1][0]+0.5);
//            pt.setStrokeColor(Color.PINK);
//            ov.add(pt);
//        }

		inimg.setOverlay(detections);
        scoreImagePlus.setOverlay(detections);

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
//        inimg.getCanvas().zoomIn(0, 0);

	}

	public double extractScore(int locIdx, boolean printVals) {

		double A0 = centrAvg.get(locIdx).get(0) - bacgr.get(locIdx);
		if (printVals) IJ.log("A0 = "+A0);
		A0 = (A0>0)? 1 : 0;//A0
		if (printVals) IJ.log("A0 = "+A0);
		//A0 = 1-Math.exp(-A0/D);

		double A1 = 0;// per cluster
		double A2 = 0;
		double A3 = 0;

		double A1P1x = 0;
		double A1P1y = 0;
		double A2P1x = 0;
		double A2P1y = 0;
		double A3P1x = 0;
		double A3P1y = 0;

		double A1P2x = 0;
		double A1P2y = 0;
		double A2P2x = 0;
		double A2P2y = 0;
		double A3P2x = 0;
		double A3P2y = 0;

		double A1P3x = 0;
		double A1P3y = 0;
		double A2P3x = 0;
		double A2P3y = 0;
		double A3P3x = 0;
		double A3P3y = 0;

		int cntCfgs = 0;

		for (int cfgIdx=0; cfgIdx<3; cfgIdx++) {

			if (angles.get(locIdx).get(cfgIdx)!=null) { // if ms gave 3 values

				// there are 3 values - check validity for each configuration

				if (cfgIdx==2 && !isProfileValid(idxs.get(locIdx).get(cfgIdx), profiles.get(locIdx).get(cfgIdx), bacgr.get(locIdx)+1)) {
					if (printVals) IJ.log("profile was not valid!");
					break;
				}

				if (cfgIdx==0) {
					double rd = neuronD.get(0)*scale.get(0);
					A1P1x = rd*Math.cos(angles.get(locIdx).get(0)[0]*(2*Math.PI/360f));
					A1P1y = -rd*Math.sin(angles.get(locIdx).get(0)[0]*(2*Math.PI/360f));
					A2P1x = rd*Math.cos(angles.get(locIdx).get(0)[1]*(2*Math.PI/360f));
					A2P1y = -rd*Math.sin(angles.get(locIdx).get(0)[1]*(2*Math.PI/360f));
					A3P1x = rd*Math.cos(angles.get(locIdx).get(0)[2]*(2*Math.PI/360f));
					A3P1y = -rd*Math.sin(angles.get(locIdx).get(0)[2]*(2*Math.PI/360f));

					double A1t = peaks.get(locIdx).get(0)[0]- bacgr.get(locIdx);
					if (printVals) IJ.log("\ndiff: "+peaks.get(locIdx).get(0)[0]+", "+bacgr.get(locIdx)+" , "+A1t);
					A1t = (A1t>D)? 1 : 0 ; // 1-Math.exp(-A1t/D)
					if (printVals) IJ.log("A1t: "+A1t);
					A1 +=  A1t;

					double A2t = peaks.get(locIdx).get(0)[1] - bacgr.get(locIdx);
					if (printVals) IJ.log("\ndiff: "+peaks.get(locIdx).get(0)[1]+", "+bacgr.get(locIdx)+" , "+A2t);
					A2t = (A2t>D)? 1 : 0 ; // 1-Math.exp(-A2t/D)
					if (printVals) IJ.log("A2t: "+A2t);
					A2 += A2t;

					double A3t = peaks.get(locIdx).get(0)[2] - bacgr.get(locIdx);
					if (printVals) IJ.log("\ndiff: "+peaks.get(locIdx).get(0)[2]+", "+bacgr.get(locIdx)+" , "+A3t);
					A3t = (A3t>D)? 1 : 0; // 1-Math.exp(-A3t/D)
					if (printVals) IJ.log("A3t: "+A3t);
					A3 += A3t;
				}

				if (cfgIdx==1) {
					double rd = neuronD.get(1)*scale.get(1);
					A1P2x = rd*Math.cos(angles.get(locIdx).get(1)[0]*(2*Math.PI/360f));
					A1P2y = -rd*Math.sin(angles.get(locIdx).get(1)[0]*(2*Math.PI/360f));
					A2P2x = rd*Math.cos(angles.get(locIdx).get(1)[1]*(2*Math.PI/360f));
					A2P2y = -rd*Math.sin(angles.get(locIdx).get(1)[1]*(2*Math.PI/360f));
					A3P2x = rd*Math.cos(angles.get(locIdx).get(1)[2]*(2*Math.PI/360f));
					A3P2y = -rd*Math.sin(angles.get(locIdx).get(1)[2]*(2*Math.PI/360f));

					double A1t = peaks.get(locIdx).get(1)[0]- bacgr.get(locIdx);
					if (printVals) IJ.log("\ndiff: "+peaks.get(locIdx).get(1)[0]+", "+bacgr.get(locIdx)+" , "+A1t);
					A1t = (A1t>D)? 1 : 0 ; // 1-Math.exp(-A1t/D)
					if (printVals) IJ.log("A1t: "+A1t);
					A1 +=  A1t;

					double A2t = peaks.get(locIdx).get(1)[1] - bacgr.get(locIdx);
					if (printVals) IJ.log("\ndiff: "+peaks.get(locIdx).get(1)[1]+", "+bacgr.get(locIdx)+" , "+A2t);
					A2t = (A2t>D)? 1 : 0 ; // 1-Math.exp(-A2t/D)
					if (printVals) IJ.log("A2t: "+A2t);
					A2 += A2t;

					double A3t = peaks.get(locIdx).get(1)[2] - bacgr.get(locIdx);
					if (printVals) IJ.log("\ndiff: "+peaks.get(locIdx).get(1)[2]+", "+bacgr.get(locIdx)+" , "+A3t);
					A3t = (A3t>D)? 1 : 0; // 1-Math.exp(-A3t/D)
					if (printVals) IJ.log("A3t: "+A3t);
					A3 += A3t;
				}

				if (cfgIdx==2) {
					double rd = neuronD.get(2)*scale.get(2);
					A1P3x = rd*Math.cos(angles.get(locIdx).get(2)[0]*(2*Math.PI/360f));
					A1P3y = -rd*Math.sin(angles.get(locIdx).get(2)[0]*(2*Math.PI/360f));
					A2P3x = rd*Math.cos(angles.get(locIdx).get(2)[1]*(2*Math.PI/360f));
					A2P3y = -rd*Math.sin(angles.get(locIdx).get(2)[1]*(2*Math.PI/360f));
					A3P3x = rd*Math.cos(angles.get(locIdx).get(2)[2]*(2*Math.PI/360f));
					A3P3y = -rd*Math.sin(angles.get(locIdx).get(2)[2]*(2*Math.PI/360f));

					double A1t = peaks.get(locIdx).get(2)[0]- bacgr.get(locIdx);
					if (printVals) IJ.log("\ndiff: "+peaks.get(locIdx).get(2)[0]+", "+bacgr.get(locIdx)+" , "+A1t);
					A1t = (A1t>D)? 1 : 0 ; // 1-Math.exp(-A1t/D)
					if (printVals) IJ.log("A1t: "+A1t);
					A1 +=  A1t;

					double A2t = peaks.get(locIdx).get(2)[1] - bacgr.get(locIdx);
					if (printVals) IJ.log("\ndiff: "+peaks.get(locIdx).get(2)[1]+", "+bacgr.get(locIdx)+" , "+A1t);
					A2t = (A2t>D)? 1 : 0 ; // 1-Math.exp(-A2t/D)
					if (printVals) IJ.log("A2t: "+A2t);
					A2 += A2t;

					double A3t = peaks.get(locIdx).get(2)[2] - bacgr.get(locIdx);
					if (printVals) IJ.log("\ndiff: "+peaks.get(locIdx).get(2)[2]+", "+bacgr.get(locIdx)+" , "+A1t);
					A3t = (A3t>D)? 1 : 0; // 1-Math.exp(-A3t/D)
					if (printVals) IJ.log("A3t: "+A3t);
					A3 += A3t;
				}

				cntCfgs++;

			}
			else {
				break;
			}

		}

		double score = 0;

		if (cntCfgs==3) {

			//A1 /= cntCfgs;
			double A1v1x, A1v1y, iA1v1, A1v2x, A1v2y, iA1v2;
			A1v1x = A1P2x-A1P1x;
			A1v1y = A1P2y-A1P1y;
			A1v2x = A1P3x-A1P2x;
			A1v2y = A1P3y-A1P2y;
			iA1v2 = Math.sqrt(A1v2x*A1v2x+A1v2y*A1v2y);
			iA1v1 = Math.sqrt(A1v1x*A1v1x+A1v1y*A1v1y);
			double cosA1 = ( (A1v1x*A1v2x+A1v1y*A1v2y) / (iA1v1 * iA1v2) );
			double align1 = ( (cosA1+1) / 2 );
			if (printVals) IJ.log("A1:" + A1 + " , how aligned? " + align1+" angle: "+ (Math.acos(cosA1)*(180f/Math.PI)) );
				/*

				 */

			//A2 /= cntCfgs;
			double A2v1x, A2v1y, iA2v1, A2v2x, A2v2y, iA2v2;
			A2v1x = A2P2x-A2P1x;
			A2v1y = A2P2y-A2P1y;
			A2v2x = A2P3x-A2P2x;
			A2v2y = A2P3y-A2P2y;
			iA2v2 = Math.sqrt(A2v2x*A2v2x+A2v2y*A2v2y);
			iA2v1 = Math.sqrt(A2v1x*A2v1x+A2v1y*A2v1y);
			double cosA2 = ((A2v1x*A2v2x+A2v1y*A2v2y) / (iA2v1 * iA2v2));
			double align2 = ( (cosA2+1) / 2 );
			if (printVals) IJ.log("A2:" + A2 + " , how aligned? " + align2+" angle: "+ (Math.acos(cosA2)*(180f/Math.PI)));
				/*

				 */

			//A3 /= cntCfgs;
			double A3v1x, A3v1y, iA3v1, A3v2x, A3v2y, iA3v2;
			A3v1x = A3P2x-A3P1x;
			A3v1y = A3P2y-A3P1y;
			A3v2x = A3P3x-A3P2x;
			A3v2y = A3P3y-A3P2y;
			iA3v2 = Math.sqrt(A3v2x*A3v2x+A3v2y*A3v2y);
			iA3v1 = Math.sqrt(A3v1x*A3v1x+A3v1y*A3v1y);
			double cosA3 = ((A3v1x*A3v2x+A3v1y*A3v2y) / (iA3v1 * iA3v2));
			double align3 = ( (cosA3+1) / 2 );
			if (printVals) IJ.log("A3:" + A3 + " , how aligned? " + align3+" angle: "+ (Math.acos(cosA3)*(180f/Math.PI)));

			//score = A0 * A1 * align1 * A2 * align2 * A3 * align3;//
			//score = A0 * Tools.min3(A1 * align1, A2 * align2, A3 * align3);
			score = ((A0 + A1 + A2 + A3)==10)? 255: 0;
			if (printVals) IJ.log("SCORE = "+score);

		}
		else {
			if (printVals) IJ.log("not enough points");
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

    private void updateList() // uses Profiler to update  profiles, profilesName, angularRes
    {
        // store it in the list for each loc
        for (int loopLocations=0; loopLocations<Profiler.locations.length; loopLocations++) {

            // profileToAdd from Profiler.profiles
            float[] profileToAdd = new float[Profiler.profiles[loopLocations].length];
            for (int k=0; k<profileToAdd.length; k++) {
                profileToAdd[k] = Profiler.profiles[loopLocations][k];
            }

            // stringToAdd from Profiler.neuronDiam, Profiler.scale, Profiler.resolDeg
            String stringToAdd = "nD_"+Profiler.neuronDiam+"_s_"+Profiler.scale+"_a_"+Profiler.resolDeg;

            // intToAdd from Profiler.resolDeg
            int intToAdd = Profiler.resolDeg;

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

		boolean profileValid = true;

		for (int i=0; i<peakIndexes.length; i++) {

			boolean intervalOK = false;

			int startInterval = Math.round(peakIndexes[i]);
			startInterval = (startInterval<0)? 0 : startInterval;
			startInterval = (startInterval>=profileToCheck.length)? profileToCheck.length-1 : startInterval;

			int endInterval = Math.round( (i==peakIndexes.length-1)? peakIndexes[0] : peakIndexes[i+1] );
			endInterval = (endInterval<0)? 0 : endInterval;
			endInterval = (endInterval>=profileToCheck.length)? profileToCheck.length-1 : endInterval ;

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

                // pick all the profiles from the list for that location

                int nrProfiles = profiles.get(i).size();    // this corresponds to nrCfgs

                /*
				overlay for the main image
			    */
                Overlay ov = detections.duplicate();//new Overlay();
                PointRoi pt = new PointRoi(atX+.5, atY+.5);
                ov.add(pt);

                for (int g=0; g<nrProfiles; g++) {  // g will loop profiles

                    int currProfileLen = profiles.get(i).get(g).length;

//                    // allocate profile median
//                    float[] profileMed = new float[currProfileLen];
//                    for (int i1=0; i1<currProfileLen; i1++) profileMed[i1] = backgrPerLocation.get(i).get(g).floatValue();

                    /*
                    plot ij  template
                     */
                    float[] angIdx = new float[currProfileLen];

                    for (int q=0; q<currProfileLen; q++) {
                        angIdx[q] = q * angularRes.get(i).get(g);
                    }

                    /*
                    plot def.
                     */

                    // profile
                    Plot p = new Plot("", "orient.[deg]", profilesName.get(i).get(g), angIdx, profiles.get(i).get(g));
					p.setSize(600, 300);
					p.draw();

                    // here add experimental H-dome profile right below just for this location
//                    // binarize it & plot those angular values that satisfy the cirteria
//                    for (int i1=0; i1<profilesPerLocation.get(i).get(g).length; i1++) {
//                        if (hdomeProfile[i1]>0.8*H) {
//
//                        }
//                    }

                    /*
                    p = new Plot("", "", "H-dome,H="+H, angIdx, hdomeProfile); // profilesPerLocation.get(i).get(g)
                    p.setSize(600, 300);
                    p.setColor(Color.CYAN);// need to make it colour to align with previous in the same stack
                    //p.addPoints(angIdx, , Plot.LINE);
                    p.draw();

                    */

                    // add convergence points
                    float[] plotConvX = new float[Analyzer.nrPoints];
                    float[] plotConvY = new float[Analyzer.nrPoints];

                    for (int i1=0; i1<plotConvX.length; i1++)
                        plotConvX[i1] = Analyzer.convIdx.get(i).get(g)[i1] * angularRes.get(i).get(g);

                    for (int i1=0; i1<plotConvY.length; i1++)
                        plotConvY[i1] = (float) Tools.interp1Darray(Analyzer.convIdx.get(i).get(g)[i1], profiles.get(i).get(g));

                    p.addPoints(plotConvX, plotConvY, Plot.X);

                    // add h-dome plot

                    float profMax = Float.MIN_VALUE;
                    float profMin = Float.MAX_VALUE;
                    for (int j1=0; j1<profiles.get(i).get(g).length; j1++) {

                        if (profiles.get(i).get(g)[j1]>profMax) profMax = profiles.get(i).get(g)[j1];

                        if (profiles.get(i).get(g)[j1]<profMin)	profMin = profiles.get(i).get(g)[j1];

                    }

/*                    // plot hdomes
                    float[] plotHdomeX = new float[hdomes.get(i).get(g).length];
                    float[] plotHdomeY = new float[plotHdomeX.length];

                    for (int i1=0; i1<plotHdomeX.length; i1++)
                        plotHdomeX[i1] = (float) i1 * angularRes.get(i).get(g);

                    for (int i1=0; i1<plotHdomeY.length; i1++)
                        plotHdomeY[i1] = (hdomes.get(i).get(g)[i1]>0)? 0.9f*profMax : 1.1f*profMin ;

                    p.setColor(Color.GREEN);
                    p.addPoints(plotHdomeX, plotHdomeY, Plot.LINE);
                    p.draw();*/

                    // backgr
                    float[] currBckgr = new float[currProfileLen];
                    for (int i1=0; i1<currProfileLen; i1++) currBckgr[i1] = bacgr.get(i);
                    p.setColor(Color.BLUE);
                    p.addPoints(angIdx, currBckgr, Plot.LINE);
//                    p.setColor(Color.GREEN);
//                    p.addPoints(angIdx, profileMed, Plot.LINE);

                    // plot dtections if there were any
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

                        // modify xAxis angles into euclidean points
                        double[][] msPeaks = new double[3][2];
                        double rd = neuronD.get(g)*scale.get(g); // plot radius

//                        float[] xAxis = new float[3];
//                        for (int i1=0; i1<3; i1++) xAxis[i1] = Analyzer.peakIdx[i][g][i1] * angularRes.get(i).get(g); // this is precalculated now

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
                    export R
                     */

                    String printProfile = "";
                    String printAngle = "";
                    String xName = "ang_"+profiles.get(i).get(g);
                    String yName = profilesName.get(i).get(g);

                    printProfile    += yName+" <- c(";
                    printAngle      += xName+" <- c(";

                    for (int i1=0; i1<currProfileLen; i1++) {

                        printProfile+=profiles.get(i).get(g)[i1]+"";
                        printAngle+=angIdx[i1]+"";

                        if(i1<profiles.get(i).get(g).length-1) {
                            printProfile+=", ";
                            printAngle+=", ";
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

                /*
                log IJ
                */

                //double D = 10;

                System.out.println("\n---\nX: "+atX+" , Y: "+atY+"");

				extractScore(i, true);

                System.out.println("\n*** CENTER  ***");
                double pk = centrAvg.get(i).get(0);
                double bg = bacgr.get(i);
                double df = (pk-bg>0) ? pk-bg : 0 ;

                System.out.println(" peak: "+pk+" ("+bg+") A0="+IJ.d2s(1-Math.exp(-df/D), 2)+" | ");

                for (int clstIdx=0; clstIdx<3; clstIdx++) {
                    System.out.println("\n*** CLUSTER "+clstIdx+" ***");

                    System.out.print(" peaks: ");
                    for (int cfgIdx=0; cfgIdx<3; cfgIdx++) {
                        if (peaks.get(i).get(cfgIdx)!=null) {
                            double pk1 = peaks.get(i).get(cfgIdx)[clstIdx];
                            double bg1 = bacgr.get(i);
                            double df1 = (pk1-bg1>0) ? pk1-bg1 : 0 ;
                            System.out.print(pk1+"("+bg1+") -> "+IJ.d2s(1-Math.exp(-df1/D), 2)+" | ");
                        }
                    }

                    System.out.print(" angles: ");
                    for (int cfgIdx=0; cfgIdx<3; cfgIdx++) {
                        if (angles.get(i).get(cfgIdx)!=null)
                            System.out.print(angles.get(i).get(cfgIdx)[clstIdx] + " |  ");
                    }

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
