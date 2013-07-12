package profile;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
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

	ImagePlus 	inimg;
	String 		inimgPath;
	ImagePlus   inmask;

    FloatProcessor  scoreimg;

    // visualization
    ImagePlus       vizProfileImage;
    ImageStack      vizProfileStack;

    ArrayList<ArrayList<float[]>>   profiles;
    ArrayList<ArrayList<Integer>>   angularRes;
    ArrayList<ArrayList<String>>    profilesName;

	ArrayList<ArrayList<float[]>>   hdomes;

    ArrayList<Double>               neuronD; // per cfg
    ArrayList<Double>               scale;   // per cfg

	// ms output
    ArrayList<ArrayList<float[]>>   angles;   	// 3 angles  (necessary for matching)
    ArrayList<ArrayList<float[]>>   idxs;   	// 3 profile indexes
    ArrayList<ArrayList<float[]>>   peaks;      // each angle one peak score

    ArrayList<Float>   				bacgr;    // per loc

    ArrayList<ArrayList<Float>>   	centrAvg;

    int CPU_NR;

	public int setup(String s, ImagePlus imagePlus) {

		if(imagePlus==null) {
			IJ.showMessage("needs image opened"); return DONE; }
		inimg = Tools.convertToFloatImage(imagePlus);
		inimg.setTitle("inimg");
		inimgPath = imagePlus.getOriginalFileInfo().directory+imagePlus.getOriginalFileInfo().fileName;
		return DOES_8G+DOES_32+NO_CHANGES;

	}

	public void run(ImageProcessor imageProcessor) {

		String inmaskPath = Tools.removeExtension(inimgPath)+".mask";

		if (new File(inmaskPath).exists()) {
			inmask = new ImagePlus(inmaskPath);
			inmask.setTitle("inmask");
		}
		else {
			byte[] array = new byte[inimg.getWidth() * inimg.getHeight()];
			for (int a = 0; a < array.length; a++) array[a] = (byte) 255;
			inmask = new ImagePlus("inmask", new ByteProcessor(inimg.getWidth(), inimg.getHeight(), array));
		}

        /*
        score image initialize
         */
        scoreimg = new FloatProcessor(inimg.getWidth(), inimg.getHeight());

        CPU_NR = 5;  // GD

		Profiler.loadTemplate(inimg.getProcessor(), (ByteProcessor) inmask.getProcessor());

        int totalLocations = Profiler.locations.length;

        profiles     	= new ArrayList<ArrayList<float[]>>(totalLocations);
		angularRes    	= new ArrayList<ArrayList<Integer>>(totalLocations);
        profilesName  	= new ArrayList<ArrayList<String>>(totalLocations);

        neuronD = new ArrayList<Double>(); // this is actually per configuration
        scale 	= new ArrayList<Double>();

        //loop parameters: neuronDiam, scale from GD

        double neuronDiamMin = 3;
        double neuronDiamMax = 3;

		IJ.log("calculating profiles... ");
		long t1 = System.currentTimeMillis();

		int totalConfs = 0;
        for (double neuronDiam = neuronDiamMin; neuronDiam<=neuronDiamMax; neuronDiam++) {
            for (double sc=1.5; sc<=2.5; sc+=1.0) {

				totalConfs++;

                neuronD.add(neuronDiam);
                scale.add(sc);

                Profiler.loadParams(neuronDiam, sc);

                int totalProfiles = Profiler.offsets.size();
                Profiler ms_jobs[] = new Profiler[CPU_NR];
                for (int i = 0; i < ms_jobs.length; i++) {
                    ms_jobs[i] = new Profiler(i*totalProfiles/CPU_NR,  (i+1)*totalProfiles/CPU_NR);
                    ms_jobs[i].start();
                }
                for (int i = 0; i < ms_jobs.length; i++) {
                    try {
                        ms_jobs[i].join();
                    } catch (InterruptedException e) {
                        e.printStackTrace();
                    }
                }

                updateList();

            }
        }

		long t2 = System.currentTimeMillis();
		IJ.log("done "+((t2-t1)/1000f)+" sec.");

		IJ.log("calculating H-domes... ");
		t1 = System.currentTimeMillis();

		float H = 5, Hmin = 3; // GD
		hdomes = new ArrayList<ArrayList<float[]>>(totalLocations);
		for (int locIdx= 0; locIdx<totalLocations; locIdx++) {

			ArrayList<float[]> hdomesPerConf = new ArrayList<float[]>(totalConfs);

			for (int cfgIdx=0; cfgIdx<totalConfs; cfgIdx++) {
				/*
				calculate configuration and add it
				 */
				hdomesPerConf.add(Tools.hdome_Circular(profiles.get(locIdx).get(cfgIdx), H, Hmin));
			}

			hdomes.add(hdomesPerConf);

		}

		t2 = System.currentTimeMillis();
		IJ.log("done "+((t2-t1)/1000f)+" sec.");

        // ms, use profiles
        IJ.log("calculating peaks... ");
        t1 = System.currentTimeMillis();

        Analyzer.loadProfiles(profiles);
        Analyzer ms_jobs[] = new Analyzer[CPU_NR];
        for (int i = 0; i < ms_jobs.length; i++) {
            ms_jobs[i] = new Analyzer(i*totalLocations/CPU_NR,  (i+1)*totalLocations/CPU_NR);
            ms_jobs[i].start();
        }
        for (int i = 0; i < ms_jobs.length; i++) {
            try {
                ms_jobs[i].join();
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

            /*
            calculate backgr at his loc. - local median
             */

            int patchR = 2*(int) Math.ceil(neuronDiamMax);

            float[] neigh = new float[(2*patchR+1)*(2*patchR+1)];

            float currBkgr = 0;

            if (atX>patchR && atY>patchR && atX<inimg.getWidth()-patchR && atY<inimg.getHeight()-patchR) {

                // fill the array in
                int idx = 0;
                for (int locX = atX-patchR; locX<=atX+patchR; locX++) {
                    for (int locY = atY-patchR; locY<=atY+patchR; locY++) {
                        neigh[idx] = inimg.getProcessor().getf(locX, locY);
                        idx++;
                    }
                }

                // take the median as a bkg estimate
                currBkgr = (float) Tools.median_Wirth(neigh);
            }

            bacgr.add(currBkgr);

        }

        IJ.log("done, extracting the scores...");

        double D = 10; // GD  element

        // loop once again and extract scores

        for (int locIdx=0; locIdx<totalLocations; locIdx++) {

            double A0 = centrAvg.get(locIdx).get(0) - bacgr.get(locIdx);
            A0 = (A0>0)? 1-Math.exp(-A0/D) : 0;

            double A1 = 0;// per cluster
            double A2 = 0;
            double A3 = 0;

            int cntCfgs = 0;

            for (int cfgIdx=0; cfgIdx<totalConfs; cfgIdx++) {

                if (angles.get(locIdx).get(cfgIdx)!=null) { // if ms gave 3+ values

                    double A1t = peaks.get(locIdx).get(cfgIdx)[0]- bacgr.get(locIdx);
                    A1t = (A1t>0)? 1-Math.exp(-A1t/D) : 0 ;
                    A1 +=  A1t;

                    double A2t = peaks.get(locIdx).get(cfgIdx)[1] - bacgr.get(locIdx);
                    A2t = (A2t>0)? 1-Math.exp(-A2t/D) : 0 ;
                    A2 += A2t;

                    double A3t = peaks.get(locIdx).get(cfgIdx)[2] - bacgr.get(locIdx);
                    A3t = (A3t>0)? 1-Math.exp(-A3t/D) : 0;
                    A3 += A3t;

                    cntCfgs++;

                }

            }

//            if (cntCfgs>0) {
//                A1 /= cntCfgs;
//                A2 /= cntCfgs;
//                A3 /= cntCfgs;
//            }

            double score = 0;

            if (A0>0 && A1>0 && A2>0 && A3>0) {



                score = A0 * Tools.min3(A1, A2, A3);// ( A1 + A2 + A3 );

//                score += Math.log(A0);
//                score += Math.log(A1);
//                score += Math.log(A2);
//                score += Math.log(A3);
//                score = Math.exp(score/4);
            }

            scoreimg.setf(Profiler.locations[locIdx][0], Profiler.locations[locIdx][1], (float) score);

        }

        /*
        detectio of local max
         */

        int h = 10;     // GD
        int maxIter = 250;
        double eps = 0.00001;

        MS2D.loadTemplate(scoreimg, (ByteProcessor) inmask.getProcessor(), h, maxIter, eps);
        int totalProfiles = MS2D.toProcess;
        MS2D ms_jobs1[] = new MS2D[CPU_NR];
        for (int i = 0; i < ms_jobs.length; i++) {
            ms_jobs1[i] = new MS2D(i*totalProfiles/CPU_NR,  (i+1)*totalProfiles/CPU_NR);
            ms_jobs1[i].start();
        }
        for (int i = 0; i < ms_jobs.length; i++) {
            try {
                ms_jobs1[i].join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }

        IJ.log("done!");

        // clustering
        double dist = 1.0;  // GD
        int M = 50;

        double[][] Tnew = MS2D.T;//extractConvPoints(dist, M);

        IJ.log(Tnew.length+" points extracted");

        Overlay ov = new Overlay();
        for (int i1=0; i1<Tnew.length; i1++) {
            PointRoi pt = new PointRoi(Tnew[i1][1]+0.5, Tnew[i1][0]+0.5);
            pt.setStrokeColor(Color.PINK);
            ov.add(pt);
        }

        inimg.show();
        inimg.setOverlay(ov);
        inimg.getCanvas().addMouseListener(this); // prepare for mouse interaction

        vizProfileImage = new ImagePlus();

        // overlay mask
        inmask.show();
        IJ.selectWindow("inimg");
        IJ.run("Add Image...", "image=inmask x="+0+" y="+0+" opacity=20");
        inmask.close();

        // overlay score
        ImagePlus scoreImg = new ImagePlus("score", scoreimg);
        scoreImg.show();
        scoreImg.setOverlay(ov);
//        IJ.selectWindow("inimg");
//        IJ.run("Add Image...", "image=score x="+0+" y="+0+" opacity=70");
//        scoreImg.close();


        IJ.selectWindow("inimg");
        IJ.setTool("hand");
        inimg.getCanvas().zoomIn(0, 0);
        inimg.getCanvas().zoomIn(0, 0);
//        inimg.getCanvas().zoomIn(0, 0);
//        inimg.getCanvas().zoomIn(0, 0);


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
                Overlay ov = new Overlay();
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

                        if (profiles.get(i).get(g)[j1]<profMin)profMin = profiles.get(i).get(g)[j1];

                    }

                    float[] plotHdomeX = new float[hdomes.get(i).get(g).length];
                    float[] plotHdomeY = new float[plotHdomeX.length];

                    for (int i1=0; i1<plotHdomeX.length; i1++)
                        plotHdomeX[i1] = (float) i1 * angularRes.get(i).get(g);

                    for (int i1=0; i1<plotHdomeY.length; i1++)
                        plotHdomeY[i1] = (hdomes.get(i).get(g)[i1]>0)? 0.9f*profMax : 1.1f*profMin ;

//                    for (int i1=1; i1<plotHdomeY.length; i1++)
//                        plotHdomeY[i1] = (hdomeProfile[(int) Math.round(indexesHdome[i1])]>0)?  plotHdomeY[0] : Float.NaN;
                      //          ((float) Tools.interp1Darray(indexesHdome[i1], hdomeProfile)>0)?
//                                        (float) Tools.interp1Darray(indexesHdome[i1], profilesPerLocation.get(i).get(g));

                    p.setColor(Color.GREEN);
                    p.addPoints(plotHdomeX, plotHdomeY, Plot.LINE);
                    p.draw();

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

                double D = 10;

                System.out.println("\n---\nX: "+atX+" , Y: "+atY+"");

                System.out.println("\n*** CENTER  ***");
                double pk = centrAvg.get(i).get(0);
                double bg = bacgr.get(i);
                double df = (pk-bg>0) ? pk-bg : 0 ;

                System.out.println(" peak: "+pk+" ("+bg+") "+IJ.d2s(1-Math.exp(-df/D), 2)+" | ");

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
