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
import java.util.Arrays;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 7/7/13
 * Time: 8:13 PM
 */
public class ProfilerDemo implements PlugInFilter, MouseListener {

	ImagePlus 	inimg;
	String 		inimgPath;
	ImagePlus   inmask;

    FloatProcessor  scoreimg;

    // visualization
    ImagePlus       vizProfileImage;
    ImageStack      vizProfileStack;

    // to store profiles, list of values per location
    ArrayList<ArrayList<float[]>>   profilesPerLocation;
    ArrayList<ArrayList<Integer>>   angResDegPerLocation;
    ArrayList<ArrayList<String>>    profileNamePerLocation;

    ArrayList<Double>               neuronDPerLocation;
    ArrayList<Double>               scalePerLocation;

    ArrayList<ArrayList<float[]>>   anglesDegPerLocation;   // 3 angles
    ArrayList<ArrayList<float[]>>   scoresPerLocation;      // each angle one score

    ArrayList<ArrayList<Double>>   backgrPerLocation;       // necessary for scoring scheme - only for the location - change and use for scoring
    ArrayList<ArrayList<Double>>   centralAvgPerLocation;

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

        CPU_NR = 4;

		Profiler.loadTemplate(inimg.getProcessor(), (ByteProcessor) inmask.getProcessor());

        int totalLocations = Profiler.locations.length;
        profilesPerLocation     = new ArrayList<ArrayList<float[]>>(totalLocations);
        angResDegPerLocation    = new ArrayList<ArrayList<Integer>>(totalLocations);
        profileNamePerLocation  = new ArrayList<ArrayList<String>>(totalLocations);

        neuronDPerLocation = new ArrayList<Double>(); // this is actually per configuration
        scalePerLocation = new ArrayList<Double>();

        //loop parameters
        for (double neuronDiam = 3; neuronDiam<=3; neuronDiam++) {
            for (double scale=1.5; scale<=2.5; scale+=.5) {

                neuronDPerLocation.add(neuronDiam);
                scalePerLocation.add(scale);

                Profiler.loadParams(neuronDiam, scale);
                IJ.log("calculating profiles... neuronDiam="+neuronDiam+", scale="+scale);
                long t1 = System.currentTimeMillis();
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
                long t2 = System.currentTimeMillis();
                IJ.log("done extracting profiles "+((t2-t1)/1000f)+" sec.");

                updateList();

            }
        }

        // ms, use profilesPerLocation
        IJ.log("calculating local peaks... ");
        long t1 = System.currentTimeMillis();

        Analyzer.loadProfiles(profilesPerLocation);
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

        long t2 = System.currentTimeMillis();
        IJ.log("done extracting peaks "+((t2-t1)/1000f)+" sec.");

        IJ.log("finding matches and corresponding scores for each location...");

        int totLocs = profilesPerLocation.size();

        anglesDegPerLocation    = new ArrayList<ArrayList<float[]>>(totLocs);
        scoresPerLocation       = new ArrayList<ArrayList<float[]>>(totLocs);
        backgrPerLocation       = new ArrayList<ArrayList<Double>>(totLocs);
        centralAvgPerLocation   = new ArrayList<ArrayList<Double>>(totLocs);

        for (int locIdx=0; locIdx<totLocs; locIdx++) {

            int totCfgs = profilesPerLocation.get(locIdx).size();

            ArrayList<float[]>  anglesDegPerCfg = new ArrayList<float[]>(totCfgs);
            ArrayList<float[]>  scoresPerCfg = new ArrayList<float[]>(totCfgs);
            ArrayList<Double>   backgrPerCfg = new ArrayList<Double>(totCfgs);
            ArrayList<Double>   centralAvgPerCfg = new ArrayList<Double>(totCfgs);

            boolean isFirst = true;
            float[] refAnglesDeg    = new float[3];

            for (int cfgIdx=0; cfgIdx<totCfgs; cfgIdx++) {

                /*
                 find central bloc score for this location
                  */

                int atX = Profiler.locations[locIdx][0];
                int atY = Profiler.locations[locIdx][1];
                centralAvgPerCfg.add(centralAvg(atX, atY, neuronDPerLocation.get(cfgIdx), (FloatProcessor) inimg.getProcessor()));

                // get indexes form Analyzer for this conf & loc
                if (Analyzer.peakIdx[locIdx][cfgIdx]!=null) {

                    /*
                    converging angles and scores
                     */

                    float[] anglesDeg       = new float[3];
                    float[] scores          = new float[3];

                    if (isFirst) {

                        // is first time, align as it is
                        refAnglesDeg[0] = anglesDeg[0] = Analyzer.peakIdx[locIdx][cfgIdx][0] * angResDegPerLocation.get(locIdx).get(cfgIdx);
                        refAnglesDeg[1] = anglesDeg[1] = Analyzer.peakIdx[locIdx][cfgIdx][1] * angResDegPerLocation.get(locIdx).get(cfgIdx);
                        refAnglesDeg[2] = anglesDeg[2] = Analyzer.peakIdx[locIdx][cfgIdx][2] * angResDegPerLocation.get(locIdx).get(cfgIdx);

                        // this is raw average, needs to be normalized
                        scores[0] = (float) Tools.interp1Darray(Analyzer.peakIdx[locIdx][cfgIdx][0], profilesPerLocation.get(locIdx).get(cfgIdx));
                        scores[1] = (float) Tools.interp1Darray(Analyzer.peakIdx[locIdx][cfgIdx][1], profilesPerLocation.get(locIdx).get(cfgIdx));
                        scores[2] = (float) Tools.interp1Darray(Analyzer.peakIdx[locIdx][cfgIdx][2], profilesPerLocation.get(locIdx).get(cfgIdx));

                        isFirst = false;

                    }
                    else {
                        // not the first time, check angular matches
                        // read them and store wrt first
                        anglesDeg[0] = Analyzer.peakIdx[locIdx][cfgIdx][0] * angResDegPerLocation.get(locIdx).get(cfgIdx);
                        anglesDeg[1] = Analyzer.peakIdx[locIdx][cfgIdx][1] * angResDegPerLocation.get(locIdx).get(cfgIdx);
                        anglesDeg[2] = Analyzer.peakIdx[locIdx][cfgIdx][2] * angResDegPerLocation.get(locIdx).get(cfgIdx);


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
                        float[] scoresTemp      = new float[3];

                        for (int t=0; t<mapping.length; t++) {
                            anglesDegTemp[mapping[t]] = anglesDeg[t];
                            scoresTemp[mapping[t]] = (float) Tools.interp1Darray(Analyzer.peakIdx[locIdx][cfgIdx][t], profilesPerLocation.get(locIdx).get(cfgIdx));
                        }

//                        IJ.log("...");
//                        IJ.log(Arrays.toString(anglesDegTemp));
//                        IJ.log("...");

                        for (int t=0; t<mapping.length; t++) {
                            anglesDeg[t] = anglesDegTemp[t];
                            scores[t] = scoresTemp[t];
                        }

                    }

//                    IJ.log("***");
//                    IJ.log(Arrays.toString(anglesDeg));
//                    IJ.log(Arrays.toString(refAnglesDeg));
//                    IJ.log("---");

                    anglesDegPerCfg.add(anglesDeg);
                    scoresPerCfg.add(scores);

                    /*
                    estimate background (median)
                     */
                    backgrPerCfg.add(Tools.median_Wirth(profilesPerLocation.get(locIdx).get(cfgIdx).clone()));

                }
                else {

                    anglesDegPerCfg.add(null);
                    scoresPerCfg.add(null);
                    backgrPerCfg.add(Double.NaN);

                }

            }

            anglesDegPerLocation.add(anglesDegPerCfg);
            scoresPerLocation.add(scoresPerCfg);
            centralAvgPerLocation.add(centralAvgPerCfg);
            backgrPerLocation.add(backgrPerCfg);

        }

        IJ.log("done!");

        // loop once again and extract scores
        for (int locIdx=0; locIdx<totLocs; locIdx++) {

            int totCfgs = profilesPerLocation.get(locIdx).size();

            double A0 = 0;
            double A1 = 0;
            double A2 = 0;
            double A3 = 0;

            int cntCfgs = 0;

            for (int cfgIdx=0; cfgIdx<totCfgs; cfgIdx++) {

                if (anglesDegPerLocation.get(locIdx).get(cfgIdx)!=null) {

                    A0 += centralAvgPerLocation.get(locIdx).get(cfgIdx);
                    A1 += scoresPerLocation.get(locIdx).get(cfgIdx)[0];// - backgrPerLocation.get(locIdx).get(cfgIdx);
                    A2 += scoresPerLocation.get(locIdx).get(cfgIdx)[1];// - backgrPerLocation.get(locIdx).get(cfgIdx);
                    A3 += scoresPerLocation.get(locIdx).get(cfgIdx)[2];// - backgrPerLocation.get(locIdx).get(cfgIdx);

                    cntCfgs++;

                }

            }

            if (cntCfgs>0) {
                A0 /= cntCfgs;
                A1 /= cntCfgs;
                A2 /= cntCfgs;
                A3 /= cntCfgs;
            }

            // geometric mean for cluster averages   sqrt4th(A0*A1*A2*A3)
            double sum = 0;

            if (A0>0 && A1>0 && A2>0 && A3>0) {
                sum += Math.log(A0);
                sum += Math.log(A1);
                sum += Math.log(A2);
                sum += Math.log(A3);
                sum = Math.exp(sum/4);
            }

            scoreimg.setf(Profiler.locations[locIdx][0], Profiler.locations[locIdx][1], (float) sum);

        }
        new ImagePlus("score", scoreimg).show();

        // prepare for clicking

        inimg.show();
        inimg.getCanvas().addMouseListener(this);

        vizProfileImage = new ImagePlus();

        inmask.show();
        IJ.selectWindow("inimg");
        IJ.run("Add Image...", "image=inmask x="+0+" y="+0+" opacity=30");
        inmask.close();

        IJ.selectWindow("inimg");
        IJ.setTool("hand");
        inimg.getCanvas().zoomIn(0, 0);
        inimg.getCanvas().zoomIn(0, 0);
        inimg.getCanvas().zoomIn(0, 0);
        inimg.getCanvas().zoomIn(0, 0);


	}

    private void updateList() // uses Profiler
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

            if (profilesPerLocation.size()<Profiler.locations.length) {

                ArrayList<float[]> A = new ArrayList<float[]>();
                A.add(profileToAdd);
                profilesPerLocation.add(A);

                ArrayList<String> B = new ArrayList<String>();
                B.add(stringToAdd);
                profileNamePerLocation.add(B);

                ArrayList<Integer> C = new ArrayList<Integer>();
                C.add(intToAdd);
                angResDegPerLocation.add(C);

            }
            else {
                profilesPerLocation.get(loopLocations).add(profileToAdd);
                profileNamePerLocation.get(loopLocations).add(stringToAdd);
                angResDegPerLocation.get(loopLocations).add(intToAdd);
            }

        }
    }

    public double centralAvg(int atX, int atY, double diam, FloatProcessor inip)
    {
        double out = 0;
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

                int nrProfiles = profilesPerLocation.get(i).size();

                /*
				overlay for the main image
			    */
                Overlay ov = new Overlay();
                PointRoi pt = new PointRoi(atX+.5, atY+.5);
                ov.add(pt);

                for (int g=0; g<nrProfiles; g++) {  // g will loop profiles

                    int currProfileLen = profilesPerLocation.get(i).get(g).length;

                    // allocate profile mean
                    float[] profileMean = new float[currProfileLen];
                    for (int i1=0; i1<currProfileLen; i1++) profileMean[i1] = centralAvgPerLocation.get(i).get(g).floatValue();

                    // allocate profile median
                    float[] profileMed = new float[currProfileLen];
                    for (int i1=0; i1<currProfileLen; i1++) profileMed[i1] = backgrPerLocation.get(i).get(g).floatValue();

                    /*
                    plot ij  template
                     */
                    float[] angIdx = new float[currProfileLen];

                    for (int q=0; q<currProfileLen; q++) {
                        angIdx[q] = q * angResDegPerLocation.get(i).get(g);
                    }

                    /*
                    plot def.
                     */

                    // profile
                    Plot p = new Plot("", "orient.[deg]", profileNamePerLocation.get(i).get(g), angIdx, profilesPerLocation.get(i).get(g));
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
                    float[] plotConvX = new float[100];
                    float[] plotConvY = new float[100];

                    for (int i1=0; i1<100; i1++)
                        plotConvX[i1] = Analyzer.convIdx.get(i).get(g)[i1] * angResDegPerLocation.get(i).get(g);

                    for (int i1=0; i1<100; i1++)
                        plotConvY[i1] = (float) Tools.interp1Darray(Analyzer.convIdx.get(i).get(g)[i1], profilesPerLocation.get(i).get(g));

                    p.addPoints(plotConvX, plotConvY, Plot.X);

                    // add h-dome plot
                    float[] hdomeProfile = Tools.hdome_Circular(profilesPerLocation.get(i).get(g), 5, 3);

                    float[] indexesHdome = new float[hdomeProfile.length];
                    float[] plotHdomeX = new float[hdomeProfile.length];
                    float[] plotHdomeY = new float[hdomeProfile.length];

                    for (int i1=0; i1<indexesHdome.length; i1++)
                        indexesHdome[i1] = (float) i1;// / indexesHdome.length;

                    for (int i1=0; i1<plotHdomeX.length; i1++)
                        plotHdomeX[i1] = indexesHdome[i1] * angResDegPerLocation.get(i).get(g);

                    plotHdomeY[0] = Float.MIN_VALUE;
                    for (int j1=0; j1<profilesPerLocation.get(i).get(g).length; j1++) {
                        if (profilesPerLocation.get(i).get(g)[j1]>plotHdomeY[0])
                            plotHdomeY[0] = profilesPerLocation.get(i).get(g)[j1];
                    }

                    plotHdomeY[0] *= 0.9;

                    for (int i1=1; i1<plotHdomeY.length; i1++)
                        plotHdomeY[i1] = (hdomeProfile[(int) Math.round(indexesHdome[i1])]>0)?  plotHdomeY[0] : Float.NaN;
                      //          ((float) Tools.interp1Darray(indexesHdome[i1], hdomeProfile)>0)?
//                                        (float) Tools.interp1Darray(indexesHdome[i1], profilesPerLocation.get(i).get(g));
                    //: Float.NaN ;

                    p.setColor(Color.RED);
                    p.addPoints(plotHdomeX, plotHdomeY, Plot.BOTTOM_MARGIN);
                    p.draw();

                    p.setColor(Color.BLUE);
                    p.addPoints(angIdx, profileMean, Plot.LINE);
                    p.setColor(Color.GREEN);
                    p.addPoints(angIdx, profileMed, Plot.LINE);

                    // plot dtections if there were any
					if (Analyzer.peakIdx[i][g]!=null) {


                         /*
                         i - location, g - configuration
                          */


                        if (anglesDegPerLocation.get(i).get(g)!=null) {

                            float pX, pY;
                            pX = anglesDegPerLocation.get(i).get(g)[0];
                            pY = scoresPerLocation.get(i).get(g)[0];
                            p.setColor(Color.RED);
                            p.addPoints(new float[]{pX}, new float[]{pY}, Plot.BOX);
                            p.draw();

                            pX = anglesDegPerLocation.get(i).get(g)[1];
                            pY = scoresPerLocation.get(i).get(g)[1];
                            p.setColor(Color.GREEN);
                            p.addPoints(new float[]{pX}, new float[]{pY}, Plot.BOX);
                            p.draw();

                            pX = anglesDegPerLocation.get(i).get(g)[2];
                            pY = scoresPerLocation.get(i).get(g)[2];
                            p.setColor(Color.BLUE);
                            p.addPoints(new float[]{pX}, new float[]{pY}, Plot.BOX);
                            p.draw();

                        }

                        // modify xAxis angles into euclidean points
                        double[][] msPeaks = new double[3][2];
                        double rd = neuronDPerLocation.get(g)*scalePerLocation.get(g);

                        float[] xAxis = new float[3];
                        for (int i1=0; i1<3; i1++) xAxis[i1] = Analyzer.peakIdx[i][g][i1] * angResDegPerLocation.get(i).get(g); // this is precalculated now

                        float pointX, pointY;
                        PointRoi point;
                        // cluster 0 - RED
                        // i - location, g-configuration, 0-cluster
                        pointX = (float) (atX+rd*Math.cos(anglesDegPerLocation.get(i).get(g)[0]*(2*Math.PI/360))+.5);
                        pointY = (float) (atY-rd*Math.sin(anglesDegPerLocation.get(i).get(g)[0]*(2*Math.PI/360))+.5);
                        point = new PointRoi(pointX, pointY);
                        point.setStrokeColor(Color.RED);
                        ov.add(point);

                        // cluster 1 - GREEN
                        // i - location, g-configuration
                        pointX = (float) (atX+rd*Math.cos(anglesDegPerLocation.get(i).get(g)[1]*(2*Math.PI/360))+.5);
                        pointY = (float) (atY-rd*Math.sin(anglesDegPerLocation.get(i).get(g)[1]*(2*Math.PI/360))+.5);
                        point = new PointRoi(pointX, pointY);
                        point.setStrokeColor(Color.GREEN);
                        ov.add(point);

                        // cluster 2 - BLUE
                        // i - location, g-configuration,
                        pointX = (float) (atX+rd*Math.cos(anglesDegPerLocation.get(i).get(g)[2]*(2*Math.PI/360))+.5);
                        pointY = (float) (atY-rd*Math.sin(anglesDegPerLocation.get(i).get(g)[2]*(2*Math.PI/360))+.5);
                        point = new PointRoi(pointX, pointY);
                        point.setStrokeColor(Color.BLUE);
                        ov.add(point);

                        ov.add(new OvalRoi(atX-rd+.5, atY-rd+.5, 2*rd, 2*rd));

                    }

                    inimg.setOverlay(ov);

                    //vizProfileStack.addSlice("", p.getProcessor());

                    vizProfileStack.addSlice("", p.getProcessor());

                    /*
                    export R
                     */

                    String printProfile = "";
                    String printAngle = "";
                    String xName = "ang_"+profileNamePerLocation.get(i).get(g);
                    String yName = profileNamePerLocation.get(i).get(g);

                    printProfile    += yName+" <- c(";
                    printAngle      += xName+" <- c(";

                    for (int i1=0; i1<currProfileLen; i1++) {

                        printProfile+=profilesPerLocation.get(i).get(g)[i1]+"";
                        printAngle+=angIdx[i1]+"";

                        if(i1<profilesPerLocation.get(i).get(g).length-1) {
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
