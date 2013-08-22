package profile;

import conn.Find_Connected_Regions;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.*;
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
	String			logFile;
	PrintWriter 	logWriter;

    float scale;
	float scale_A;
	float scale_B;               // only one scale now
    double neuronDiamMax, D;
    ArrayList<ArrayList<float[]>>   profiles;

    PlotWindow pwP_A, pwI_A, pwP_B, pwI_B;

    Overlay     currOvl = new Overlay();
    Overlay 	detectionOverlay;

    // DEBUG
    ByteProcessor ringAip;
    ByteProcessor ringBip;

    // some variables (used by file export (buttonPressed) and mouse click plots (mouseCLicked))
    float[] ang_A, ang_B;  // they cover [0, 360) but with different sampling rate
    float[] i_A, i_B;
    float[] iTh_A, iTh_B;
    float[] peakAng_A, peakAng_B;
    float[] peakIdx_A, peakIdx_B;
    float[] profile_A = null, profile_B = null;

    int CPU_NR;

    private static float Deg2Rad = (float) (Math.PI/180f);
    private static float Rad2Deg = (float) (180f/Math.PI);
    private static float MIN_COS_ANG = 0.85f;

	public int setup(String s, ImagePlus imagePlus) {
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

		int totalJobs; // used in paralelization lines
		long t1, t2;

		neuronDiamMax       =  Prefs.get("advantra.critpoint.neuronDiamMax", 3);
		scale               = (float) Prefs.get("advantra.critpoint.scale", 1.5);

		GenericDialog gd = new GenericDialog("JUNCTION DET.");
		gd.addNumericField("neuronD ",  neuronDiamMax, 0, 10, " MAX");
		gd.addNumericField("scale ",  scale, 1, 10, " x(neuronD)");

		gd.showDialog();
		if (gd.wasCanceled()) return;

		neuronDiamMax       =  gd.getNextNumber();
		Prefs.set("advantra.critpoint.neuronDiamMax",   neuronDiamMax);

		scale               = (float) gd.getNextNumber();
		Prefs.set("advantra.critpoint.scale",   scale);

		CPU_NR = Runtime.getRuntime().availableProcessors();

		logFile = "profile_v11_output.log";
		// empty the log file
		try {
			PrintWriter writer = new PrintWriter(logFile);
			writer.print("");
			writer.close();
		}
		catch (FileNotFoundException ex) {}
//		try {
//			//out = new PrintWriter(logFile);
//		}
//		catch (FileNotFoundException ex) {}

        /*
        ***********************************************************
         */

		int neighbourhoodR = (int) Math.ceil(scale*neuronDiamMax);
		IJ.log("extracting background... r = "+neighbourhoodR);

		t1 = System.currentTimeMillis();

		Masker.loadTemplate(imp.getProcessor(), neighbourhoodR);

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

//		ImagePlus inmask = new ImagePlus("inmask", Masker.mask);
//		inmask.show();

		int nr = 0;
		for (int i=0; i<totalJobs; i++) if (Masker.mask.getf(i)==255) nr++;

		IJ.log("done "+((t2-t1)/1000f)+" sec. total "+nr+" check locations given by mask ("+(100*(float)nr/(totalJobs))+" % kept)");

		/*
        ***********************************************************
         */

		IJ.log("calculating profiles... ");
		t1 = System.currentTimeMillis();

		// set the template first
		Profiler.loadTemplate(imp.getProcessor(), Masker.mask);

		// radiuses and scales for rings
		scale_A = scale;
		scale_B = scale_A + 1;

		Profiler.loadParams(neuronDiamMax, scale_A, false);
		totalJobs = Profiler.offsets.size();

		Profiler[] profiler_jobs;
		profiles     	= new ArrayList<ArrayList<float[]>>(Profiler.locations.length);

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

		addProfilerList(profiles, Profiler.profiles);

		Profiler.loadParams(neuronDiamMax, scale_B, false);
		totalJobs = Profiler.offsets.size();
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

		addProfilerList(profiles, Profiler.profiles);

		t2 = System.currentTimeMillis();
		IJ.log("done "+((t2-t1)/1000f)+" sec.");

        /*
        ***********************************************************
         */

		IJ.log("calculating peaks for selected locations...");
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

		// loop once again and extract scores
		IJ.log("logging scores...");

		// initialize file
		try {
			logWriter = new PrintWriter(new BufferedWriter(new FileWriter(logFile, true)));
			DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
			Date date = new Date();
			logWriter.println("BIFURCATION DETECTION LOG\n"+dateFormat.format(date));
		} catch (IOException e) {}

		ByteProcessor scoreimg = new ByteProcessor(imp.getWidth(), imp.getHeight());
        detectionOverlay = new Overlay();
		for (int locIdx=0; locIdx<Profiler.locations.length; locIdx++) {
			int sc =
					score(
                    	Profiler.locations[locIdx][0], Profiler.locations[locIdx][1],
                    	(FloatProcessor) imp.getProcessor(),
						Masker.backgr,
						Masker.margin,
                    	Analyzer.peakIdx.get(locIdx),
						scale_A, scale_B);


			scoreimg.set(Profiler.locations[locIdx][0], Profiler.locations[locIdx][1], sc);
            if (sc>0) {
                detectionOverlay.add(new PointRoi(Profiler.locations[locIdx][0]+.5, Profiler.locations[locIdx][1]+.5));
            }

		}

        /*
        -----------------------------------------------
         */

		logWriter.close();
		IJ.log("log saved to "+new File(logFile).getAbsolutePath());

		ImagePlus scoreImagePlus = new ImagePlus("score", scoreimg);
		//scoreImagePlus.show();
		IJ.selectWindow(inimgTitle);
        IJ.setTool("hand");

        /*
        -----------------------------------------------
         */

		t1 = System.currentTimeMillis();
		Find_Connected_Regions conn_reg = new Find_Connected_Regions(scoreImagePlus, true);
		conn_reg.run("");
		t2 = System.currentTimeMillis();
		int nr_regions = conn_reg.getNrConnectedRegions();
		IJ.log(nr_regions+" connected regions extracted.\n" + "elapsed: "+((t2-t1)/1000f)+ " seconds.");

		//ImagePlus imageLabels = conn_reg.showLabels();
		//imageLabels.show();

        /*
        -----------------------------------------------
         */

        ringAip = new ByteProcessor(imp.getWidth(), imp.getHeight());
        ringBip = new ByteProcessor(imp.getWidth(), imp.getHeight());
        for (int ii=0; ii<imp.getWidth()*imp.getHeight(); ii++) ringAip.set(ii,99);

        for (int ii=0; ii<Profiler.locations.length; ii++) {
            int x = Profiler.locations[ii][0];
            int y = Profiler.locations[ii][1];
            ringAip.set(x,y,Analyzer.peakIdx.get(ii).get(0).size());
            ringBip.set(x,y,Analyzer.peakIdx.get(ii).get(1).size());
        }

        //detectionOverlay = formOverlay(conn_reg.getConnectedRegions());
		canvas.setOverlay(detectionOverlay);
		canvas.addMouseListener(this);
		canvas.addMouseMotionListener(this);
		canvas.addKeyListener(this);

		ImagePlus outIm1 = new ImagePlus("viz1", imp.getProcessor());
		outIm1.setOverlay(formPointOverlay(conn_reg.getConnectedRegions(), 0));
		outIm1.show();


//        ImagePlus outIm = new ImagePlus("viz", imp.getProcessor());
//        outIm.setOverlay(formOverlay(conn_reg.getConnectedRegions()));
//        outIm.show();

        IJ.selectWindow(inimgTitle);
        IJ.setTool("hand");

	}

    public void mouseClicked(MouseEvent e)
    {

		// on mouse click extract locations, plot profile, intensities along profile and overlay points
        currOvl.clear();
        currOvl = detectionOverlay.duplicate();

        int offscreenX = canvas.offScreenX(e.getX());
        int offscreenY = canvas.offScreenY(e.getY());

        IJ.log(""+offscreenX+", "+offscreenY+" \n");

        //float cI = Interpolator.interpolateAt(offscreenX, offscreenY, (FloatProcessor) imp.getProcessor());
        //float cB = Interpolator.interpolateAt(offscreenX, offscreenY, (FloatProcessor) Masker.backgr);
        PointRoi pt = new PointRoi(offscreenX+.5, offscreenY+.5);
        currOvl.add(pt);

        // define ring A
        double rd_A = neuronDiamMax*scale_A;
        OvalRoi ring_A = new OvalRoi(offscreenX-rd_A+.5, offscreenY-rd_A+.5, 2*rd_A, 2*rd_A);
        currOvl.add(ring_A);

        profile_A = Profiler.extractProfile(neuronDiamMax, scale_A, offscreenX, offscreenY, (FloatProcessor) imp.getProcessor());
        peakIdx_A = Analyzer.extractPeakIdxs(profile_A); // MS returns values range [0, length)

        ArrayList<Float> A_Ang  = new ArrayList<Float>(4); // store those that were selected as ok
        ArrayList<Float> A_x    = new ArrayList<Float>(4);
        ArrayList<Float> A_y    = new ArrayList<Float>(4);

		peakAng_A   = null;
        if (peakIdx_A!=null) {
            float[] peakX_A, peakY_A, peakI_A, peakB_A;
            peakAng_A = new float[peakIdx_A.length];
                    peakX_A = new float[peakIdx_A.length];
                            peakY_A = new float[peakIdx_A.length];
                                    peakI_A = new float[peakIdx_A.length];
                                            peakB_A = new float[peakIdx_A.length];

            for (int i=0; i<peakIdx_A.length; i++) {
                peakAng_A[i] = peakIdx_A[i] * Profiler.getResolDeg(scale_A) * Deg2Rad;
                peakX_A[i] = (float) (offscreenX + rd_A * Math.cos( peakAng_A[i] ));
                peakY_A[i] = (float) (offscreenY - rd_A * Math.sin( peakAng_A[i] ));
                peakI_A[i] = Interpolator.interpolateAt(peakX_A[i], peakY_A[i], (FloatProcessor) imp.getProcessor());
                peakB_A[i] = Interpolator.interpolateAt(peakX_A[i], peakY_A[i], Masker.backgr);

//					PointRoi pt = new PointRoi(peakX_A[i]+.5, peakY_A[i]+.5);
//                  pt.setStrokeColor(getColor(i));
//                  currOvl.add(pt);

				A_Ang.add(peakIdx_A[i] * Profiler.getResolDeg(scale_A)); // because hungarian matching method will take angle differences in degrees
                A_x.add(peakX_A[i]);
                A_y.add(peakY_A[i]);

            }
        }

        /*
        visualizations
         */

		Plot chartP_A, chartI_A;
        float[] profile_A_MinMax = Tools.getMinMax(profile_A);

        if (pwI_A == null) {
            ang_A   = new float[profile_A.length];
            i_A     = new float[profile_A.length];
            iTh_A   = new float[profile_A.length];
        }

        for (int i = 1; i <= ang_A.length; i++) {
            ang_A[i - 1] = (i-1)*Profiler.getResolDeg(scale_A);
            float pX = (float) (offscreenX + rd_A * Math.cos( ang_A[i - 1] * Deg2Rad ));
            float pY = (float) (offscreenY - rd_A * Math.sin( ang_A[i - 1] * Deg2Rad ));
            i_A[i-1]    = Interpolator.interpolateAt(pX, pY, (FloatProcessor) imp.getProcessor());
            iTh_A[i-1]  = Interpolator.interpolateAt(pX, pY, Masker.backgr) + Interpolator.interpolateAt(pX, pY, Masker.margin);
        }

        chartP_A = new Plot("", "", "", ang_A, profile_A);
        chartP_A.setSize(500, 250);
        chartP_A.addPoints(ang_A, profile_A, PlotWindow.CIRCLE);
        chartP_A.draw();  chartP_A.setColor(Color.RED);
        // draw peak detections
        if (peakAng_A!=null)
            for (int i = 0; i < peakIdx_A.length; i++) {
                chartP_A.drawLine(peakAng_A[i] * Rad2Deg, profile_A_MinMax[0], peakAng_A[i] * Rad2Deg, Tools.interp1Darray(peakIdx_A[i], profile_A));
            }

        chartI_A = new Plot("", "", "", ang_A, i_A);
        chartI_A.setSize(500, 250);
		chartI_A.addPoints(ang_A, i_A, PlotWindow.CIRCLE);
        chartI_A.draw(); chartI_A.setColor(Color.BLUE);
        chartI_A.addPoints(ang_A, iTh_A, PlotWindow.LINE);
        chartI_A.draw();  chartI_A.setColor(Color.RED);
        // draw peak detections
        if (peakAng_A!=null)
            for (int i = 0; i < peakIdx_A.length; i++)
                chartI_A.drawLine(peakAng_A[i] * Rad2Deg, Tools.interp1Darray(peakIdx_A[i], iTh_A), peakAng_A[i] * Rad2Deg, Tools.interp1Darray(peakIdx_A[i], i_A));

        if (pwP_A == null) pwP_A = chartP_A.show();
        pwP_A.drawPlot(chartP_A);
        pwP_A.setTitle("RING A: profile, x = " + offscreenX + ", y = " + offscreenY);

        if (pwI_A == null) pwI_A = chartI_A.show();
        pwI_A.drawPlot(chartI_A);
        pwI_A.setTitle("RING A: intenst, x = " + offscreenX + ", y = " + offscreenY);

        /*
        **************************
         */

        // define ring B, similarly as A
//        scale_B = scale_A + 1;
        double rd_B = neuronDiamMax*scale_B;
        OvalRoi ringB = new OvalRoi(offscreenX-rd_B+.5, offscreenY-rd_B+.5, 2*rd_B, 2*rd_B);
        currOvl.add(ringB);

        profile_B = Profiler.extractProfile(neuronDiamMax, scale_B, offscreenX, offscreenY, (FloatProcessor) imp.getProcessor());
        peakIdx_B = Analyzer.extractPeakIdxs(profile_B); // MS values range [0, length)

        ArrayList<Float> B_Ang  = new ArrayList<Float>(4); // store those that were selected as ok
        ArrayList<Float> B_x    = new ArrayList<Float>(4);
        ArrayList<Float> B_y    = new ArrayList<Float>(4);

        peakAng_B   = null;
        if (peakIdx_B!=null) {
            float[] peakX_B, peakY_B, peakI_B, peakB_B;
            peakAng_B = new float[peakIdx_B.length];
            peakX_B = new float[peakIdx_B.length];
            peakY_B = new float[peakIdx_B.length];
            peakI_B = new float[peakIdx_B.length];
            peakB_B = new float[peakIdx_B.length];

            for (int i=0; i<peakIdx_B.length; i++) {
                peakAng_B[i] = peakIdx_B[i] * Profiler.getResolDeg(scale_B) * Deg2Rad;
                peakX_B[i] = (float) (offscreenX + rd_B * Math.cos( peakAng_B[i] ));
                peakY_B[i] = (float) (offscreenY - rd_B * Math.sin( peakAng_B[i] ));
                peakI_B[i] = Interpolator.interpolateAt(peakX_B[i], peakY_B[i], (FloatProcessor) imp.getProcessor());
                peakB_B[i] = Interpolator.interpolateAt(peakX_B[i], peakY_B[i], Masker.backgr);

//				PointRoi pt = new PointRoi(peakX_B[i]+.5, peakY_B[i]+.5);
//              pt.setStrokeColor(getColor(i));
//              currOvl.add(pt);

				B_Ang.add(peakIdx_B[i] * Profiler.getResolDeg(scale_B)); // because hungarian matching method will take angle differences in degrees
                B_x.add(peakX_B[i]);
                B_y.add(peakY_B[i]);

            }
        }

        /*
        visualizations
         */

        Plot chartP_B, chartI_B;
        float[] profile_B_MinMax = Tools.getMinMax(profile_B);

        if (pwI_B == null) {
            ang_B   = new float[profile_B.length];
            i_B     = new float[profile_B.length];
            iTh_B   = new float[profile_B.length];
        }

        for (int i = 1; i <= ang_B.length; i++) {
            ang_B[i - 1] = (i-1)*Profiler.getResolDeg(scale_B);
            float pX = (float) (offscreenX + rd_B * Math.cos( ang_B[i - 1] * Deg2Rad ));
            float pY = (float) (offscreenY - rd_B * Math.sin( ang_B[i - 1] * Deg2Rad ));
            i_B[i-1]  = Interpolator.interpolateAt(pX, pY, (FloatProcessor) imp.getProcessor());
            iTh_B[i-1]  = Interpolator.interpolateAt(pX, pY, Masker.backgr) + Interpolator.interpolateAt(pX, pY, Masker.margin);
        }

        chartP_B = new Plot("", "", "", ang_B, profile_B);
        chartP_B.setSize(500, 250);
        chartP_B.addPoints(ang_B, profile_B, PlotWindow.CIRCLE);
        chartP_B.draw();  chartP_B.setColor(Color.RED);
        // draw peak detections
        if (peakAng_B!=null)
            for (int i = 0; i < peakIdx_B.length; i++)
                chartP_B.drawLine(peakAng_B[i] * Rad2Deg, profile_B_MinMax[0], peakAng_B[i] * Rad2Deg, Tools.interp1Darray(peakIdx_B[i], profile_B));

        chartI_B = new Plot("", "", "", ang_B, i_B);
        chartI_B.setSize(500, 250);
		chartI_B.addPoints(ang_B, i_B, PlotWindow.CIRCLE);
        chartI_B.draw(); chartI_B.setColor(Color.BLUE);
        chartI_B.addPoints(ang_B, iTh_B, PlotWindow.LINE);
        chartI_B.draw();  chartI_B.setColor(Color.RED);
        // draw peak detections
        if (peakAng_B!=null)
            for (int i = 0; i < peakIdx_B.length; i++)
                chartI_B.drawLine(peakAng_B[i] * Rad2Deg, Tools.interp1Darray(peakIdx_B[i], iTh_B), peakAng_B[i] * Rad2Deg, Tools.interp1Darray(peakIdx_B[i], i_B));

        if (pwP_B == null) pwP_B = chartP_B.show();
        pwP_B.drawPlot(chartP_B);
        pwP_B.setTitle("RING B: profile, x = " + offscreenX + ", y = " + offscreenY);

        if (pwI_B == null) pwI_B = chartI_B.show();
        pwI_B.drawPlot(chartI_B);
        pwI_B.setTitle("RING B: intenst, x = " + offscreenX + ", y = " + offscreenY);

        /*
        **************************
         */

        int[][] map = Tools.hungarianMappingAnglesDeg(A_Ang, B_Ang); // map(A ring index, B ring index)

        // after mapping
        float[] angDiv = new float[map.length];

        for (int k=0; k<map.length; k++) {

            // point in A : A_x.get(map[k][0], A_y.get(map[k][0])
            float ax = A_x.get(map[k][0]);
            float ay = A_y.get(map[k][0]);
			// add it
            pt = new PointRoi(ax+.5, ay+.5);
            pt.setStrokeColor(getColor(k));
            currOvl.add(pt);

            // point in B : B_x.get(map[k][1]), B_y.get(map[k][1])
            float bx = B_x.get(map[k][1]);
            float by = B_y.get(map[k][1]);
			// add it
            pt = new PointRoi(bx+.5, by+.5);
            pt.setStrokeColor(getColor(k));
            currOvl.add(pt);

			IJ.log("mapping # "+offscreenX+" , "+offscreenY+" , "+ax+" , "+ay+" , "+bx+" , "+by+"\n");

            angDiv[k] = Tools.angularDeviation(offscreenX, offscreenY, ax, ay, bx, by);

            IJ.log("mapping # "+k+" cos(ang) = "+angDiv[k]+" : "+A_Ang.get(map[k][0])+" \t "+B_Ang.get(map[k][1]));

			if (angDiv[k]>MIN_COS_ANG) {

				PolygonRoi pr = new PolygonRoi(
													  new float[]{offscreenX+.5f, ax+.5f, bx+.5f},
													  new float[]{offscreenY+.5f, ay+.5f, by+.5f},
													  3,
													  PolygonRoi.FREELINE
													  );
				currOvl.add(pr);

			}

        }

		//System.out.println("\nCLUSTER ANG DIVS "+Arrays.toString(angDiv));
		//IJ.log("\nCLUSTER ANG DIVS "+Arrays.toString(angDiv));

        /*
        **************************
         */

        canvas.setOverlay(currOvl);


    }

    private Overlay formOverlay(ArrayList<ArrayList<int[]>> regs)
    {

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

                if (!(B>Tools.VERY_SMALL_POSITIVE)) {

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
                ovroi.setStrokeWidth(3);
                detections.add(ovroi);

            }
        }

        return detections;

    }

    private int score(
							 float locX, float locY,
							 FloatProcessor inip1,
							 FloatProcessor backgr1,
							 FloatProcessor margin1,
							 ArrayList<ArrayList<Float>> peakIdxsAtLoc, //[2][nr. conv. points]
							 float scale_A, float scale_B) 				// 2 rings specific
	{
		logWriter.println(IJ.d2s(locX, 0) + "," + IJ.d2s(locY, 0) + " -> ");

		/*******************/

		// precheck
		int nrPointsA = peakIdxsAtLoc.get(0).size();
		if (nrPointsA < 3) {
			logWriter.println("ring A : "+peakIdxsAtLoc.get(0).size()+" peaks => NO");
			return 0;
		}

		int nrPointsB = peakIdxsAtLoc.get(1).size();
		if (nrPointsB < 3) {
			logWriter.println("ring B : "+peakIdxsAtLoc.get(1).size()+" peaks => NO");
			return 0;
		}

		// before mapping (need to have them as lists to match)
		ArrayList<ArrayList<Float>> peakAngsAtLoc = new ArrayList<ArrayList<Float>>(2);

		ArrayList<Float> peakAngA 	= new ArrayList<Float>(nrPointsA);
		for (int ii=0; ii<nrPointsA; ii++)
			peakAngA.add(ii, peakIdxsAtLoc.get(0).get(ii) * Profiler.getResolDeg(scale_A));

		peakAngsAtLoc.add(peakAngA);

		ArrayList<Float> peakAngB   = new ArrayList<Float>(nrPointsB);
		for (int ii=0; ii<nrPointsB; ii++)
			peakAngB.add(ii, peakIdxsAtLoc.get(1).get(ii) * Profiler.getResolDeg(scale_B));

		peakAngsAtLoc.add(peakAngB);

//		log += "ring angles: ";
//		for (int i=0; i<peakAngsAtLoc.size(); i++) {
//			log += "ring "+i+" : ";
//			for (int j=0; j<peakAngsAtLoc.get(i).size(); j++) {
//				log += peakAngsAtLoc.get(i).get(j)+" ";
//			}
//		}
//		log += "\n";

		// map[k][0] is index in peakAngsAtLoc.get(0), map[k][1] is index in peakAngsAtLoc.get(1)
		int[][] map = Tools.hungarianMappingAnglesDeg(peakAngsAtLoc.get(0), peakAngsAtLoc.get(1));

		logWriter.println(map.length + " matches"); // map.length >= 3

		// remap points in arrays
		float[][][] p = new float[2][map.length][2]; // [ring idx][cluster idx][xy]

		for (int k=0; k<map.length; k++) { // k loops clusters
				// ring A   take map[k][0] element
				p[0][k][0] =
				(float) (locX + neuronDiamMax * scale_A * Math.cos(peakAngsAtLoc.get(0).get(map[k][0]) * Deg2Rad));
				p[0][k][1] =
				(float) (locY - neuronDiamMax * scale_A * Math.sin(peakAngsAtLoc.get(0).get(map[k][0]) * Deg2Rad));
				// ring B   take map[k][1] element
				p[1][k][0] =
				(float) (locX + neuronDiamMax * scale_B * Math.cos(peakAngsAtLoc.get(1).get(map[k][1]) * Deg2Rad));
				p[1][k][1] =
				(float) (locY - neuronDiamMax * scale_B * Math.sin(peakAngsAtLoc.get(1).get(map[k][1]) * Deg2Rad));
		}

		// check each cluster of the 2 ring combination for directionality ( mapping does not mean they're aligned)
		ArrayList<Integer> idxDirectedClusters = new ArrayList<Integer>(map.length);
		for (int k=0; k<map.length; k++) {  // 3 or 4 clusters

			logWriter.println("mapping # "+locX+" , "+locY+" , "+p[0][k][0]+" , "+p[0][k][1]+" , "+p[1][k][0]+" , "+p[1][k][1]);
			float cosAngle = Tools.angularDeviation(locX, locY, p[0][k][0], p[0][k][1], p[1][k][0], p[1][k][1]);
			logWriter.println("mapping # "+k+" cos(ang) = "+cosAngle+" : "+peakAngsAtLoc.get(0).get(map[k][0])+" \t "+peakAngsAtLoc.get(1).get(map[k][1]));

			if (cosAngle>MIN_COS_ANG &&
						isForeground(p[0][k], inip1, backgr1, margin1) &&
						isForeground(p[1][k], inip1, backgr1, margin1)) {
				idxDirectedClusters.add(k); // contains mapping index k
			}
		}

		logWriter.println("found "+idxDirectedClusters.size()+" aligned foreground clusters out of matched "+map.length+" ... ");

		if (idxDirectedClusters.size()<3) {
			return 0; // idxDirectedClusters.size() how many matched & directed
		}

		// check directions for each ring
		float[][] directionsRad = new float[2][idxDirectedClusters.size()];  // at least 2x3, 2x4
		for (int d=0; d<idxDirectedClusters.size(); d++) {
				// ring A
				int k = idxDirectedClusters.get(d); // map[k][0] is index in A, map[k][1] is index in B
				directionsRad[0][d] = peakAngsAtLoc.get(0).get(map[k][0]) * Deg2Rad;
				directionsRad[1][d] = peakAngsAtLoc.get(1).get(map[k][1]) * Deg2Rad;
		}
		// don't allow any inter-angle to be more than 180 deg.
		Arrays.sort(directionsRad[0]);
		Arrays.sort(directionsRad[1]);

		logWriter.println("checked directions in [rad] : ");
		for (int i=0; i< directionsRad.length; i++) {
			logWriter.println("ring "+i+" : ");
			for (int j=0; j< directionsRad[0].length; j++) {
				logWriter.print(directionsRad[i][j]+" ");
			}
			logWriter.println();
		}

		if (!sortedDirectionsAnglesAreLessThanPI(directionsRad[0])) {
			logWriter.println("ring A had two directions with illegal angle => NO");
			return 0;
		}
		if (!sortedDirectionsAnglesAreLessThanPI(directionsRad[1])) {
			logWriter.println("ring B had two directions with illegal angle => NO");
			return 0;
		}

		// check each ring profile
		boolean chk;
		float r_A = (float) (neuronDiamMax*scale_A);
		chk = profileOK(locX, locY, r_A, directionsRad[0], inip1, backgr1, margin1);
		if (!chk) {
			logWriter.println("profile A was illegal");
			return 0;
		}
		float r_B = (float) (neuronDiamMax*scale_B);
		chk = profileOK(locX, locY, r_B, directionsRad[1], inip1, backgr1, margin1);
		if (!chk) {
			logWriter.println("profile B was illegal");
			return 0;
		}

		logWriter.println(" YES! \n");

		/*******************/

		return 255;

    }

	private boolean isForeground(float[] loc, FloatProcessor inip1, FloatProcessor back1, FloatProcessor marg1)
	{
		return Interpolator.interpolateAt(loc[0], loc[1], inip1) > Interpolator.interpolateAt(loc[0], loc[1], back1) + Interpolator.interpolateAt(loc[0], loc[1], marg1);
	}

	private boolean profileOK(
										float posx,
										float posy,
										float posr,
										float[] dirAngRadSorted,
										FloatProcessor inip1,
										FloatProcessor backgr1,
										FloatProcessor margin1)
	{
		float DF = 0.001f;
		float oStep = (float) (2*Math.PI)/100f;
		int totalPeaks = dirAngRadSorted.length;

		boolean up, down; // check if it went below th. between peaks
		float currX, currY, currI, currB, currM;

		// directionAnglesRadSorted[0]
		up = false;
		down = true;
		for (float aa=dirAngRadSorted[0]; aa>=dirAngRadSorted[dirAngRadSorted.length-1]-2*Math.PI+DF; aa-=oStep) {
			currX = (float) (posx + posr * Math.cos(aa));
			currY = (float) (posy - posr * Math.sin(aa));
			currI = Interpolator.interpolateAt(currX, currY, inip1);
			currB = Interpolator.interpolateAt(currX, currY, backgr1);
			currM = Interpolator.interpolateAt(currX, currY, margin1);
			up = up || (currI>currB+currM);
			down = up && (currI<currB+currM); //Masker.HYSTERESIS
			if (down) break;
		}
		if (!down) {
			logWriter.println("between "+(dirAngRadSorted[0]*Rad2Deg)+" and "+(dirAngRadSorted[dirAngRadSorted.length-1]*Rad2Deg));
			return false; // don't go further
		}

		// the rest
		for (int peakIdx=1; peakIdx<totalPeaks; peakIdx++) {
			up = false;
			down = true;
			for (float aa=dirAngRadSorted[peakIdx]; aa>=dirAngRadSorted[peakIdx-1]+DF; aa-=oStep) {
				currX = (float) (posx + posr * Math.cos(aa));
				currY = (float) (posy - posr * Math.sin(aa));
				currI = Interpolator.interpolateAt(currX, currY, inip1);
				currB = Interpolator.interpolateAt(currX, currY, backgr1);
				currM = Interpolator.interpolateAt(currX, currY, margin1);
				up = up || (currI>currB+currM);
				down = up && (currI<currB+currM); // Masker.HYSTERESIS
				if (down) break;
			}
			if (!down) {
				logWriter.println("between "+(dirAngRadSorted[peakIdx]*Rad2Deg)+" and "+(dirAngRadSorted[peakIdx-1]*Rad2Deg));
				return false; // don't go further
			}
		}

		return true;

	}

	private boolean sortedDirectionsAnglesAreLessThanPI(float[] angsRad)
	{

		if (angsRad[0]-angsRad[angsRad.length-1]+(2*Math.PI)>Math.PI) {
			return false;
		}

		for (int i=1; i<angsRad.length; i++) {
			if (angsRad[i]-angsRad[i-1]>Math.PI) {
				return false;
			}
		}

		return true;

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

        // pressing key 'U' will export data to profile.log

        if (e.getKeyCode() == KeyEvent.VK_U) {

            if (pwP_A!=null) {

                // export to csv when clicked
                String profileFile = "profile_v11.log";

                // empty the file
                PrintWriter writer = null;
                try {
                    writer = new PrintWriter(profileFile);  writer.print("");   writer.close();
                } catch (FileNotFoundException ex) {}

                try {
                    PrintWriter out;

                    out = new PrintWriter(new BufferedWriter(new FileWriter(profileFile, true)));

                    // profile A
                    out.print("angle[deg]");for (int idx=0; idx<profile_A.length; idx++)    out.print(" " + ang_A[idx]);                out.print("\n");
                    out.print("profile_A"); for (int idx=0; idx<profile_A.length; idx++)    out.print(" " + profile_A[idx]);            out.print("\n");
                    out.print("peak_A_ang");for (int idx=0; idx<peakIdx_A.length; idx++)    out.print(" " + peakAng_A[idx] * Rad2Deg);  out.print("\n");
                    out.print("peak_A_val");for (int idx=0; idx<peakIdx_A.length; idx++)    out.print(" " + Tools.interp1Darray(peakIdx_A[idx], profile_A)); out.print("\n");

                    // profile B
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
