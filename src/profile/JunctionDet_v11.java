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
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 8/7/13
 * Time: 9:45 AM
 */
public class JunctionDet_v11 implements PlugInFilter, MouseListener, MouseMotionListener, KeyListener {

    ImagePlus imp;
    ImageCanvas canvas;
	String inimgTitle;
    double scale;               // only one scale now
    double neuronDiamMax, D;
    ArrayList<ArrayList<float[]>>   profiles;

    PlotWindow pwP_A, pwI_A, pwP_B, pwI_B;

    Overlay     currOvl = new Overlay();

    // some variables (used by .csv export (buttonPressed) and mouse click plots (mouseCLicked))
    float[] profile_A, profile_B;

    int CPU_NR;
	String bgExtractionMode;

    private static float Deg2Rad = (float) (Math.PI/180f);
    private static float Rad2Deg = (float) (180f/Math.PI);

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
		scale = 2.0;

		neuronDiamMax       =  Prefs.get("advantra.critpoint.neuronDiamMax", 3);
		scale               =  Prefs.get("advantra.critpoint.scale", 1.5);
		bgExtractionMode        =           Prefs.get("advantra.critpoint.mask.bgExtractionMode", "MEAN");

		GenericDialog gd = new GenericDialog("JUNCTION DET.");
		gd.addNumericField("neuronD ",  neuronDiamMax, 0, 10, " MAX");
		gd.addNumericField("scale ",  scale, 1, 10, " x(neuronD)");
		gd.addChoice("background extraction", new String[]{"MEAN", "MEDIAN"}, bgExtractionMode);

		gd.showDialog();
		if (gd.wasCanceled()) return;

		neuronDiamMax       =  gd.getNextNumber();
		Prefs.set("advantra.critpoint.neuronDiamMax",   neuronDiamMax);

		scale               = gd.getNextNumber();
		Prefs.set("advantra.critpoint.scale",   scale);

		bgExtractionMode      = gd.getNextChoice();
		Prefs.set("advantra.critpoint.mask.bgExtractionMode",   bgExtractionMode);

		CPU_NR = Runtime.getRuntime().availableProcessors();

        boolean doCalculations = false;

        /*
        ***********************************************************
         */
		IJ.log("extracting background...");

		t1 = System.currentTimeMillis();

		int neighbourhoodR = (int) Math.ceil(scale*neuronDiamMax);

		if (bgExtractionMode=="MEAN") {
			Masker.loadTemplate(imp.getProcessor(), neighbourhoodR, 0, true);
		}
		else if (bgExtractionMode=="MEDIAN") {
			Masker.loadTemplate(imp.getProcessor(), neighbourhoodR, 1, true);
		}
		else {
			IJ.log("extraction mode was wrong");
			return;
		}

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
		IJ.log("done "+((t2-t1)/1000f)+" sec.");

		ImagePlus inmask = new ImagePlus("inmask", Masker.maskip);
		//inmask.show();

		int nr = 0;
		for (int i=0; i<totalJobs; i++) if (Masker.maskip.getf(i)==255) nr++;

		IJ.log("done. total "+nr+" check locations given by mask ("+(100*(float)nr/(totalJobs))+" % kept)");

		/*
        ***********************************************************
         */
        if (doCalculations) {

		GenericDialog gd1 = new GenericDialog("PRUNE REGIONS?");
		//gd1.addMessage("extracted "+nr+" out of "+totalProfiles+"  ("+IJ.d2s(nr*100f/totalProfiles, 2)+"%)");
		gd1.addMessage("PRUNE REGIONS?");

		gd1.showDialog();
		if (!gd1.wasCanceled()) {
			IJ.log("pruning...");

			/*
			***********************************
			 */


			// get connected regions
			t1 = System.currentTimeMillis();
			IJ.log("conn. regions...");
			Find_Connected_Regions conn_reg = new Find_Connected_Regions(inmask, true);
			conn_reg.run("");
			t2 = System.currentTimeMillis();

			int nr_regions = conn_reg.getNrConnectedRegions();

			IJ.log("\n"+nr_regions+" connected regions extracted.\n" +
						   "elapsed: "+((t2-t1)/1000f)+ " seconds.");

			ArrayList<ArrayList<int[]>> regs = conn_reg.getConnectedRegions();

			IJ.log("pruning..."); // TODO: mask pruning inside Masker class, here call method

			for (int i=0; i<regs.size(); i++) {

				if (regs.get(i).size()>Masker.VERY_SMALL_REGION_SIZE) {

					float regAvg = 0;
					float bkgAvg = 0;
					for (int w=0; w<regs.get(i).size(); w++) {
						int takeValX = regs.get(i).get(w)[1];
						int takeValY = regs.get(i).get(w)[0];
						regAvg += imp.getProcessor().getf(takeValX, takeValY);
						bkgAvg += Masker.back.getf(takeValX, takeValY);
					}
					regAvg /= regs.get(i).size();
					bkgAvg /= regs.get(i).size();

					if (regAvg < bkgAvg + Masker.VISIBLE_INTENSITY_DIFF ) {
						for (int q=0; q<regs.get(i).size(); q++) {
							int removeX = regs.get(i).get(q)[1];
							int removeY = regs.get(i).get(q)[0];
							inmask.getProcessor().set(removeX, removeY, 0);
						}
					}
				}
				else {  // prune very small regions
					for (int q=0; q<regs.get(i).size(); q++) {
						int removeX = regs.get(i).get(q)[1];
						int removeY = regs.get(i).get(q)[0];
						inmask.getProcessor().set(removeX, removeY, 0);
					}
				}

			}


			/*
			***********************************
			 */


			IJ.log("done.");
		}

		inmask.updateAndDraw();
		IJ.selectWindow(inimgTitle);
		//IJ.run("Add Image...", "image=inmask x="+0+" y="+0+" opacity=40");
		//IJ.selectWindow(inimgTitle);

		canvas.zoomIn(0, 0);

		// check amount of extracted locations
		nr = 0;
		for (int i=0; i<totalJobs; i++) if (Masker.maskip.getf(i)==255) nr++;
		IJ.log("extracted "+nr+" out of "+totalJobs+"  ("+IJ.d2s(nr*100f/totalJobs, 2)+"%)");

        /*
        ***********************************************************
         */

		Profiler.loadTemplate(imp.getProcessor(), Masker.maskip);
		int totalLocations = Profiler.locations.length;
		profiles     	= new ArrayList<ArrayList<float[]>>(totalLocations);

		IJ.log("calculating profiles... ");
		t1 = System.currentTimeMillis();
//        double R = neuronDiamMax*scale;
		Profiler.loadParams(neuronDiamMax, scale, false);

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

		updateProfilesList();

		t2 = System.currentTimeMillis();
		IJ.log("done "+((t2-t1)/1000f)+" sec.");

        /*
        ***********************************************************
         */

		IJ.log("calculating peaks for selected locations... "+CPU_NR);
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

		ByteProcessor scoreimg = new ByteProcessor(imp.getWidth(), imp.getHeight());

		for (int locIdx=0; locIdx<totalLocations; locIdx++) {

			double score = 0;

			if (Analyzer.peakIdx[locIdx][0] != null) {

				//score =  extractScore(locIdx, false);

			}

			scoreimg.set(Profiler.locations[locIdx][0], Profiler.locations[locIdx][1], ((int) score));

		}

		//ImagePlus scoreImagePlus = new ImagePlus("score", scoreimg);
		//scoreImagePlus.show();
		IJ.selectWindow(inimgTitle);
            IJ.setTool("hand");

        /*
        -----------------------------------------------
         */


/*		t1 = System.currentTimeMillis();
		Find_Connected_Regions conn_reg = new Find_Connected_Regions(scoreImagePlus, true);
		conn_reg.run("");
		t2 = System.currentTimeMillis();
		int nr_regions = conn_reg.getNrConnectedRegions();

		IJ.log(nr_regions+" connected regions extracted.\n" +
					   "elapsed: "+((t2-t1)/1000f)+ " seconds.");

		//ImagePlus imageLabels = conn_reg.showLabels();
		//imageLabels.show();*/


        }//skipCalculations
        /*
        -----------------------------------------------
         */

		//canvas.setOverlay(formOverlay(conn_reg.getConnectedRegions()));
		canvas.addMouseListener(this);
		canvas.addMouseMotionListener(this);
		canvas.addKeyListener(this);

        IJ.selectWindow(inimgTitle);
        IJ.setTool("hand");
	}

    public void mouseClicked(MouseEvent e) {

		// on mouse click extract locations, plot profile, intensities along profile and overlay points
        currOvl.clear();

        int offscreenX = canvas.offScreenX(e.getX());
        int offscreenY = canvas.offScreenY(e.getY());

//        ArrayList<Float> C_x    = new ArrayList<Float>();
//        ArrayList<Float> C_y    = new ArrayList<Float>();

        float cI = Interpolator.interpolateAt(offscreenX, offscreenY, (FloatProcessor) imp.getProcessor());
        float cB = Interpolator.interpolateAt(offscreenX, offscreenY, (FloatProcessor) Masker.back);
        if (cI-cB > Masker.VISIBLE_INTENSITY_DIFF) {
            PointRoi pt = new PointRoi(offscreenX+.5, offscreenY+.5);
            //pt.setStrokeColor(Color.YELLOW);
            currOvl.add(pt);
//            C_x.add(Float.valueOf(offscreenX));
//            C_y.add(Float.valueOf(offscreenY));
        }

        // define ring A
        double scale_A = scale;
        double rd_A = neuronDiamMax*scale_A;
        OvalRoi ring_A = new OvalRoi(offscreenX-rd_A+.5, offscreenY-rd_A+.5, 2*rd_A, 2*rd_A);
        currOvl.add(ring_A);

        profile_A = Profiler.extractProfile(neuronDiamMax, scale_A, offscreenX, offscreenY, (FloatProcessor) imp.getProcessor());
        float[] peakIdx_A = Analyzer.extractPeakIdxs(profile_A); // MS

        ArrayList<Float> A_Ang  = new ArrayList<Float>(); // store those that were selected as ok
        ArrayList<Float> A_x    = new ArrayList<Float>();
        ArrayList<Float> A_y    = new ArrayList<Float>();

		float[] peakAng_A   = null;
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
                peakB_A[i] = Interpolator.interpolateAt(peakX_A[i], peakY_A[i], (FloatProcessor) Masker.back);

                if (peakI_A[i]-peakB_A[i] > Masker.VISIBLE_INTENSITY_DIFF) {  // criteria!!
                    //PointRoi pt = new PointRoi(peakX_A[i]+.5, peakY_A[i]+.5);
                    //pt.setStrokeColor(getColor(i));
                    //currOvl.add(pt);
                    A_Ang.add(peakAng_A[i] * Rad2Deg); // because hungarian matching method will take angle differences in degrees
                    A_x.add(peakX_A[i]);
                    A_y.add(peakY_A[i]);
                }

            }
        }

        /*
        visualizations
         */

		Plot chartP_A, chartI_A;
        float[] profile_A_MinMax = Tools.getMinMax(profile_A);

        // intensity profile
        float[] x   = new float[profile_A.length];
        float[] xI  = new float[profile_A.length];
        float[] xB  = new float[profile_A.length];
        for (int i = 1; i <= x.length; i++) {
            x[i - 1] = (i-1)*Profiler.getResolDeg(scale_A);
            float pX = (float) (offscreenX + rd_A * Math.cos( x[i - 1] * ((float)Math.PI/180f) ));
            float pY = (float) (offscreenY - rd_A * Math.sin( x[i - 1] * ((float)Math.PI/180f) ));
            xI[i-1]  = Interpolator.interpolateAt(pX, pY, (FloatProcessor) imp.getProcessor());
            xB[i-1]  = Interpolator.interpolateAt(pX, pY, (FloatProcessor) Masker.back);
        }

        chartP_A = new Plot("", "", "", x, profile_A);
        chartP_A.setSize(600, 300);
        chartP_A.addPoints(x, profile_A, PlotWindow.CIRCLE);
        chartP_A.draw();  chartP_A.setColor(Color.RED);
        // draw peak detections
        if (peakAng_A!=null)
            for (int i = 0; i < peakIdx_A.length; i++)
                chartP_A.drawLine(peakAng_A[i] * Rad2Deg, profile_A_MinMax[0], peakAng_A[i] * Rad2Deg, Tools.interp1Darray(peakIdx_A[i], profile_A));

        chartI_A = new Plot("", "", "", x, xI);
        chartI_A.setSize(600, 300);
        chartI_A.draw(); chartI_A.setColor(Color.BLUE);
        chartI_A.addPoints(x, xB, PlotWindow.LINE);
        chartI_A.draw();  chartI_A.setColor(Color.RED);
        // draw peak detections
        if (peakAng_A!=null)
            for (int i = 0; i < peakIdx_A.length; i++)
                chartI_A.drawLine(peakAng_A[i] * Rad2Deg, Tools.interp1Darray(peakIdx_A[i], xB), peakAng_A[i] * Rad2Deg, Tools.interp1Darray(peakIdx_A[i], xI));

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
        double scale_B = scale_A + 1;
        double rd_B = neuronDiamMax*scale_B;
        OvalRoi ringB = new OvalRoi(offscreenX-rd_B+.5, offscreenY-rd_B+.5, 2*rd_B, 2*rd_B);
        currOvl.add(ringB);

        profile_B = Profiler.extractProfile(neuronDiamMax, scale_B, offscreenX, offscreenY, (FloatProcessor) imp.getProcessor());
        float[] peakIdx_B = Analyzer.extractPeakIdxs(profile_B); // MS

        ArrayList<Float> B_Ang  = new ArrayList<Float>(); // store those that were selected as ok
        ArrayList<Float> B_x    = new ArrayList<Float>();
        ArrayList<Float> B_y    = new ArrayList<Float>();

        float[] peakAng_B   = null;
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
                peakB_B[i] = Interpolator.interpolateAt(peakX_B[i], peakY_B[i], (FloatProcessor) Masker.back);

                if (peakI_B[i]-peakB_B[i] > Masker.VISIBLE_INTENSITY_DIFF) {  // criteria!! maybe use xB averages here
                    //PointRoi pt = new PointRoi(peakX_B[i]+.5, peakY_B[i]+.5);
                    //pt.setStrokeColor(getColor(i));
                    //currOvl.add(pt);
                    B_Ang.add(peakAng_B[i] * Rad2Deg); // because hungarian matching method will take angle differences in degrees
                    B_x.add(peakX_B[i]);
                    B_y.add(peakY_B[i]);
                }

            }
        }

        /*
        visualizations
         */

        Plot chartP_B, chartI_B;
        float[] profile_B_MinMax = Tools.getMinMax(profile_B);

        // intensity profile
        x   = new float[profile_B.length]; // use the same variable as for ring A
        xI  = new float[profile_B.length];
        xB  = new float[profile_B.length];
        for (int i = 1; i <= x.length; i++) {
            x[i - 1] = (i-1)*Profiler.getResolDeg(scale_B);
            float pX = (float) (offscreenX + rd_B * Math.cos( x[i - 1] * Deg2Rad ));
            float pY = (float) (offscreenY - rd_B * Math.sin( x[i - 1] * Deg2Rad ));
            xI[i-1]  = Interpolator.interpolateAt(pX, pY, (FloatProcessor) imp.getProcessor());
            xB[i-1]  = Interpolator.interpolateAt(pX, pY, (FloatProcessor) Masker.back);
        }

        chartP_B = new Plot("", "", "", x, profile_B);
        chartP_B.setSize(600, 300);
        chartP_B.addPoints(x, profile_B, PlotWindow.CIRCLE);
        chartP_B.draw();  chartP_B.setColor(Color.RED);
        // draw peak detections
        if (peakAng_B!=null)
            for (int i = 0; i < peakIdx_B.length; i++)
                chartP_B.drawLine(peakAng_B[i] * Rad2Deg, profile_B_MinMax[0], peakAng_B[i] * Rad2Deg, Tools.interp1Darray(peakIdx_B[i], profile_B));

        chartI_B = new Plot("", "", "", x, xI);
        chartI_B.setSize(600, 300);
        chartI_B.draw(); chartI_B.setColor(Color.BLUE);
        chartI_B.addPoints(x, xB, PlotWindow.LINE);
        chartI_B.draw();  chartI_B.setColor(Color.RED);
        // draw peak detections
        if (peakAng_B!=null)
            for (int i = 0; i < peakIdx_B.length; i++)
                chartI_B.drawLine(peakAng_B[i] * Rad2Deg, Tools.interp1Darray(peakIdx_B[i], xB), peakAng_B[i] * Rad2Deg, Tools.interp1Darray(peakIdx_B[i], xI));

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

//        IJ.log("***"+A_Ang.size()+" <->" + B_Ang.size());
//        for (int w=0; w<map.length; w++) {
//            IJ.log(Arrays.toString(map[w]));
//        }
//        IJ.log("***");

        // after mapping
        ArrayList<float[][]> clustersXY = new ArrayList<float[][]>(map.length);
        float[] angDiv = new float[map.length];
        // map[index from A][index from B]

        for (int k=0; k<map.length; k++) {

            float[][] ringLocsXY = new float[2][2];

            ringLocsXY[0][0] = A_x.get(map[k][0]);
            ringLocsXY[0][1] = A_y.get(map[k][0]); // point in A

            PointRoi pt = new PointRoi(ringLocsXY[0][0]+.5, ringLocsXY[0][1]+.5);
            pt.setStrokeColor(getColor(k));
            currOvl.add(pt);

            ringLocsXY[1][0] = B_x.get(map[k][1]);
            ringLocsXY[1][1] = B_y.get(map[k][1]); // point in B

            pt = new PointRoi(ringLocsXY[1][0]+.5, ringLocsXY[1][1]+.5);
            pt.setStrokeColor(getColor(k));
            currOvl.add(pt);

            clustersXY.add(k, ringLocsXY);

            angDiv[k] = Tools.angularDeviation(new float[]{offscreenX, offscreenY}, clustersXY.get(k)[0], clustersXY.get(k)[1]);

        }



        IJ.log(Tools.angularDeviation(new float[]{0,0}, new float[]{.5f,.5f}, new float[]{1,1})+" "+Arrays.toString(angDiv));

        /*
        **************************
         */

        canvas.setOverlay(currOvl);

    }

    private Overlay formOverlay(ArrayList<ArrayList<int[]>> regs)
    {

        Overlay detections = new Overlay();

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

//                IJ.log("A: "+A+" , B: "+B+"");
//                if(regs.get(i).get(0)[1]==584 && regs.get(i).get(0)[0]==254) {
//                    IJ.log("debug");
//                    IJ.log("x="+regs.get(i).get(0)[1]+", y="+regs.get(i).get(0)[0]+", A="+A+", B="+B+", angle="+ellipseParams[4]);
//                }

                float scalePlot = 2;
                detections.add(Tools.drawEllipse(ellipseParams[1]+0.5f, ellipseParams[0]+0.5f, scalePlot*A, scalePlot*B, ellipseParams[4], Color.YELLOW, 2, 50));

            }
        }

        return detections;

    }

    private double extractScore(int locIdx, boolean printVals) {

        double score = 0;

//		Analyzer.peakIdx[locIdx]

        //if (Analyzer.peakIdx[locIdx][0] != null) {

            float[] domes = new float[3];
            int[] roundedPks = new int[]{Math.round(Analyzer.peakIdx[locIdx][0][0]), Math.round(Analyzer.peakIdx[locIdx][0][1]), Math.round(Analyzer.peakIdx[locIdx][0][2])};

            // p1, p2, p3, w1, w2, w3 add them
            float p1 = (float) Tools.interp1Darray(Analyzer.peakIdx[locIdx][0][0], Profiler.profiles[locIdx]);
            float p2 = (float) Tools.interp1Darray(Analyzer.peakIdx[locIdx][0][1], Profiler.profiles[locIdx]);
            float p3 = (float) Tools.interp1Darray(Analyzer.peakIdx[locIdx][0][2], Profiler.profiles[locIdx]);

            float pL1 = Tools.findNextStationaryValue(Analyzer.peakIdx[locIdx][0][0], roundedPks, Profiler.profiles[locIdx]);
            float pL2 = Tools.findNextStationaryValue(Analyzer.peakIdx[locIdx][0][1], roundedPks, Profiler.profiles[locIdx]);
            float pL3 = Tools.findNextStationaryValue(Analyzer.peakIdx[locIdx][0][2], roundedPks, Profiler.profiles[locIdx]);

            float pLmin = Math.min(pL1, Math.min(pL2, pL3));

            float w1 = p1/(p1+p2+p3);
            float w2 = p2/(p1+p2+p3);
            float w3 = p3/(p1+p2+p3);

            score = w1*((p1-pL1)/(pL1-pLmin)) + w2*((p2-pL2)/(pL2-pLmin)) + w3*((p3-pL3)/(pL3-pLmin));

            /*

            for (int i = 0; i < 3; i++) {
                float pH, pL;
                pH = (float) Tools.interp1Darray(Analyzer.peakIdx[locIdx][0][i], Profiler.profiles[locIdx]);
                pL = Tools.findNextStationaryValue(Analyzer.peakIdx[locIdx][0][i], roundedPks, Profiler.profiles[locIdx]);
                domes[i] = (pH>pL)?(pH-pL):0;
            }

            if (Tools.min3(domes)>D) {  // th depends on the thickness
                score = 255;
            }
            */


        //}

        return score;

    }

    private void updateProfilesList() // uses Profiler to update  profiles, profilesName, angularRes
    {
        // store it in the list for each loc
        for (int loopLocations=0; loopLocations<Profiler.locations.length; loopLocations++) {

            // profileToAdd from Profiler.profiles
            float[] profileToAdd = new float[Profiler.profiles[loopLocations].length];
            for (int k=0; k<profileToAdd.length; k++) {
                profileToAdd[k] = Profiler.profiles[loopLocations][k];
            }

                ArrayList<float[]> dummy = new ArrayList<float[]>();
                dummy.add(profileToAdd);
                profiles.add(dummy);

        }
    }

    public void mouseMoved(MouseEvent e) {
        mouseClicked(e);
    }

    public void keyPressed(KeyEvent e) {


        if (e.getKeyCode() == KeyEvent.VK_U) {

            if (pwP_A!=null) {

                // export to csv
                String profileFile = "profile.csv", msFile = "ms.csv", csFile = "cs.csv";

                // empty the file
                PrintWriter writer = null;
                try {
                    writer = new PrintWriter(profileFile);  writer.print("");   writer.close();
                    writer = new PrintWriter(msFile);       writer.print("");   writer.close();
                    writer = new PrintWriter(csFile);    writer.print("");   writer.close();
                } catch (FileNotFoundException ex) {}

                try {
                    PrintWriter out;

                    out = new PrintWriter(new BufferedWriter(new FileWriter(profileFile, true)));
                    out.println("Angle, Response");
                    int idx = 0;
                    for (int aDeg = 0; aDeg<360; aDeg+=Profiler.resolDeg) {
                        out.println(aDeg+", "+profile_A[idx++]);
                    }
                    out.close();

//                    out = new PrintWriter(new BufferedWriter(new FileWriter(msFile, true)));
//                    out.println("Angle, Response");

//                    for (int i = 0; i<xMS.length; i++) {
//                        out.println(xMS[i]+", "+yMS[i]);
//                    }
//                    out.close();

//                    out = new PrintWriter(new BufferedWriter(new FileWriter(csFile, true)));
//                    out.println("Angle, Lower, Higher");

//                    for (int i = 0; i<xCS.length; i++) {
//                        out.println(xCS[i]+", "+yCSlower[i]+", "+yCSupper[i]);
//                    }
//                    out.close();

                } catch (IOException e1) {}

                IJ.log("exported to : \n" + new File(profileFile).getAbsolutePath()+ " \n " + new File(msFile).getAbsolutePath() + "\n" + new File(csFile).getAbsolutePath());

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
