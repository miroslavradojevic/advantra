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
    PlotWindow pwP_A;
	PlotWindow pwI_B;

    Overlay     currOvl = new Overlay();
    PointRoi    pks     = new PointRoi(0, 0);

    // some variables
    //float[] exProf;
//    float[] xMS;
//    float[] yMS;
//    float[] xCS;
//    float[] yCSlower;
//    float[] yCSupper;

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
//		D                   =  Prefs.get("advantra.critpoint.D", 10);

		GenericDialog gd = new GenericDialog("JUNCTION DET.");
		gd.addNumericField("neuronD ",  neuronDiamMax, 0, 10, " MAX");
		gd.addNumericField("scale ",  scale, 1, 10, " x(neuronD)");
		gd.addChoice("background extraction", new String[]{"MEAN", "MEDIAN"}, bgExtractionMode);
//		gd.addNumericField("D",         D, 1, 10, "score param.");

		gd.showDialog();
		if (gd.wasCanceled()) return;

		neuronDiamMax       =  gd.getNextNumber();
		Prefs.set("advantra.critpoint.neuronDiamMax",   neuronDiamMax);

		scale               = gd.getNextNumber();
		Prefs.set("advantra.critpoint.scale",   scale);

		bgExtractionMode      = gd.getNextChoice();
		Prefs.set("advantra.critpoint.mask.bgExtractionMode",   bgExtractionMode);

//		D                   =  gd.getNextNumber();
//		Prefs.set("advantra.critpoint.D",               D);

		CPU_NR = Runtime.getRuntime().availableProcessors();


        boolean doCalculations = false;

        /*
        ***********************************************************
         */
		IJ.log("excluding background...");

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
		inmask.show();

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
        //canvas.setOverlay(currOvl);
        //canvas.getImage().updateImage();
	}

    public void mouseClicked(MouseEvent e) {

		// on mouse click extract locations, plot profile, intensities along profile and overlay points
        currOvl.clear();
        int offscreenX = canvas.offScreenX(e.getX());
        int offscreenY = canvas.offScreenY(e.getY());

        // define ring A
        double rd_A = neuronDiamMax*scale;

        float[] profile_A = Profiler.extractProfile(neuronDiamMax, scale, offscreenX, offscreenY, (FloatProcessor) imp.getProcessor());
        float[] profile_A_MinMax = Tools.getMinMax(profile_A);
        float[] peakIdx_A = Analyzer.extractPeakIdxs(profile_A);

		float[] peakAng_A, peakX_A, peakY_A, peakI_A, peakB_A;
        if (peakIdx_A==null) { peakAng_A = peakX_A = peakY_A = peakI_A = peakB_A = null; }
        else {
            peakAng_A = peakX_A = peakY_A = peakI_A = peakB_A = new float[peakIdx_A.length];

            for (int i=0; i<peakIdx_A.length; i++) {
                peakAng_A[i] = peakIdx_A[i] * Profiler.getResolDeg(scale) * Deg2Rad;//(float)(Math.PI/180f);
                peakX_A[i] = (float) (offscreenX + rd_A * Math.cos( peakAng_A[i] ));
                peakY_A[i] = (float) (offscreenY - rd_A * Math.sin( peakAng_A[i] ));
                peakI_A[i] = Interpolator.interpolateAt(peakX_A[i], peakY_A[i], (FloatProcessor) imp.getProcessor());
                peakB_A[i] = Interpolator.interpolateAt(peakX_A[i], peakY_A[i], (FloatProcessor) Masker.back);
                if (peakI_A[i]-peakB_A[i] > Masker.VISIBLE_INTENSITY_DIFF) {
                    currOvl.add(new PointRoi(peakX_A[i]+.5, peakY_A[i]+.5));
                }

            }
        }

		Plot chartP_A, chartI_B;



        // define ring B

        /*
        plots
         */

        // Fill in X axis (frame number)
        float[] x   = new float[profile_A.length];
        float[] xI  = new float[profile_A.length];
        float[] xB  = new float[profile_A.length];
        for (int i = 1; i <= x.length; i++) {
            x[i - 1] = (i-1)*Profiler.getResolDeg(scale);
            float pX = (float) (offscreenX + rd * Math.cos( x[i - 1] * ((float)Math.PI/180f) ));
            float pY = (float) (offscreenY - rd * Math.sin( x[i - 1] * ((float)Math.PI/180f) ));
            xI[i-1]  = Interpolator.interpolateAt(pX, pY, (FloatProcessor) imp.getProcessor());
            xB[i-1]  = Interpolator.interpolateAt(pX, pY, (FloatProcessor) Masker.back);
        }

        Plot chart1 = new Plot("", "", "", x, xI);
        chart1.setSize(900, 450);
        chart1.draw();  chart1.setColor(Color.BLUE);
        chart1.addPoints(x, xB, PlotWindow.LINE);

        chart1.draw();  chart1.setColor(Color.RED);
        if (peakAng!=null) {
            for (int i=0; i<peakAng.length; i++) {
                chart1.drawLine(peakAng[i]*Rad2Deg, peakI[i], peakAng[i]*Rad2Deg, peakB[i]);
            }
        }

        if (pwI == null) {
            pwI = chart1.show();
        }

        // Prepare plot window
        Plot chart = new Plot("", "", "", x, exProf);
        chart.setSize(900, 450);

        if (pwP == null) pwP = chart.show();
        pwP.drawPlot(chart);
        pwP.setTitle("FilterProfile, x = " + offscreenX + ", y = " + offscreenY);

        chart.addPoints(x, exProf, PlotWindow.CIRCLE);
        chart.draw();  chart1.setColor(Color.RED);
        if (peakAng!=null) {
            for (int i=0; i<peakAng.length; i++) {
                chart.drawLine(peakAng[i]*Rad2Deg, profileMinMax[0], peakAng[i]*Rad2Deg, profileMinMax[1]);
            }
        }


        currOvl.add(new OvalRoi(offscreenX-rd+.5, offscreenY-rd+.5, 2*rd, 2*rd));
        float cI = Interpolator.interpolateAt(offscreenX, offscreenY, (FloatProcessor) imp.getProcessor());
        float cB = Interpolator.interpolateAt(offscreenX, offscreenX, (FloatProcessor) Masker.back);
        if (cI-cB>Masker.VISIBLE_INTENSITY_DIFF) {
            PointRoi pt = new PointRoi(offscreenX+.5, offscreenY+.5);
            pt.setStrokeColor(Color.RED);
            currOvl.add(pt);

        }
        else
            currOvl.add(new PointRoi(offscreenX + .5, offscreenY + .5));

        canvas.setOverlay(currOvl);


/*
        //Add MS plot
        String domes = "";
        if (Analyzer.peakIdx[0][0]!=null) { // >=3 peaks
            chart.draw();
            chart.setColor(Color.RED);
            chart.setLineWidth(5);
            xCS         = new float[3];
            yCSlower    = new float[3];
            yCSupper    = new float[3];
            int[] roundedPeaks = new int[]{
                    Math.round(Analyzer.peakIdx[0][0][0]),
                    Math.round(Analyzer.peakIdx[0][0][1]),
                    Math.round(Analyzer.peakIdx[0][0][2])
            };
            for (int i=0; i<3; i++) {
                xCS[i] = Analyzer.peakIdx[0][0][i] * Profiler.resolDeg;
                yCSupper[i] = (float) Tools.interp1Darray(Analyzer.peakIdx[0][0][i], exProf);
                yCSlower[i] = Tools.findNextStationaryValue(Analyzer.peakIdx[0][0][i], roundedPeaks, exProf);
                chart.drawLine(xCS[i], yCSlower[i], xCS[i], yCSupper[i]);
                domes += "   " + IJ.d2s((yCSupper[i]>yCSlower[i])?(yCSupper[i]-yCSlower[i]):0, 2) + "  ";
            }
            //chart.addPoints(xMS, yMS, Plot.BOX);

        }


        */


        pwI.drawPlot(chart1);
        pwI.setTitle("IntensityProfile, x = " + offscreenX + ", y = " + offscreenY);

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

            if (pwP!=null) {

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
                        out.println(aDeg+", "+exProf[idx++]);
                    }
                    out.close();

                    out = new PrintWriter(new BufferedWriter(new FileWriter(msFile, true)));
                    out.println("Angle, Response");

                    for (int i = 0; i<xMS.length; i++) {
                        out.println(xMS[i]+", "+yMS[i]);
                    }
                    out.close();

                    out = new PrintWriter(new BufferedWriter(new FileWriter(csFile, true)));
                    out.println("Angle, Lower, Higher");

                    for (int i = 0; i<xCS.length; i++) {
                        out.println(xCS[i]+", "+yCSlower[i]+", "+yCSupper[i]);
                    }
                    out.close();

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
