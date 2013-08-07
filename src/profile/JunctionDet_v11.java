package profile;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.*;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
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
    double scale;               // only one scale now
    double neuronDiamMax, D;
    ArrayList<ArrayList<float[]>>   profiles;
    PlotWindow pw;

    // some variables
    float[] exProf;
    float[] xMS;
    float[] yMS;
    float[] xCS;
    float[] yCSlower;
    float[] yCSupper;

    int CPU_NR;

    public void mouseClicked(MouseEvent e) {

        int offscreenX = canvas.offScreenX(e.getX());
        int offscreenY = canvas.offScreenY(e.getY());

        exProf = Profiler.extractProfile(offscreenX, offscreenY);
        Analyzer.extractPeakIdxs(exProf);

        // Fill in X axis (frame number)
        float[] x = new float[exProf.length];
        for (int i = 1; i <= x.length; i++) x[i - 1] = (i-1)*Profiler.resolDeg;

        // Prepare plot window
        Plot chart = new Plot("", "", "", x, exProf);
        chart.setSize(900, 450);

        if (pw == null) {
            pw = chart.show();
            //pw.addWindowListener(this);
        } else
            pw.setTitle("Profile, x = " + offscreenX + ", y = " + offscreenY);

        // Add the points for prettier plots
        chart.addPoints(x, exProf, PlotWindow.CIRCLE);

        // Add MS convergence points
        xMS = new float[Analyzer.nrPoints];
        yMS = new float[Analyzer.nrPoints];
        float[] yLW;
        for (int i=0; i<Analyzer.nrPoints; i++) {
            xMS[i] = Analyzer.convIdx.get(0).get(0)[i] * Profiler.resolDeg;
            yMS[i] = (float) Tools.interp1Darray(Analyzer.convIdx.get(0).get(0)[i], exProf);
        }
        chart.draw();
        chart.setColor(Color.GREEN);
        chart.addPoints(xMS, yMS, Plot.X);

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

        pw.drawPlot(chart);
        pw.setTitle("Profile, x = " + offscreenX + ", y = " + offscreenY + " : domes -> " +domes);

    }

    @Override
    public void mousePressed(MouseEvent e) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void mouseReleased(MouseEvent e) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void mouseEntered(MouseEvent e) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void mouseExited(MouseEvent e) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public int setup(String s, ImagePlus imagePlus) {
        if(imagePlus==null) return DONE;
        imp = Tools.convertToFloatImage(imagePlus);
        imp.setTitle("inimg");
        canvas = imagePlus.getCanvas();
        return DOES_8G+DOES_32+NO_CHANGES;
    }

    public void run(ImageProcessor imageProcessor) {

        int totalJobs; // used in paralelization lines
        long t1, t2;
        scale = 2.0;

        neuronDiamMax       =  Prefs.get("advantra.critpoint.neuronDiamMax", 3);
        scale               =  Prefs.get("advantra.critpoint.scale", 1.5);
        D                   =  Prefs.get("advantra.critpoint.D", 10);
        CPU_NR              = (int) Prefs.get("advantra.critpoint.CPU_NR", 4);

        GenericDialog gd = new GenericDialog("JUNCTION DET.");
        gd.addNumericField("neuronD ",  neuronDiamMax, 0, 10, " MAX");
        gd.addNumericField("scale ",  scale, 1, 10, " x(neuronD)");
        gd.addNumericField("D",         D, 1, 10, "score param.");
        gd.addNumericField("CPU_NR ",   CPU_NR, 0, 10, "");

        gd.showDialog();
        if (gd.wasCanceled()) return;

        neuronDiamMax       =  gd.getNextNumber();
        Prefs.set("advantra.critpoint.neuronDiamMax",   neuronDiamMax);
        scale               = gd.getNextNumber();
        Prefs.set("advantra.critpoint.scale",   scale);
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
        Masker.loadTemplate(imp.getProcessor(), neighbourhoodR, (float) D);
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

        //scoreimg = new ByteProcessor(imp.getWidth(), imp.getHeight());
        Profiler.loadTemplate(imp.getProcessor(), Masker.maskip);
        // can see the locations now
        int totalLocations = Profiler.locations.length;

        IJ.log("total "+totalLocations+" locations given by mask ("+(100*(float)totalLocations/(imp.getWidth()*imp.getHeight()))+" % kept)");

        profiles     	= new ArrayList<ArrayList<float[]>>(totalLocations);

        /*
        ***********************************************************
         */

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

        IJ.log("calculating peaks for selected locations... ");
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

        ImagePlus imageLabels = conn_reg.showLabels();
        imageLabels.show();

        /*
        -----------------------------------------------
         */

        canvas.setOverlay(formOverlay(conn_reg.getConnectedRegions()));
        canvas.addMouseListener(this);
        canvas.addMouseMotionListener(this);
        canvas.addKeyListener(this);
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

        // Analyzer, Profiler
        //IJ.log("extracting score, analyser has " + Analyzer.peakIdx.length + " points " + Analyzer.peakIdx[0][0].length);
        //IJ.log("extracting score, profiler has " + Profiler.profiles.length + " locations with " + Profiler.profiles[0].length + " points");
        double score = 0;

        if (Analyzer.peakIdx[locIdx][0] != null) {

            float[] domes = new float[3];
            int[] roundedPks = new int[]{Math.round(Analyzer.peakIdx[locIdx][0][0]), Math.round(Analyzer.peakIdx[locIdx][0][1]), Math.round(Analyzer.peakIdx[locIdx][0][2])};

            for (int i = 0; i < 3; i++) {
                float pH, pL;
                pH = (float) Tools.interp1Darray(Analyzer.peakIdx[locIdx][0][i], Profiler.profiles[locIdx]);
                pL = Tools.findNextStationaryValue(Analyzer.peakIdx[locIdx][0][i], roundedPks, Profiler.profiles[locIdx]);
                domes[i] = (pH>pL)?(pH-pL):0;
            }

            if (Tools.min3(domes)>D) {
                score = 255;
            }

        }

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

    @Override
    public void mouseDragged(MouseEvent e) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public void mouseMoved(MouseEvent e) {
        mouseClicked(e);
    }

    public void keyTyped(KeyEvent e) {

    }

    public void keyPressed(KeyEvent e) {


        if (e.getKeyCode() == KeyEvent.VK_U) {

            if (pw!=null) {

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
    public void keyReleased(KeyEvent e) {
        //To change body of implemented methods use File | Settings | File Templates.
    }
}
