package detection2d;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.*;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;

/**
 * Created by miroslav on 1/6/14.
 * Demo usage of PeakAnalyzer2D: uses also Masker2D, Profiler2D, and PeakExtractor2D shows delineated
 * branching structure extracted on all foreground locations as a joint Overlay on the image
 */
public class PeakAnalyzer2DDemo implements PlugInFilter, MouseListener, MouseMotionListener {

    float[][] 	inimg_xy;               // store input image as an array

    /*
    parameters
     */
    float       s = 1.5f;               // scale is fixed
    float       iDiff, D, minCos, scatterDist;
    int         M = 2;

    int         CPU_NR;

    /*
    interface things
     */
    ImagePlus pfl_im = new ImagePlus();
    ImageProcessor pfl_ip = null;
    ImageWindow pfl_iw;
    ImageCanvas cnv;

    public int setup(String s, ImagePlus imagePlus) {

        if(imagePlus==null) return DONE;

        inimg_xy = new float[imagePlus.getWidth()][imagePlus.getHeight()]; // x~column, y~row

        if (imagePlus.getType()== ImagePlus.GRAY8) {
            byte[] read = (byte[]) imagePlus.getProcessor().getPixels();
            for (int idx=0; idx<read.length; idx++) {
                inimg_xy[idx%imagePlus.getWidth()][idx/imagePlus.getWidth()] = (float) (read[idx] & 0xff);
            }

        }
        else if (imagePlus.getType()==ImagePlus.GRAY32) {
            float[] read = (float[]) imagePlus.getProcessor().getPixels();
            for (int idx=0; idx<read.length; idx++) {
                inimg_xy[idx%imagePlus.getWidth()][idx/imagePlus.getWidth()] = read[idx];
            }
        }
        else {
            IJ.log("image type not recognized");
            return DONE;
        }

        /******************************
         Generic Dialog
         *****************************/
        this.iDiff 					    	= (float)   Prefs.get("advantra.critpoint.mask.iDiff", 5);
        this.D 					        	= (float)   Prefs.get("advantra.critpoint.profile.d", 4);
        this.M 					        	= (int)     Prefs.get("advantra.critpoint.analyze.m", 2);
        this.minCos 					    = (float)   Prefs.get("advantra.critpoint.analyze.min_cos", 0.5f);
        this.scatterDist                    = (float)   Prefs.get("advantra.critpoint.analyze.scatter_d", 5f);

        GenericDialog gd = new GenericDialog("PROFILER2DDEMO");
        gd.addNumericField("iDiff ", 	iDiff, 			0, 10, "intensity margin");
        gd.addNumericField("D ", 	    D, 			    0, 10, "neuron diameter[pix]");
        gd.addNumericField("M ", 	    M, 			    0, 10, "branch len.");
        gd.addNumericField("min cos ", 	minCos, 	    1, 10, "alignment parameter.");
        gd.addNumericField("scatter d ",scatterDist, 	1, 10, "max allowed scatter (robustness test)");

        gd.showDialog();
        if (gd.wasCanceled()) return DONE;

        iDiff       	= (float) gd.getNextNumber();
        Prefs.set("advantra.critpoint.mask.iDiff", 	    	iDiff);

        D       	    = (float) gd.getNextNumber();
        Prefs.set("advantra.critpoint.profile.d", 	    	D);

        M       	    = (int) gd.getNextNumber();
        Prefs.set("advantra.critpoint.analyze.m", 	    	M);

        minCos       	= (float) gd.getNextNumber();
        Prefs.set("advantra.critpoint.analyze.min_cos", 	minCos);

        scatterDist     = (float) gd.getNextNumber();
        Prefs.set("advantra.critpoint.analyze.scatter_d", 	scatterDist);

        CPU_NR = Runtime.getRuntime().availableProcessors();

        cnv = imagePlus.getCanvas();

        return DOES_8G+DOES_32+NO_CHANGES;

    }

    public void run(ImageProcessor imageProcessor) {

        Sphere2D sph2d = new Sphere2D(D, s);

        ImagePlus samplingScheme =  sph2d.showSampling();
        samplingScheme.show();
        samplingScheme.getWindow().setSize(600, 600);
        samplingScheme.getCanvas().fitToWindow();

        ImagePlus weightScheme = sph2d.showWeights();
        weightScheme.show();
        weightScheme.getWindow().setSize(600, 600);
        weightScheme.getCanvas().fitToWindow();

        /*
        main
         */
        long t1, t2;
        t1 = System.currentTimeMillis();


        Masker2D.loadTemplate(inimg_xy, 0, sph2d.getOuterRadius(), iDiff);
        int totalLocs = inimg_xy.length * inimg_xy[0].length;

        Masker2D ms_jobs[] = new Masker2D[CPU_NR];
        for (int i = 0; i < ms_jobs.length; i++) {
            ms_jobs[i] = new Masker2D(i*totalLocs/CPU_NR,  (i+1)*totalLocs/CPU_NR);
            ms_jobs[i].start();
        }
        for (int i = 0; i < ms_jobs.length; i++) {
            try {
                ms_jobs[i].join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }

        Masker2D.formRemainingOutputs();


        Profiler2D.loadTemplate(sph2d, Masker2D.i2xy, inimg_xy);
        int totalProfileComponents = sph2d.getProfileLength();

        Profiler2D pf_jobs[] = new Profiler2D[CPU_NR];
        for (int i = 0; i < pf_jobs.length; i++) {
            pf_jobs[i] = new Profiler2D(i*totalProfileComponents/CPU_NR,  (i+1)*totalProfileComponents/CPU_NR);
            pf_jobs[i].start();
        }
        for (int i = 0; i < pf_jobs.length; i++) {
            try {
                pf_jobs[i].join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }


        PeakExtractor2D.loadTemplate(sph2d, Masker2D.i2xy, Profiler2D.prof2, inimg_xy, Masker2D.xy2i);
        int totalPeakExtrComponents = Profiler2D.prof2.length; // number of profiles == number of locations i2xy.length

        PeakExtractor2D pe_jobs[] = new PeakExtractor2D[CPU_NR];
        for (int i = 0; i < pe_jobs.length; i++) {
            pe_jobs[i] = new PeakExtractor2D(i*totalPeakExtrComponents/CPU_NR, (i+1)*totalPeakExtrComponents/CPU_NR);
            pe_jobs[i].start();
        }
        for (int i = 0; i < pe_jobs.length; i++) {
            try {
                pe_jobs[i].join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }


        PeakAnalyzer2D.loadTemplate(Masker2D.i2xy, Masker2D.xy2i, PeakExtractor2D.peaks_xy, M, minCos, scatterDist); // initialize peak analyzer parameters
        int totalPeakAnalyzeComponents = Masker2D.i2xy.length; // number of locations

        PeakAnalyzer2D pa_jobs[] = new PeakAnalyzer2D[CPU_NR];
        for (int i = 0; i < pa_jobs.length; i++) {
            pa_jobs[i] = new PeakAnalyzer2D(i*totalPeakAnalyzeComponents/CPU_NR, (i+1)*totalPeakAnalyzeComponents/CPU_NR);
            pa_jobs[i].start();
        }
        for (int i = 0; i < pa_jobs.length; i++) {
            try {
                pa_jobs[i].join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }

        t2 = System.currentTimeMillis();
        IJ.log("done. "+((t2-t1)/1000f)+"sec.");

        /*
			mouse listeners after processing
		 */
        cnv.addMouseListener(this);
        cnv.addMouseMotionListener(this);

        // add foreground mask on top
		//IJ.log("window title: "+cnv.getImage().getWindow().getTitle());

		// overlay foreground
		ImagePlus foregroundMask = new ImagePlus("foreground", Masker2D.getMask());

		//foregroundMask.show();
		//String foregroundMaskName = foregroundMask.getWindow().getTitle();

		IJ.selectWindow(cnv.getImage().getWindow().getTitle());
		IJ.setTool("hand");
		//IJ.run("Add Image...", "image="+"foreground"+" x=0 y=0 opacity=25");
		//foregroundMask.close();


    }


    public void mouseClicked(MouseEvent e)
    {

        mouseMoved(e);
        IJ.log("TRAIN! POSSIBLE");

    }

    public void mouseMoved(MouseEvent e)
    {
        //mouseClicked(e);

        int clickX = cnv.offScreenX(e.getX());
        int clickY = cnv.offScreenY(e.getY());

        Overlay ov = PeakAnalyzer2D.getDelineation(clickX, clickY);
		ImageRoi fgroi = new ImageRoi(0, 0, Masker2D.getMask());// !!! not very efficient to be done each time
		fgroi.setOpacity(0.15);
		ov.add(fgroi);

        cnv.setOverlay(ov);

    }

    public void mousePressed(MouseEvent e)  {}
    public void mouseReleased(MouseEvent e) {}
    public void mouseEntered(MouseEvent e)  {}
    public void mouseExited(MouseEvent e)   {}
    public void mouseDragged(MouseEvent e)  {}

}
