package detection2d;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 12/17/13
 * Time: 11:28 AM
 * Tests Profiler2D, extract the mask using Masker2D first and extract the profiles
 */
public class Profiler2DDemo implements PlugInFilter, MouseListener, MouseMotionListener {

    float[][] 	inimg_xy; // store input image as an array

    /*
    parameters
     */
    float       iDiff, D, s=1.5f; // maskNhoodRadius
    int         CPU_NR;

	/*
	interface things
	 */
	ImagePlus       pfl_im = new ImagePlus();
	ImageStack      pfl_ip = null;
	ImageWindow     pfl_iw;
	ImageCanvas 	cnv;

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
        //this.maskNhoodRadius             	= (float)   Prefs.get("advantra.critpoint.mask.nhoodRadius", 5);
        this.iDiff 					    	= (float)   Prefs.get("advantra.critpoint.mask.iDiff", 5);
        this.D 					        	= (float)   Prefs.get("advantra.critpoint.profile.d", 4);

        GenericDialog gd = new GenericDialog("PROFILER2DDEMO");
        //gd.addNumericField("radius ", 	maskNhoodRadius, 	0, 10, "spatial neighbourhood");
        gd.addNumericField("iDiff ", 	iDiff, 			0, 10, "intensity margin");
        gd.addNumericField("D ", 	D, 			0, 10, "neuron diameter[pix]");

        gd.showDialog();
        if (gd.wasCanceled()) return DONE;

        //maskNhoodRadius      = (float) gd.getNextNumber();
        //Prefs.set("advantra.critpoint.mask.nhoodRadius",	maskNhoodRadius);

        iDiff       	= (float) gd.getNextNumber();
        Prefs.set("advantra.critpoint.mask.iDiff", 	    	iDiff);

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
        sph2d.showWeights().show();

        long t1, t2;

        /*
        main
         */
        t1 = System.currentTimeMillis();

        Masker2D.loadTemplate(inimg_xy, 0, sph2d.getOuterRadius());
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

        Profiler2D.loadTemplate(sph2d, Masker2D.i2xy, Masker2D.xy2i, inimg_xy);
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

        t2 = System.currentTimeMillis();

        IJ.log("done. "+((t2-t1)/1000f)+"sec.");

        ImagePlus outmask = new ImagePlus("mask", Masker2D.getMask());
        outmask.show();

		IJ.log("click on the original image...");

		cnv.addMouseListener(this);
		cnv.addMouseMotionListener(this);
    }

	public void mouseClicked(MouseEvent e) {

		int clickX = cnv.offScreenX(e.getX());
		int clickY = cnv.offScreenY(e.getY());

		pfl_ip = Profiler2D.getProfile(clickX, clickY);
        pfl_im.setStack(pfl_ip);

		if (pfl_iw==null) pfl_iw = new ImageWindow(pfl_im);
		else pfl_iw.setImage(pfl_im);

		//pfl_iw.setSize(600, 300);
		//pfl_iw.getCanvas().fitToWindow();

	}

	public void mouseMoved(MouseEvent e) {
		mouseClicked(e);
	}

	public void mousePressed(MouseEvent e) {}
	public void mouseReleased(MouseEvent e) {}
	public void mouseEntered(MouseEvent e) {}
	public void mouseExited(MouseEvent e) {}
	public void mouseDragged(MouseEvent e) {}

}
