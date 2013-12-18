package detection2d;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 12/17/13
 * Time: 4:10 PM
 * Demo usage of PeakExtractor2D: uses Masker2D, Profiler2D, and shows peaks extracted on all foreground locations as a joint Overlay on the image
 */
public class PeakExtractor2DDemo implements PlugInFilter {

    float[][] 	inimg_xy; // store input image as an array

    /*
    parameters
     */
    float       iDiff, D, s=1.5f; // maskNhoodRadius
    int         CPU_NR;

    /*
    interface things
     */
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

        /*
        main
         */
        long t1, t2;

        /*
        main
         */
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
        System.out.println("managed to initialize");
        int totalPeakExtrComponents = Profiler2D.prof2.length; // number of profiles == number of locations

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


        t2 = System.currentTimeMillis();

        IJ.log("done. "+((t2-t1)/1000f)+"sec.");

        // show the peaks (all peaks in one overlay)
        Overlay ov = new Overlay();
        float R_MAX = 2; // ~ corresponds highest weight
        float R = R_MAX;
        for (int l = 0; l<PeakExtractor2D.peaks2.length; l++) { // loop locations

            for (int t = 0; t<PeakExtractor2D.peaks2[l].length; t++) {  // loop threads

                int pk_x = PeakExtractor2D.peaks2[l][t][0];
                int pk_y = PeakExtractor2D.peaks2[l][t][1];

                if (pk_x != -1 && pk_y != -1) {
                    // R ~ size - weight of the peak
                    OvalRoi ovroi = new OvalRoi(pk_x+.5 - (R/2), pk_y+.5 - (R/2), R, R);
                    ov.add(ovroi);

                }

            }


        }

        cnv.getImage().setOverlay(ov);

    }
}
