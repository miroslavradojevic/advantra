package detection2d;

import ij.IJ;
import ij.plugin.PlugIn;
import tracing2d.BayesianTracerMulti;

/**
 * Created by miroslav on 1-12-14.
 * extract critpoint score for foreground selection of points - threaded
 */
public class NS_LandmarkDet implements PlugIn {

    int R               = 4;     // will cover radiuses R * BayesianTracerMulti.sstep
    int percentile = 95;    // how many to keep as the foreground/background, neighbourhood will be defined with radius and
    int nbhood = (int) Math.ceil(2*R* BayesianTracerMulti.sstep[BayesianTracerMulti.sstep.length-1]);
    int CPU_NR = Runtime.getRuntime().availableProcessors() + 1;


    public void run(String s) {

                /*
            extract the foreground -> Masker2D.xy2i, Masker2D.i2xy
         */
        IJ.log("extracting foreground...");
//        Masker2D.loadTemplate(
//                likelihood_xy, nbhood, nbhood, percentile); //image, margin, check, percentile
//        int totalLocs = likelihood_xy.length * likelihood_xy[0].length;
//        Masker2D ms_jobs[] = new Masker2D[CPU_NR];
//        for (int i = 0; i < ms_jobs.length; i++) {
//            ms_jobs[i] = new Masker2D(i*totalLocs/CPU_NR,  (i+1)*totalLocs/CPU_NR);
//            ms_jobs[i].start();
//        }
//        for (int i = 0; i < ms_jobs.length; i++) {
//            try {
//                ms_jobs[i].join();
//            } catch (InterruptedException e) {
//                e.printStackTrace();
//            }
//        }
//        Masker2D.defineThreshold();
//        Masker2D.formRemainingOutputs();
        IJ.log(" done. ");

//        ImagePlus mask = new ImagePlus("mask", Masker2D.getMask());
//        mask.show();//      IJ.saveAs(mask, "Tiff", midresults_dir+"mask_"+D[didx]+".tif");


    }
}
