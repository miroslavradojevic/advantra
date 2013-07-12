package profile;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 7/12/13
 * Time: 11:45 AM
 */
public class MaskerDemo implements PlugInFilter {

    ImagePlus 	inimg;

    /*
    parameters
     */
    int N, CPU_NR;
    float th;


    public int setup(String s, ImagePlus imagePlus) {
        if(imagePlus==null) return DONE;
        inimg = Tools.convertToFloatImage(imagePlus);
        inimg.setTitle("inimg");
        return DOES_8G+DOES_32+NO_CHANGES;
    }

    public void run(ImageProcessor imageProcessor) {

        /*
        Generic Dialog
         */
        N       = (int) Prefs.get("advantra.critpoint.mask.N", 9);
        th      = (float) Prefs.get("advantra.critpoint.mask.th", 10);
        CPU_NR  = (int) Prefs.get("advantra.critpoint.CPU_NR", 4);

        GenericDialog gd = new GenericDialog("MASK CREATOR");
        gd.addNumericField("N ", N, 0, 10, "spatial neighbourhood 2N+1x2N+1 to get median");
        gd.addNumericField("th", th, 1, 10, "higher than background");
        gd.addNumericField("CPU_NR ",   CPU_NR, 1, 10, "spatial neighbourhood");

        gd.showDialog();
        if (gd.wasCanceled()) return;

        N       =  (int) gd.getNextNumber();
        Prefs.set("advantra.critpoint.mask.N", 	N);
        th      = (float) gd.getNextNumber();
        Prefs.set("advantra.critpoint.mask.th", th);
        CPU_NR       =  (int) gd.getNextNumber();
        Prefs.set("advantra.critpoint.CPU_NR", 	CPU_NR);

        /*
        main
         */

        IJ.log("extracting background...");
        long t1 = System.currentTimeMillis();
        Masker.loadTemplate(inimg.getProcessor(), N, th);
        int totalProfiles = Masker.image_height*Masker.image_width;
        Masker ms_jobs[] = new Masker[CPU_NR];
        for (int i = 0; i < ms_jobs.length; i++) {
            ms_jobs[i] = new Masker(i*totalProfiles/CPU_NR,  (i+1)*totalProfiles/CPU_NR);
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
        IJ.log("done "+((t2-t1)/1000f)+" sec.");

        inimg.show();
        ImagePlus inmask = new ImagePlus("inmask", Masker.maskip);

        // overlay mask
        inmask.show();
        IJ.selectWindow("inimg");
        IJ.run("Add Image...", "image=inmask x="+0+" y="+0+" opacity=40");
        //inmask.close();

        // show extracted backgr
        new ImagePlus("backgr", Masker.back).show();

        IJ.selectWindow("inimg");
        IJ.setTool("hand");
        inimg.getCanvas().zoomIn(0, 0);
        inimg.getCanvas().zoomIn(0, 0);

    }
}
