package profile;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;

import java.awt.*;
import java.io.File;
import java.util.ArrayList;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 7/11/13
 * Time: 1:46 PM
 */
public class MS2DDemo implements PlugInFilter {

    ImagePlus 	inimg;
    ImagePlus   inmask;
    String      inmaskPath;
    /*
    parameters
     */
    int H, MaxIter, M, CPU_NR;
    double eps, d;
    boolean useMask;


    public int setup(String s, ImagePlus imagePlus) {

        if(imagePlus==null) return DONE;
        inimg = Tools.convertToFloatImage(imagePlus);
        inimg.setTitle("inimg");

        if (imagePlus==null || imagePlus.getOriginalFileInfo()==null) { //  || imagePlus.getOriginalFileInfo().directory==""
            inmaskPath = null;
        }
        else {
            inmaskPath = Tools.removeExtension(imagePlus.getOriginalFileInfo().directory+imagePlus.getOriginalFileInfo().fileName)+".mask";
        }
        return DOES_8G+DOES_32+NO_CHANGES;
    }

    public void run(ImageProcessor imageProcessor) {

        H       = (int) Prefs.get("advantra.critpoint.ms2d.H", 4);
        MaxIter = (int) Prefs.get("advantra.critpoint.ms2d.MaxIter", 200);
        eps     =       Prefs.get("advantra.critpoint.ms2d.eps", 0.0001);
        d       =       Prefs.get("advantra.critpoint.ms2d.d", 0.5);
        M       = (int) Prefs.get("advantra.critpoint.ms2d.M", 5);
        useMask =       Prefs.get("advantra.critpoint.ms2d.useMask", true);
        CPU_NR  = (int) Prefs.get("advantra.critpoint.CPU_NR", 4);

        GenericDialog gd = new GenericDialog("MEAN SHIFT");
        gd.addNumericField("H ",        H,      0, 10, "spatial neighbourhood");
        gd.addNumericField("Max Iter",  MaxIter,0, 10, "max #iterations");
        gd.addNumericField("Eps",       eps,    6, 10, "convergence threshold");
        gd.addMessage("---");
        gd.addNumericField("d",         d,      1, 10, "min intra-cluster distance");
        gd.addNumericField("M",         M,      0, 10, "min # points in cluster");
        gd.addMessage("suggested mask: "+inmaskPath);
        gd.addCheckbox("try suggested mask", true);
        gd.addNumericField("CPU_NR ",   CPU_NR, 1, 10, "spatial neighbourhood");

        gd.showDialog();
        if (gd.wasCanceled()) return;

        H       =  (int) gd.getNextNumber();
        Prefs.set("advantra.critpoint.ms2d.H", 	H);
        MaxIter =  (int) gd.getNextNumber();
        Prefs.set("advantra.critpoint.ms2d.MaxIter", MaxIter);
        eps     =  gd.getNextNumber();
        Prefs.set("advantra.critpoint.ms2d.eps", eps);
        d       =  gd.getNextNumber();
        Prefs.set("advantra.critpoint.ms2d.d", d);
        M       =  (int) gd.getNextNumber();
        Prefs.set("advantra.critpoint.ms2d.M", M);
        useMask = gd.getNextBoolean();
        Prefs.set("advantra.critpoint.ms2d.useMask", useMask);
        CPU_NR       =  (int) gd.getNextNumber();
        Prefs.set("advantra.critpoint.CPU_NR", 	CPU_NR);

        if (useMask && inmaskPath!=null && new File(inmaskPath).exists()) {
            IJ.log("suggested mask found...");
            inmask = new ImagePlus(inmaskPath);
            inmask.setTitle("inmask");
        }
        else { // take all pixels
            IJ.log("suggested mask not found...");
            byte[] array = new byte[inimg.getWidth() * inimg.getHeight()];
            for (int a = 0; a < array.length; a++) array[a] = (byte) 255;
            inmask = new ImagePlus("inmask", new ByteProcessor(inimg.getWidth(), inimg.getHeight(), array));
        }



        IJ.log("mean-shifting...");
		long t1 = System.currentTimeMillis();
        MS2D.loadTemplate(inimg.getProcessor(), (ByteProcessor) inmask.getProcessor(), H, MaxIter, eps);
        int totalProfiles = MS2D.toProcess;
        MS2D ms_jobs[] = new MS2D[CPU_NR];
        for (int i = 0; i < ms_jobs.length; i++) {
            ms_jobs[i] = new MS2D(i*totalProfiles/CPU_NR,  (i+1)*totalProfiles/CPU_NR);
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




        IJ.log("clustering...");
		t1=System.currentTimeMillis();
		ArrayList<ArrayList<double[]>> clusters = MS2D.extractClusters(d);
		t2=System.currentTimeMillis();
		IJ.log("done " + ((t2 - t1) / 1000f) + " s.");



        IJ.log("cluster centroids...");
        t1=System.currentTimeMillis();
        ArrayList<double[]> centroids = MS2D.extractClusterCentroids(clusters, M);
        t2=System.currentTimeMillis();
        IJ.log("done " + ((t2 - t1) / 1000f) + " sec.");


        Overlay ov = new Overlay();
        for (int i1=0; i1<centroids.size(); i1++) {

            //IJ.log("centroid "+i1+" : size -> "+centroids.get(i1)[2]);

            double atX = centroids.get(i1)[0];
            double atY = centroids.get(i1)[1];

            //double rds = Math.sqrt(centroids.get(i1)[2]/Math.PI);
            double rds = centroids.get(i1)[3]/40f;

            OvalRoi oval = new OvalRoi(atX+0.5-rds, atY+0.5-rds, 2*rds, 2*rds);
            oval.setStrokeColor(Color.RED);
            ov.add(oval);

        }


        inimg.show();
        inimg.setOverlay(ov);

//        inmask.show();
//        IJ.selectWindow("inimg");
//        IJ.run("Add Image...", "image=inmask x="+0+" y="+0+" opacity=30");
//        inmask.close();

        IJ.selectWindow("inimg");
        IJ.setTool("hand");
        inimg.getCanvas().zoomIn(0, 0);
        inimg.getCanvas().zoomIn(0, 0);
//        inimg.getCanvas().zoomIn(0, 0);
//        inimg.getCanvas().zoomIn(0, 0);


    }
}
