package profile;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;

import java.io.File;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 7/11/13
 * Time: 1:46 PM
 */
public class MS2DDemo implements PlugInFilter {

    ImagePlus 	inimg;
    String 		inimgPath;
    ImagePlus   inmask;

    public int setup(String s, ImagePlus imagePlus) {
        if(imagePlus==null) {
            IJ.showMessage("needs image opened"); return DONE; }
        inimg = Tools.convertToFloatImage(imagePlus);
        inimg.setTitle("inimg");
        inimgPath = imagePlus.getOriginalFileInfo().directory+imagePlus.getOriginalFileInfo().fileName;
        return DOES_8G+DOES_32+NO_CHANGES;
    }

    public void run(ImageProcessor imageProcessor) {

        String inmaskPath = Tools.removeExtension(inimgPath)+".mask";

        if (new File(inmaskPath).exists() && true) {
            inmask = new ImagePlus(inmaskPath);
            inmask.setTitle("inmask");
        }
        else {
            byte[] array = new byte[inimg.getWidth() * inimg.getHeight()];
            for (int a = 0; a < array.length; a++) array[a] = (byte) 255;
            inmask = new ImagePlus("inmask", new ByteProcessor(inimg.getWidth(), inimg.getHeight(), array));
        }

        IJ.log("mean-shifting...");
		long t1 = System.currentTimeMillis();
        int CPU_NR = 8;
        MS2D.loadTemplate(inimg.getProcessor(), (ByteProcessor) inmask.getProcessor(), 5, 200, 0.0001);
        int totalProfiles = MS2D.toProcess;//offsets.size();
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
        IJ.log("done "+((t2-t1)/1000f)+" s.");

        IJ.log("clustering...");
		t1=System.currentTimeMillis();

		//MS2D.extractConvPoints(0.5, 50);
		MS2D.extractClusters(0.5, 50);

		t2=System.currentTimeMillis();
		IJ.log("done "+((t2-t1)/1000f)+" s.");


		double[][] Tnew = MS2D.T;
        IJ.log(Tnew.length+" points ");

        Overlay ov = new Overlay();
        for (int i1=0; i1<Tnew.length; i1++) {
            PointRoi pt = new PointRoi(Tnew[i1][1]+0.5, Tnew[i1][0]+0.5);
            ov.add(pt);
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
