package profile;

import conn.Find_Connected_Regions;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;

import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 7/12/13
 * Time: 11:45 AM
 */
public class MaskerDemo implements PlugInFilter {

    ImagePlus 	inimg;
    String      inimgTitle;
    ImageCanvas incanvas;
    ImagePlus   inmask;

    /*
    parameters
     */
    float       nhoodRadius;
    int         CPU_NR;
//    float   th;
    String      bgExtractionMode;
    boolean     bgLocal;

    long t1, t2;

    public int setup(String s, ImagePlus imagePlus) {
        if(imagePlus==null) return DONE;
        inimg = Tools.convertToFloatImage(imagePlus);
        inimgTitle = imagePlus.getTitle();
        incanvas = imagePlus.getCanvas();
        return DOES_8G+DOES_32+NO_CHANGES;
    }

    public void run(ImageProcessor imageProcessor) {

        /*
        Generic Dialog
         */
        nhoodRadius             = (float)   Prefs.get("advantra.critpoint.mask.nhoodRadius", 5);
        bgExtractionMode        =           Prefs.get("advantra.critpoint.mask.bgExtractionMode", "MEAN");
        bgLocal                 =           Prefs.get("advantra.critpoint.mask.bgLocal", true);


        GenericDialog gd = new GenericDialog("MASK EXTRACTOR");
        gd.addNumericField("radius ", nhoodRadius, 0, 10, "spatial neighbourhood");
        gd.addChoice("background extraction", new String[]{"MEAN", "MEDIAN"}, bgExtractionMode);
        gd.addCheckbox("local neighbourhood", true);

        gd.showDialog();
        if (gd.wasCanceled()) return;

        nhoodRadius       = (float) gd.getNextNumber();
        Prefs.set("advantra.critpoint.mask.nhoodRadius", 	    nhoodRadius);
        bgExtractionMode      = gd.getNextChoice();
        Prefs.set("advantra.critpoint.mask.bgExtractionMode",   bgExtractionMode);
        bgLocal       =  gd.getNextBoolean();
        Prefs.set("advantra.critpoint.mask.bgLocal", 	        bgLocal);

		CPU_NR = Runtime.getRuntime().availableProcessors();

        /*
        main
         */

        IJ.log("extracting background...    "+bgExtractionMode+", local? "+bgLocal);

        t1 = System.currentTimeMillis();

        if (bgExtractionMode=="MEAN") {
            Masker.loadTemplate(inimg.getProcessor(), nhoodRadius, 0, bgLocal);
        }
        if (bgExtractionMode=="MEDIAN") {
            Masker.loadTemplate(inimg.getProcessor(), nhoodRadius, 1, bgLocal);
        }

        int totalProfiles = inimg.getHeight()*inimg.getWidth();

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
        t2 = System.currentTimeMillis();
        IJ.log("done "+((t2-t1)/1000f)+" sec.");

        inmask = new ImagePlus("inmask", Masker.maskip);
        inmask.show();

        // show extracted backgr
        //new ImagePlus("backgr", Masker.back).show();

        // check amount of extracted locations
        int nr = 0;
        for (int i=0; i<totalProfiles; i++) if (Masker.maskip.getf(i)==255) nr++;

        GenericDialog gd1 = new GenericDialog("PRUNE REGIONS?");
        gd1.addMessage("extracted "+nr+" out of "+totalProfiles+"  ("+IJ.d2s(nr*100f/totalProfiles, 2)+"%)");
        gd1.addMessage("PRUNE REGIONS?");

        gd1.showDialog();
        if (!gd1.wasCanceled()) {

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

			IJ.log("pruning...");

			for (int i=0; i<regs.size(); i++) {

				if (regs.get(i).size()>Masker.VERY_SMALL_REGION_SIZE) {

					float regAvg = 0;
					float bkgAvg = 0;
					for (int w=0; w<regs.get(i).size(); w++) {
						int takeValX = regs.get(i).get(w)[1];
						int takeValY = regs.get(i).get(w)[0];
						regAvg += inimg.getProcessor().getf(takeValX, takeValY);
						bkgAvg += Masker.back.getf(takeValX, takeValY);
					}
					regAvg /= regs.get(i).size();
					bkgAvg /= regs.get(i).size();

					if (regAvg < bkgAvg + Masker.VISIBLE_INTENSITY_DIFF ) { // Masker.back.getf((int)centroidX, (int)centroidY)
						for (int q=0; q<regs.get(i).size(); q++) {
							int removeX = regs.get(i).get(q)[1];
							int removeY = regs.get(i).get(q)[0];
							inmask.getProcessor().set(removeX, removeY, 0);
						}
					}
//                    else {
//                        //IJ.log("KEPT: region average: "+regAvg+", size "+regs.get(i).size());
//                    }
				}
				else {  // prune small regions
					for (int q=0; q<regs.get(i).size(); q++) {
						int removeX = regs.get(i).get(q)[1];
						int removeY = regs.get(i).get(q)[0];
						inmask.getProcessor().set(removeX, removeY, 0);
					}
				}

			}

        }

		IJ.log("done.");

		inmask.updateAndDraw();
		IJ.selectWindow(inimgTitle);
		IJ.run("Add Image...", "image=inmask x="+0+" y="+0+" opacity=40");
		IJ.selectWindow(inimgTitle);
		IJ.setTool("hand");
		incanvas.zoomIn(0, 0);

        // check amount of extracted locations
        nr = 0;
        for (int i=0; i<totalProfiles; i++) if (Masker.maskip.getf(i)==255) nr++;
        IJ.log("extracted "+nr+" out of "+totalProfiles+"  ("+IJ.d2s(nr*100f/totalProfiles, 2)+"%)");

//        inmask.close();

    }
}
