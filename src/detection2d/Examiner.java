package detection2d;

import graphcut.Graph_Cut;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.*;
import ij.io.FileSaver;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.io.File;
import java.util.Random;

/**
 * Created by miroslav on 1/6/14.
 * Demo usage of PeakAnalyzer2D: uses also Masker2D, Profiler2D, and PeakExtractor2D shows delineated
 * branching structure extracted on all foreground locations as a joint Overlay on the image
 */
public class Examiner implements PlugInFilter, MouseListener, MouseMotionListener {

    float[][] 	inimg_xy;               // store input image as an array
	String 		image_dir;
	String		image_name;

    /*
        parameters
     */
    float       s = 1.3f;
    float       D;
    float       minCos;
    int         M;
    float       mu_ON;
    float       mu_OFF;
    float       sig_ON;
    float       sig_OFF;
    int         L;
    float       sampling_crosswise;
    String      mode;

    int         CPU_NR;

    /*
    interface elements - all the windows that pop up as you click/move with mouse
     */
    ImagePlus       pfl_im  = new ImagePlus();    // used with live inspections (viz patches)
    ImagePlus       pfl_im1 = new ImagePlus();    // used with live inspections (plot)
    ImagePlus       pfl_im2 = new ImagePlus();
    ImagePlus       pfl_im3 = new ImagePlus();
    ImagePlus       pfl_im4 = new ImagePlus();

    ImageStack      pfl_is  = null;
    ImageStack      pfl_is1 = null;
    ImageStack      pfl_is2 = null;
    ImageStack      pfl_is3 = null;
    ImageStack      pfl_is4 = null;

	// output
	static int[][] 		    thr_number; // keep the maximum here always
	static float[][][]		thr_fit;
	static float[][][]		circ_stat; // (R_, kurtosis, circ std, sample circ dispersion)

    ImageCanvas cnv;

    public int setup(String s, ImagePlus imagePlus) {

        if(imagePlus==null) return DONE;

        inimg_xy = new float[imagePlus.getWidth()][imagePlus.getHeight()]; // x~column, y~row

		image_dir = imagePlus.getOriginalFileInfo().directory; //  + File.separator  + image_name
		image_name = imagePlus.getShortTitle();

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
        this.D 					        	= (float)   Prefs.get("critpoint.detection2d.d", 4);
        this.M 					        	= (int)     Prefs.get("critpoint.detection2d.m", 2);
        this.minCos 					    = (float)   Prefs.get("critpoint.detection2d.min_cos", -0.5f);
        this.mu_ON                          = (float)   Prefs.get("critpoint.detection2d.mu_on", 0f);
        this.mu_OFF 					    = (float)   Prefs.get("critpoint.detection2d.mu_off", 0.5f);
		this.sig_ON						    = (float) 	Prefs.get("critpoint.detection2d.sig_on", 0.1f);
		this.sig_OFF						= (float) 	Prefs.get("critpoint.detection2d.sig_off", 0.1f);
		this.L						        = (int) 	Prefs.get("critpoint.detection2d.l", 4);
		this.sampling_crosswise			    = (float) 	Prefs.get("critpoint.detection2d.sampling_crosswise", 0.3f);
		this.mode			                = (String) 	Prefs.get("critpoint.detection2d.mode", "NCC");

        GenericDialog gd = new GenericDialog("EXAMINER");
        gd.addNumericField("D ", 	                D, 			        0,  10, "neuron diameter[pix]");
        gd.addNumericField("M ", 	                M, 			        0,  10, "");
        gd.addNumericField("min cos ", 	            minCos, 	        1,  10, "");
        gd.addNumericField("mu_ON ", 	            mu_ON, 			    1,  10, "");
        gd.addNumericField("mu_OFF ",               mu_OFF, 	        1,  10, "");
		gd.addNumericField("sig_ON ",               sig_ON, 		    1,  10, "");
		gd.addNumericField("sig_OFF ",              sig_OFF, 		    1,  10, "");
		gd.addNumericField("L ",                    L, 		            0,  10, "");
		gd.addNumericField("sampling_crosswise ",   sampling_crosswise, 1,  10, "");
        gd.addStringField("mode ",                  mode,                   10);

        gd.showDialog();
        if (gd.wasCanceled()) return DONE;

        D       	    = (float) gd.getNextNumber();
        Prefs.set("critpoint.detection2d.d", 	    	D);

        M       	    = (int) gd.getNextNumber();
        Prefs.set("critpoint.detection2d.m", 	    	M);

        minCos       	= (float) gd.getNextNumber();
        Prefs.set("critpoint.detection2d.min_cos", 	minCos);

        mu_ON       	= (float) gd.getNextNumber();
        Prefs.set("critpoint.detection2d.mu_on", 	    mu_ON);

        mu_OFF          = (float) gd.getNextNumber();
        Prefs.set("critpoint.detection2d.mu_off", 	    mu_OFF);

		sig_ON		= (float) gd.getNextNumber();
		Prefs.set("critpoint.detection2d.sig_on", 	    sig_ON);

        sig_OFF		= (float) gd.getNextNumber();
        Prefs.set("critpoint.detection2d.sig_off", 	sig_OFF);

        L		    = (int) gd.getNextNumber();
        Prefs.set("critpoint.detection2d.l", 	        L);

        sampling_crosswise		= (float) gd.getNextNumber();
        Prefs.set("critpoint.detection2d.sampling_crosswise", 	sampling_crosswise);

        mode		            = gd.getNextString();
        Prefs.set("critpoint.detection2d.mode", 	    mode);

        CPU_NR = Runtime.getRuntime().availableProcessors() + 1;

        cnv = imagePlus.getCanvas();

        return DOES_8G+DOES_32+NO_CHANGES;

    }

    public void run(ImageProcessor imageProcessor) {

		// todo add D loop here and store the outputs into static outputs
		for (int d_scale=0; d_scale<3; d_scale++) {
			float DD = D + d_scale * 1;

            if (true) {

                Sphere2D sph2d1 = new Sphere2D(DD, s);

                ImagePlus samplingScheme =  sph2d1.showSampling();
                samplingScheme.show();
                samplingScheme.getWindow().setSize(200, 200);
                samplingScheme.getCanvas().fitToWindow();

                ImagePlus weightScheme = sph2d1.showWeights();
                weightScheme.show();
                weightScheme.getWindow().setSize(200, 200);
                weightScheme.getCanvas().fitToWindow();

            }
		}

        Sphere2D sph2d = new Sphere2D(D, s);

        IJ.log("PROCESSING D="+D+", s="+s);

        /*
        main
         */
        long t1, t2;
        t1 = System.currentTimeMillis();
        //************************************************
        IJ.log("extracting mask...");
        Masker2D.loadTemplate(inimg_xy, sph2d.getOuterRadius(), 1.5f*sph2d.getOuterRadius()); //margin = sph2d.getOuterRadius()
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
		Masker2D.defineThreshold();
        Masker2D.formRemainingOutputs();
		new ImagePlus("MASK", Masker2D.getMask()).show();
        //************************************************
        IJ.log("calculating profiles...");
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
        //************************************************
        IJ.log("extracting peaks...");
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
        PeakExtractor2D.getCircStat().show();
        //************************************************
        IJ.log("analyzing peaks + extracting features & descriptors...");
        Delineator2D.loadTemplate(Masker2D.i2xy, Masker2D.xy2i, PeakExtractor2D.peaks_i, PeakExtractor2D.peaks_w, PeakExtractor2D.peaks_lhood, inimg_xy, //Masker2D.back_xy,
                M,
                minCos,
                D,
                mu_ON,
                mu_OFF,
                sig_ON,
                sig_OFF,
                L,
                sampling_crosswise,
                mode
                );
        int totalPeakAnalyzeComponents = Masker2D.i2xy.length; // number of locations
        Delineator2D pa_jobs[] = new Delineator2D[CPU_NR];
        for (int i = 0; i < pa_jobs.length; i++) {
            pa_jobs[i] = new Delineator2D(i*totalPeakAnalyzeComponents/CPU_NR, (i+1)*totalPeakAnalyzeComponents/CPU_NR);
            pa_jobs[i].start();
        }
        for (int i = 0; i < pa_jobs.length; i++) {
            try {
                pa_jobs[i].join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
        //************************************************
        /*
        IJ.log("simple detector...");
        SimpleDetector2D.loadTemplate(inimg_xy.length, inimg_xy[0].length, Masker2D.i2xy, Masker2D.xy2i, Delineator2D.delin2, Delineator2D.lhood2);
        int totalSimpleDetectComponents = Masker2D.i2xy.length;
        SimpleDetector2D sd_jobs[] = new SimpleDetector2D[CPU_NR];
        for (int i=0; i<sd_jobs.length; i++) {
            sd_jobs[i] = new SimpleDetector2D(i*totalSimpleDetectComponents/CPU_NR, (i+1)*totalSimpleDetectComponents/CPU_NR);
            sd_jobs[i].start();
        }
        for (int i=0; i<sd_jobs.length; i++) {
            try {
                sd_jobs[i].join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
        */

        ImageStack is_lhoods = Delineator2D.exportLikelihoods(new int[]{1});
        new ImagePlus("ENDPOINTS", is_lhoods).show();
        is_lhoods = Delineator2D.exportLikelihoods(new int[]{3, 4});
        new ImagePlus("JUNCTIONS", is_lhoods).show();


        t2 = System.currentTimeMillis();
		IJ.log("done. " + ((t2 - t1) / 1000f) + "sec.");

		//        Overlay overlay_with_detections = SimpleDetector2D.drawDetections();
//        ImagePlus final_det = cnv.getImage().duplicate();
//        final_det.show();
//        final_det.setOverlay(overlay_with_detections); // add the mto the original image

//        ImagePlus end_lhood = SimpleDetector2D.showEndLhood();
//        end_lhood.show();
//        ImagePlus jun_lhood = SimpleDetector2D.showJunctionLhood();
//        jun_lhood.show();
//        ImagePlus model_grad = SimpleDetector2D.showModelGradient();
//        model_grad.show();

        // use Graphcut to get segment out stable regions // TODO put it in SimpleDetector2D
//        Graph_Cut gc = new Graph_Cut();
//        ImagePlus out = gc.processSingleChannelImage(jun_lhood, null, 500, 500, 500);
//        out.setTitle("GCsegmentation");
//        out.show();

        /*
        if (false) {

            IJ.log("export features ");
            t1 = System.currentTimeMillis();

            String export_path = image_dir+image_name+".feat";
            Delineator2D.exportFeats(export_path);
            IJ.log("done exporting to: \t" + export_path + "\n");

            export_path = image_dir+image_name+".desc";
            Delineator2D.exportDescripts(export_path);
            IJ.log("done exporting to: \t" + export_path + "\n");

            export_path = image_dir + image_name + ".frame";
            Delineator2D.exportFrames(export_path);
            IJ.log("done exporting to: \t" + export_path + "\n");

            export_path = image_dir + image_name + ".i2xy";
            Masker2D.exportI2xyCsv(export_path);
            IJ.log("done exporting to: \t" + export_path + "\n");

            t2 = System.currentTimeMillis();
            IJ.log("done. " + ((t2 - t1) / 1000f) + "sec.");

        }
        */

		//************************************************

        /*
			mouse listeners after processing
		 */
        cnv.addMouseListener(this);
        cnv.addMouseMotionListener(this);

		IJ.selectWindow(cnv.getImage().getWindow().getTitle());
		IJ.setTool("hand");
        IJ.selectWindow(cnv.getImage().getTitle());

    }

    public void mouseClicked(MouseEvent e)
    {

//        mouseMoved(e);

        int clickX = cnv.offScreenX(e.getX());
        int clickY = cnv.offScreenY(e.getY());

        /*
            output Stack with local patches & update window
         */
        /*
        pfl_is = PeakAnalyzer2D.getDelineationPatches(clickX, clickY);
        pfl_im.setStack("patches", pfl_is);
        pfl_im.show();
        pfl_im.getWindow().setSize(250, 500);
        pfl_im.getCanvas().fitToWindow();

        pfl_is1 = PeakAnalyzer2D.plotDelineationProfiles(clickX, clickY);
        pfl_im1.setStack("cross-profiles", pfl_is1);
        pfl_im1.show();
        */


        IJ.setTool("hand");

    }

    public void mouseMoved(MouseEvent e)
    {
        //mouseClicked(e);

        int clickX = cnv.offScreenX(e.getX());
        int clickY = cnv.offScreenY(e.getY());

        /*
            output Overlay & update canvas with the original
         */
        Overlay ov_to_add = Delineator2D.getDelineationOverlay(clickX, clickY);

        // add the whole image on top as overlay each time mouse is moved
        //ImageRoi fgroi = new ImageRoi(0, 0, Masker2D.getMask());    // !!! not very efficient to be done each time
        //ImageRoi simple_det_roi = new ImageRoi(0, 0, new ByteProcessor(inimg_xy.length, inimg_xy[0].length, SimpleDetsector2D.score2));

        //simple_det_roi.setOpacity(0.1);     // add foreground to know what was removed always
        //ov_to_add.add(simple_det_roi);
        //simple_det_roi.setOpacity(0.1);     // add foreground to know what was removed always
		//ov.add(simple_det_roi);

        cnv.setOverlay(ov_to_add);

		// print extracted features
		Delineator2D.print(clickX, clickY);

        /*
         output stack with local profile
        */
//        pfl_is3 = Profiler2D.getProfile(clickX, clickY);
//        pfl_im3.setStack("local_profile", pfl_is3);
//        pfl_im3.show();
		/*
			output stack: local profile with peaks detected
		 */
        pfl_is4 = PeakExtractor2D.getProfileWithPeaks(clickX, clickY);
        pfl_im4.setStack("local_profile_with_peaks", pfl_is4);
        pfl_im4.show();

		pfl_is2 = Delineator2D.plotDelineationFeatures(clickX, clickY);
		pfl_im2.setStack("features", pfl_is2);
		pfl_im2.show();

		IJ.setTool("hand");

	}

    public void mousePressed(MouseEvent e)  {}
    public void mouseReleased(MouseEvent e) {}
    public void mouseEntered(MouseEvent e)  {}
    public void mouseExited(MouseEvent e)   {}
    public void mouseDragged(MouseEvent e)  {}

}