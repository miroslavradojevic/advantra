package detection2d;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.*;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by miroslav on 1/6/14.
 * example terminal call:
 * java -Xmx4096m -jar ~/jarlib/ij.jar -ijpath ~/ImageJ/plugins/  ~/fuzzy_tests/n1_blue.tif -run "Detector2D"
 *
 */

// hint on zooming: imageplus.show(); imageplus.getWindow().setSize(200, 200); imageplus.getCanvas().fitToWindow();
public class Detector2D implements PlugInFilter, MouseListener, MouseMotionListener {

	String 		image_dir;
	String		image_name;

	float[][] 	inimg_xy;               // store input image as an array

    /*
        parameters
     */
	//
	int			critpoint_type;
	// model parameters
    float       s;
    String      Dlist;
    float       minCos = -.5f;   // hardcoded
    int         M;
	int         L;
	float       sampling_crosswise = .3f; // hardcoded
	//
	int 		Npoints = 40; // hardcoded, should not be too small but higher values can be computation costly
	float       ncc_high_mean;
    float       ncc_high_sigma;
    float       ncc_low_mean;
    float       ncc_low_sigma;
	float 		likelihood_high_mean;
	float 		likelihood_high_sigma;
	float 		likelihood_low_mean;
	float 		likelihood_low_sigma;
	float		output_sigma; // Fuzzy2D that will be used assumes out ranges from 0 to 1

	//
	float[]  	D;
	String 		image_gndtth;			// ground truth swc file with critical points to evaluate, same name as input image
	float		k = 1; // hardcoded, no need to tweak it each time, defines the stricktness of the threshold for out membership function
	float 		output_membership_th; // based on k and output_sigma
	int			min_size=3; // smallest connected region to be valid critical point
	String		output_dir_name; 		// parameter coded output folder name
	String 		output_log_name;
	PrintWriter logWriter = null;

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
	static float[][][] critpoint_det; // final off/on fuzzy memberships are stored and updated here

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
		this.critpoint_type					= (int)		Prefs.get("critpoint.detection2d.critpoint_type", 0); // 0-endpoint, 1-junction
        this.s								= (float)	Prefs.get("critpoint.detection2d.s", 1.2f);
		this.Dlist 					        =    		Prefs.get("critpoint.detection2d.d", "4");
		this.M 					        	= (int)     Prefs.get("critpoint.detection2d.m", 1);
		this.L						        = (int) 	Prefs.get("critpoint.detection2d.l", 8);

		this.ncc_high_mean                  = (float)   Prefs.get("critpoint.detection2d.ncc_high_mean", 1f);
        this.ncc_high_sigma 				= (float)   Prefs.get("critpoint.detection2d.ncc_high_sigma", 0.15f);
		this.ncc_low_mean					= (float) 	Prefs.get("critpoint.detection2d.ncc_low_mean", 0.7f);
		this.ncc_low_sigma					= (float) 	Prefs.get("critpoint.detection2d.ncc_low_sigma", 0.15f);

		this.likelihood_high_mean           = (float)   Prefs.get("critpoint.detection2d.likelihood_high_mean", 0.9f);
		this.likelihood_high_sigma 			= (float)   Prefs.get("critpoint.detection2d.likelihood_high_sigma", 0.2f);
		this.likelihood_low_mean			= (float) 	Prefs.get("critpoint.detection2d.likelihood_low_mean", 0.0f);
		this.likelihood_low_sigma			= (float) 	Prefs.get("critpoint.detection2d.likelihood_low_sigma", 0.2f);

		this.output_sigma					= (float) 	Prefs.get("critpoint.detection2d.output_sigma", 0.4f);

		GenericDialog gd = new GenericDialog("DETECTOR2D");

		gd.addNumericField("type: ", 				this.critpoint_type,	0,	2,	"(0-end, 1-junc)");
		gd.addNumericField("s: ", 					this.s,					1,	10,	"(scale)");
        gd.addStringField("Dlist", Dlist);
        gd.addNumericField("M ", 	                M, 			        	0,  10, "");
		gd.addNumericField("L ",                    L, 		            	0,  10, "");

        gd.addNumericField("NCC HIGH", 	        ncc_high_mean, 			1,  5, "MEAN");
        gd.addNumericField("NCC HIGH",        	ncc_high_sigma, 	    1,  5, "SIGMA");
		gd.addNumericField("NCC LOW",           ncc_low_mean, 		    1,  10, "MEAN");
		gd.addNumericField("NCC LOW",           ncc_low_sigma, 		    1,  10, "SIGMA");

		gd.addNumericField("LIKELIHOOD HIGH", 	likelihood_high_mean, 			1,  10, "MEAN");
		gd.addNumericField("LIKELIHOOD HIGH",   likelihood_high_sigma, 	    	1,  10, "SIGMA");
		gd.addNumericField("LIKELIHOOD LOW",    likelihood_low_mean, 		    1,  10, "MEAN");
		gd.addNumericField("LIKELIHOOD LOW",    likelihood_low_sigma, 		    1,  10, "SIGMA");

		gd.addNumericField("OUT ",    output_sigma, 		    1,  10, "SIGMA");

		gd.showDialog();
        if (gd.wasCanceled()) return DONE;

		critpoint_type = (int) gd.getNextNumber(); 			Prefs.set("critpoint.detection2d.critpoint_type", critpoint_type);
		this.s = (float) gd.getNextNumber(); 				Prefs.set("critpoint.detection2d.s", this.s);
        Dlist       	= gd.getNextString(); 				Prefs.set("critpoint.detection2d.d", Dlist);
        M       	    = (int) gd.getNextNumber(); 		Prefs.set("critpoint.detection2d.m", M);
		L		    	= (int) gd.getNextNumber();			Prefs.set("critpoint.detection2d.l", L);
		ncc_high_mean   = (float) gd.getNextNumber();   	Prefs.set("critpoint.detection2d.ncc_high_mean", ncc_high_mean);
		ncc_high_sigma  = (float) gd.getNextNumber();  		Prefs.set("critpoint.detection2d.ncc_high_sigma", ncc_high_sigma);
        ncc_low_mean    = (float) gd.getNextNumber();  		Prefs.set("critpoint.detection2d.ncc_low_mean", ncc_low_mean);
		ncc_low_sigma	= (float) gd.getNextNumber(); 		Prefs.set("critpoint.detection2d.ncc_low_sigma", ncc_low_sigma);
        likelihood_high_mean  	= (float) gd.getNextNumber(); Prefs.set("critpoint.detection2d.likelihood_high_mean", 	likelihood_high_mean);
        likelihood_high_sigma 	= (float) gd.getNextNumber(); Prefs.set("critpoint.detection2d.likelihood_high_sigma", 	likelihood_high_sigma);
        likelihood_low_mean 	= (float) gd.getNextNumber(); Prefs.set("critpoint.detection2d.likelihood_low_mean", 	likelihood_low_mean);
        likelihood_low_sigma 	= (float) gd.getNextNumber(); Prefs.set("critpoint.detection2d.likelihood_low_sigma", 	likelihood_low_sigma);

		output_sigma 	= (float) gd.getNextNumber(); Prefs.set("critpoint.detection2d.output_sigma", 	output_sigma);

		String[] dd = Dlist.split(",");
		D = new float[dd.length];
		for (int i=0; i<dd.length; i++) D[i] = Float.valueOf(dd[i]);

		if (critpoint_type==0) image_gndtth = image_name + ".end";
		else if (critpoint_type==1) image_gndtth = image_name + ".bif";
		else return DONE;

		output_membership_th = (float) Math.exp(-(0.5f*0.5f)/(2*output_sigma*output_sigma)); // 0.5 is middle of the output range
		output_membership_th = (float) (1 - Math.pow(output_membership_th,k) * (1-output_membership_th)); // h1 = 1 - h^k * (1-h)

		output_dir_name = image_dir+String.format("det_%d_%1.1f_%s_%d_%d__%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f",
														 critpoint_type, this.s, Dlist, M, L,
														 ncc_high_mean, ncc_high_sigma, ncc_low_mean, ncc_low_sigma,
														 likelihood_high_mean, likelihood_high_sigma, likelihood_low_mean, likelihood_low_sigma,
														 output_sigma);
		output_log_name = output_dir_name+File.separator+"detector2d.log";

		File f = new File(output_dir_name);
		if (!f.exists()) {

			f.mkdirs();

			try {
				logWriter = new PrintWriter(output_log_name);
				logWriter.print("");
				logWriter.close();
			} catch (FileNotFoundException ex) {}


		}
		// if it exists, just prepare to append on the existing file
		try {
			logWriter = new PrintWriter(new BufferedWriter(new FileWriter(output_log_name, true)));
		} catch (IOException e) {}

		CPU_NR = Runtime.getRuntime().availableProcessors() + 1;

		// output
		critpoint_det = new float[2][imagePlus.getWidth()][imagePlus.getHeight()]; // off and on values for selected critpoint (static - to be updated throughout the process)
        cnv = imagePlus.getCanvas();

        return DOES_8G+DOES_32+NO_CHANGES;

    }

	public void run(ImageProcessor imageProcessor){

		for (int didx = 0; didx<D.length; didx++) {

			IJ.log("scale " + D[didx]);
			long t1, t2;
        	t1 = System.currentTimeMillis();

			/********/
			Sphere2D sph2d = new Sphere2D(D[didx], s); //sph2d.showSampling().show(); sph2d.showWeights().show();
			/********/
        	IJ.log("extracting mask...");
			float nbhood_radius = 1.5f*sph2d.getOuterRadius();
        	Masker2D.loadTemplate(inimg_xy, sph2d.getOuterRadius(), nbhood_radius); //image, margin, check
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
			new ImagePlus("MASK,r="+IJ.d2s(nbhood_radius,2), Masker2D.getMask()).show();
			/********/
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
			/********/
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
			//PeakExtractor2D.getCircStat().show();
			/********/
			IJ.log("fitting model, extracting features, fuzzy detection...");
			Delineator2D.loadTemplate(Masker2D.i2xy, Masker2D.xy2i, PeakExtractor2D.peaks_i, PeakExtractor2D.peaks_w, PeakExtractor2D.peaks_lhood, inimg_xy,
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
//        int totalPeakAnalyzeComponents = Masker2D.i2xy.length; // number of locations
//        Delineator2D pa_jobs[] = new Delineator2D[CPU_NR];
//        for (int i = 0; i < pa_jobs.length; i++) {
//            pa_jobs[i] = new Delineator2D(i*totalPeakAnalyzeComponents/CPU_NR, (i+1)*totalPeakAnalyzeComponents/CPU_NR);
//            pa_jobs[i].start();
//        }
//        for (int i = 0; i < pa_jobs.length; i++) {
//            try {
//                pa_jobs[i].join();
//            } catch (InterruptedException e) {
//                e.printStackTrace();
//            }
//        }



			t2 = System.currentTimeMillis();
			IJ.log("done. " + ((t2 - t1) / 1000f) + "sec.");

		}

		logWriter.close();


	}

//    public void run(ImageProcessor imageProcessor) {
//        /*
//        IJ.log("simple detector...");
//        SimpleDetector2D.loadTemplate(inimg_xy.length, inimg_xy[0].length, Masker2D.i2xy, Masker2D.xy2i, Delineator2D.delin2, Delineator2D.lhood2);
//        int totalSimpleDetectComponents = Masker2D.i2xy.length;
//        SimpleDetector2D sd_jobs[] = new SimpleDetector2D[CPU_NR];
//        for (int i=0; i<sd_jobs.length; i++) {
//            sd_jobs[i] = new SimpleDetector2D(i*totalSimpleDetectComponents/CPU_NR, (i+1)*totalSimpleDetectComponents/CPU_NR);
//            sd_jobs[i].start();
//        }
//        for (int i=0; i<sd_jobs.length; i++) {
//            try {
//                sd_jobs[i].join();
//            } catch (InterruptedException e) {
//                e.printStackTrace();
//            }
//        }
//        */
//
////        ImageStack is_lhoods = Delineator2D.exportLikelihoods();
////        new ImagePlus("LHOODS2", is_lhoods).show();
//        // done extracting features and detecting
//        ImagePlus critpoint_lhood_img = Delineator2D.exportDetection();
//        critpoint_lhood_img.show();
//
//		//        Overlay overlay_with_detections = SimpleDetector2D.drawDetections();
////        ImagePlus final_det = cnv.getImage().duplicate();
////        final_det.show();
////        final_det.setOverlay(overlay_with_detections); // add the mto the original image
//
////        ImagePlus end_lhood = SimpleDetector2D.showEndLhood();
////        end_lhood.show();
////        ImagePlus jun_lhood = SimpleDetector2D.showJunctionLhood();
////        jun_lhood.show();
////        ImagePlus model_grad = SimpleDetector2D.showModelGradient();
////        model_grad.show();
//
//        // use Graphcut to get segment out stable regions // TODO put it in SimpleDetector2D
////        Graph_Cut gc = new Graph_Cut();
////        ImagePlus out = gc.processSingleChannelImage(jun_lhood, null, 500, 500, 500);
////        out.setTitle("GCsegmentation");
////        out.show();
//
//        /*
//        if (false) {
//
//            IJ.log("export features ");
//            t1 = System.currentTimeMillis();
//
//            String export_path = image_dir+image_name+".feat";
//            Delineator2D.exportFeats(export_path);
//            IJ.log("done exporting to: \t" + export_path + "\n");
//
//            export_path = image_dir+image_name+".desc";
//            Delineator2D.exportDescripts(export_path);
//            IJ.log("done exporting to: \t" + export_path + "\n");
//
//            export_path = image_dir + image_name + ".frame";
//            Delineator2D.exportFrames(export_path);
//            IJ.log("done exporting to: \t" + export_path + "\n");
//
//            export_path = image_dir + image_name + ".i2xy";
//            Masker2D.exportI2xyCsv(export_path);
//            IJ.log("done exporting to: \t" + export_path + "\n");
//
//            t2 = System.currentTimeMillis();
//            IJ.log("done. " + ((t2 - t1) / 1000f) + "sec.");
//
//        }
//        */
//
//		//************************************************
//
//        /*
//			mouse listeners after processing
//		 */
//        cnv.addMouseListener(this);
//        cnv.addMouseMotionListener(this);
//
//		IJ.selectWindow(cnv.getImage().getWindow().getTitle());
//		IJ.setTool("hand");
//        IJ.selectWindow(cnv.getImage().getTitle());
//
//    }

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