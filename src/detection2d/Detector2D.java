package detection2d;

import aux.ReadSWC;
import conn.Find_Connected_Regions;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.*;
import ij.io.FileSaver;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.*;
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
	static boolean show_bifpoints;
	static boolean show_endpoints;
	// model parameters
    float       s;
    String      Dlist;
    float       minCos = -.5f;   // hardcoded
    int         M;
	int         L;

//    float       masker_radius;
//    float       masker_percentile;

	float       sampling_crosswise = .3f; // hardcoded
	//
	float       ncc_high_mean;
    float       ncc_low_mean;
	float 		likelihood_high_mean;
	float 		likelihood_low_mean;
	float		output_sigma;   // Fuzzy2D that will be used assumes out ranges from 0 to 1

	//
	static float[]  	D;
	String 		image_gndtth_endpoints;		// ground truth swc file with critical points to evaluate, same name as input image
	String 		image_gndtth_bifurcations;	// ground truth swc file with critical points to evaluate, same name as input image
	float		k = 0.5f;                   // hardcoded, no need to tweak it each time, defines the stricktness of the threshold for out membership function
                                            // k higher => means more strictness
	float 		output_membership_th;       // based on k and output_sigma
	int			min_size_bif;                 // smallest connected region to be valid critical point
	int			max_size_bif;                 // smallest connected region to be valid critical point
	int			min_size_end;                 // smallest connected region to be valid critical point
	int			max_size_end;                 // smallest connected region to be valid critical point
	String		output_dir_name; 		    // parameter coded output folder name
	String 		output_log_name;
	PrintWriter logWriter = null;

    static  float   eval_scale  = 1.5f;
	String          eval_string = "";

    static float       min_D = Float.POSITIVE_INFINITY;
    static float       max_D = Float.NEGATIVE_INFINITY;

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
        this.show_bifpoints =               Prefs.get("critpoint.detection2d.show_bifpoints", true);
        this.show_endpoints =               Prefs.get("critpoint.detection2d.show_endpoints", true);

        this.s								= (float)	Prefs.get("critpoint.detection2d.s", 1.2f);
		this.Dlist 					        =    		Prefs.get("critpoint.detection2d.d", "4");
		this.M 					        	= (int)     Prefs.get("critpoint.detection2d.m", 1);
		this.L						        = (int) 	Prefs.get("critpoint.detection2d.l", 8);

//        this.masker_radius                  = (float)   Prefs.get("critpoint.detection2d.masker_radius", 4);
//        this.masker_percentile              = (float)   Prefs.get("critpoint.detection2d.masker_percentile", 50);

		this.ncc_high_mean                  = (float)   Prefs.get("critpoint.detection2d.ncc_high_mean", 1f);
		this.ncc_low_mean					= (float) 	Prefs.get("critpoint.detection2d.ncc_low_mean", 0.7f);
		this.likelihood_high_mean           = (float)   Prefs.get("critpoint.detection2d.likelihood_high_mean", 0.9f);
		this.likelihood_low_mean			= (float) 	Prefs.get("critpoint.detection2d.likelihood_low_mean", 0.0f);
		this.output_sigma					= (float) 	Prefs.get("critpoint.detection2d.output_sigma", 0.4f);

		GenericDialog gd = new GenericDialog("DETECTOR2D");

        gd.addCheckbox("BIFURCATIONS ", show_bifpoints);
        gd.addCheckbox("ENDPOINTS    ", show_endpoints);

        gd.addMessage("-- MODEL --");
        gd.addNumericField("s:", 					this.s,					1,	10,	"(scale)");
        gd.addStringField("Dlist", Dlist);
        gd.addNumericField("M", 	                M, 			        	0,  10, "");
		gd.addNumericField("L",                    L, 		            	0,  10, "");

//        gd.addMessage("-- MASKER --");
//        gd.addNumericField("radius:     ",          this.masker_radius,     1,  10, "pix");
//        gd.addNumericField("percentile: ",          this.masker_percentile, 1,  10, "[0-100]");

        gd.addMessage("-- FUZZY LOGIC SYSTEM --");
        gd.addNumericField("NCC_HIGH", 	        ncc_high_mean, 			2,  10, "");
		gd.addNumericField("NCC_LOW",           ncc_low_mean, 		    2,  10, "");
		gd.addNumericField("LIKELIHOOD_HIGH", 	likelihood_high_mean, 	2,  10, "");
		gd.addNumericField("LIKELIHOOD_LOW",    likelihood_low_mean, 	2,  10, "");
		gd.addNumericField("OUT",    			output_sigma, 		    2,  10, "SIGMA");

		gd.showDialog();
        if (gd.wasCanceled()) return DONE;

        this.show_bifpoints = gd.getNextBoolean();      	    Prefs.set("critpoint.detection2d.show_bifpoints", this.show_bifpoints);
        this.show_endpoints = gd.getNextBoolean();      	    Prefs.set("critpoint.detection2d.show_endpoints", this.show_endpoints);

		this.s = (float) gd.getNextNumber(); 				    Prefs.set("critpoint.detection2d.s", this.s);
        Dlist       	= gd.getNextString(); 				    Prefs.set("critpoint.detection2d.d", Dlist);
        M       	    = (int) gd.getNextNumber(); 		    Prefs.set("critpoint.detection2d.m", M);
		L		    	= (int) gd.getNextNumber();			    Prefs.set("critpoint.detection2d.l", L);

//        this.masker_radius = (float) gd.getNextNumber();        Prefs.set("critpoint.detection2d.masker_radius", this.masker_radius);
//        this.masker_percentile = (float) gd.getNextNumber();    Prefs.set("critpoint.detection2d.masker_percentile", this.masker_percentile);

		ncc_high_mean   = (float) gd.getNextNumber();   	    Prefs.set("critpoint.detection2d.ncc_high_mean", ncc_high_mean);
        ncc_low_mean    = (float) gd.getNextNumber();  		    Prefs.set("critpoint.detection2d.ncc_low_mean", ncc_low_mean);
        likelihood_high_mean  	= (float) gd.getNextNumber();   Prefs.set("critpoint.detection2d.likelihood_high_mean", 	likelihood_high_mean);
        likelihood_low_mean 	= (float) gd.getNextNumber();   Prefs.set("critpoint.detection2d.likelihood_low_mean", 	likelihood_low_mean);
		output_sigma 	= (float) gd.getNextNumber();           Prefs.set("critpoint.detection2d.output_sigma", 	output_sigma);

		String[] dd = Dlist.split(",");
		D = new float[dd.length];

		for (int i=0; i<dd.length; i++) {
            D[i] = Float.valueOf(dd[i]);
            if (D[i]<min_D) min_D = D[i];
            if (D[i]>max_D) max_D = D[i];
        }
        IJ.log(""+min_D+" <-> "+max_D);

		min_size_bif = Math.round(min_D);
		max_size_bif = (int) Math.round(min_D*min_D);
		//IJ.log("BIF "+min_size_bif+" <-> "+max_size_bif);


		min_size_end = Math.round(min_D);//(int) Math.round(0.5 * min_D);
        max_size_end = (int) Math.round(min_D * min_D);
		//IJ.log("END "+min_size_end+" <-> "+max_size_end);

		image_gndtth_endpoints = image_name + ".end";
		image_gndtth_bifurcations = image_name + ".bif";

		output_membership_th = (float) Math.exp(-(0.5f*0.5f)/(2*output_sigma*output_sigma)); // 0.5 is middle of the output range
		output_membership_th = (float) (1 - Math.pow(output_membership_th,k) * (1-output_membership_th)); // h1 = 1 - h^k * (1-h)

        //IJ.log(""+output_membership_th+ "   ");

		output_dir_name = image_dir+String.format(
														 "det(s,Dlist,M,L,nncH,nccL,lhoodH,lhoodL,outSig)_%.1f_%s_%d_%d_%.2f_%.2f_%.2f_%.2f_%.2f", //_%1.2f_%1.2f_%1.2f_%1.2f
														 this.s,
														 Dlist,
														 M, L,
														 ncc_high_mean,
														 ncc_low_mean,
														 likelihood_high_mean,
														 likelihood_low_mean,
														 output_sigma);
		output_log_name = output_dir_name+File.separator+"det.csv";

		File f = new File(output_dir_name);
		if (!f.exists()) {

			f.mkdirs();

			try {
				logWriter = new PrintWriter(output_log_name);
				logWriter.print("name,\tTP_BIF,\tFP_BIF,\tFN_BIF,\tP_BIF,\tR_BIF,\tTP_END,\tFP_END,\tFN_END,\tP_END,\tR_END\n");
				logWriter.close();
			} catch (FileNotFoundException ex) {}

		}
		// if it exists already in the folder, just prepare to append on the existing file
		try {
			logWriter = new PrintWriter(new BufferedWriter(new FileWriter(output_log_name, true)));
		} catch (IOException e) {}

		CPU_NR = Runtime.getRuntime().availableProcessors() + 1;

		// output
		critpoint_det = new float[3][imagePlus.getWidth()][imagePlus.getHeight()]; // off and on values for selected critpoint (static - to be updated throughout the process)
        cnv = imagePlus.getCanvas();

        return DOES_8G+DOES_32+NO_CHANGES;

    }

	public void run(ImageProcessor imageProcessor){

		for (int didx = 0; didx<D.length; didx++) {

			long t1, t2;
        	t1 = System.currentTimeMillis();

			/********/
			Sphere2D sph2d = new Sphere2D(D[didx], s); //sph2d.showSampling().show(); sph2d.showWeights().show();
			/********/
//			float nbhood_radius = 2f * D[didx];//1.0f*sph2d.getOuterRadius();
//			int   nbhood_margin = (int) Math.ceil(nbhood_radius);//2*(int) Math.ceil(nbhood_radius);
            float new_masker_radius = sph2d.getOuterRadius();   // important that it is outer radius of the sphere
            float new_masker_percentile = 50;                   // used to have these two as argument but not necessary
			Masker2D.loadTemplate(inimg_xy, (int)Math.ceil(new_masker_radius), new_masker_radius, new_masker_percentile); //image, margin, check, percentile
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
//			new ImagePlus("MASK,r="+IJ.d2s(new_masker_radius,2)+",per="+IJ.d2s(new_masker_percentile,0), Masker2D.getMask()).show();
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
			IJ.log("fitting model... extracting features... fuzzy detection...");
			Delineator2D.loadTemplate(
											 	Masker2D.i2xy,
											 	Masker2D.xy2i,
											 	PeakExtractor2D.peaks_i,
											 	PeakExtractor2D.peaks_w,
											 	PeakExtractor2D.peaks_lhood,
											 	inimg_xy,
												M,
												minCos,
												D[didx],
												ncc_high_mean,
												ncc_low_mean,
											 	likelihood_high_mean,
											 	likelihood_low_mean,
					                         	output_sigma,
											 	L,
												sampling_crosswise
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
			/********/
			IJ.log("append scores...");
			// take critpoint on/off scores from last class and append them to the critpoint_det[][][] - scale invariance
			for (int ll=0; ll<Delineator2D.critpoint_score2.length; ll++) {
				if (Delineator2D.critpoint_score2[ll]!=null) {

					float curr_endpoint 		= Delineator2D.critpoint_score2[ll][0];
					float curr_nonepoint 	    = Float.NEGATIVE_INFINITY;//Delineator2D.critpoint_score2[ll][1];
					float curr_bifpoint 		= Delineator2D.critpoint_score2[ll][2];
					float curr_max 		        = Math.max(curr_endpoint, Math.max(curr_nonepoint, curr_bifpoint));

					int corresponding_x = Masker2D.i2xy[ll][0];
					int corresponding_y = Masker2D.i2xy[ll][1];

					float stored_endpoint 	= critpoint_det[0][corresponding_x][corresponding_y];
					float stored_nonepoint 	= Float.NEGATIVE_INFINITY;//critpoint_det[1][corresponding_x][corresponding_y];
					float stored_bifpoint 	= critpoint_det[2][corresponding_x][corresponding_y];
					float stored_max	= Math.max(stored_endpoint, Math.max(stored_nonepoint, stored_bifpoint));

					if (curr_max>stored_max) { // it happened that one scale had better score - take it
						critpoint_det[0][corresponding_x][corresponding_y] = curr_endpoint;
						critpoint_det[1][corresponding_x][corresponding_y] = curr_nonepoint;
						critpoint_det[2][corresponding_x][corresponding_y] = curr_bifpoint;
					}
				}

			}

			t2 = System.currentTimeMillis();
			IJ.log("done. " + ((t2 - t1) / 1000f) + "sec.");

		}

        /* extract endpoints */
		eval_string = ""; // initialize before appending in following calls
		Overlay ov_bifurcations = exportBifurcationDetection();    //sequence bif-end is important for the logger
		Overlay ov_endpoints = exportEndpointDetection();
		logWriter.println(eval_string);
		IJ.log(""+eval_string);

		logWriter.close();

		ImagePlus final_det = cnv.getImage().duplicate();
		final_det.setTitle("det");

        Overlay ov = new Overlay();
        if (show_endpoints) for (int i_ov=0; i_ov<ov_endpoints.size(); i_ov++) ov.add(ov_endpoints.get(i_ov));
        if (show_bifpoints) for (int i_ov=0; i_ov<ov_bifurcations.size(); i_ov++) ov.add(ov_bifurcations.get(i_ov));

		final_det.setOverlay(ov);
//        final_det.show();

		// save tif with detection overlay
		FileSaver fs = new FileSaver(final_det);
		String save_path = output_dir_name + File.separator + image_name+".tif";
		fs.saveAsTiff(save_path);

        //IJ.log(save_path+ " saved.");

		// enable - disable inspection tools
		cnv.addMouseListener(this);
        cnv.addMouseMotionListener(this);
		IJ.setTool("hand");

	}

    public Overlay exportBifurcationDetection()
    {

        /* create yellow Overlay with BIFs */
        // create detection image
        int w = critpoint_det[2].length;
        int h = critpoint_det[2][0].length;

        byte[] t = new byte[w*h];
        for (int x=0; x<w; x++) {
            for (int y=0; y<h; y++) {
                float curr_on = critpoint_det[2][x][y];
                if (curr_on>=output_membership_th) {
                    t[y*w+x] = (byte) 255;
                }
            }
        }

        // take detections (binary image), find connected regions, and extract out the overlay with the detections
        ByteProcessor bp = new ByteProcessor(w, h, t);
        Find_Connected_Regions conn_reg = new Find_Connected_Regions(new ImagePlus("", bp), true);
        conn_reg.run("");
        //FileSaver fs = new FileSaver(conn_reg.showLabels());
        //fs.saveAsTiff(image_dir+"conn_bifs.tif");
        //IJ.log("connected components: "+image_dir+"conn_bifs.tif");

        Overlay ov = formPointOverlay(conn_reg.getConnectedRegions(), min_size_bif, max_size_bif, 1, 0, 0); // add red intensity based on average score

        /* before returning Overlay, check if there is a ground truth to compare the results */

        // compare dt with ground truth if it exists and print the results in a log file
        ArrayList<float[]> dt = formDetectionCentroids(conn_reg.getConnectedRegions(), min_size_bif);// xy stored

        if (new File(image_dir+image_gndtth_bifurcations).exists()) {

//            IJ.log("evaluating BIFURCATION detection using " + image_dir+image_gndtth_bifurcations);
            int tp = 0, fp = 0, fn = 0;

            ReadSWC reader = new ReadSWC(image_dir+image_gndtth_bifurcations);

            boolean[] annots = new boolean[reader.nodes.size()];

            for (int a=0; a<dt.size(); a++) {

                boolean found = false;

                for (int b=0; b<reader.nodes.size(); b++) {
                    float ax = dt.get(a)[0];
                    float ay = dt.get(a)[1];
                    float bx = reader.nodes.get(b)[reader.XCOORD];
                    float by = reader.nodes.get(b)[reader.YCOORD];
                    float dist = (float) Math.sqrt(Math.pow(bx-ax,2)+Math.pow(by-ay,2));
                    if (dist<=eval_scale*D[D.length-1]) {// compare the distance
                        found = true;
                        tp++;
                        annots[b] = true;
                    }
                }

                if (!found) {
                    fp++;
                }

            }

            for (int a=0; a<annots.length; a++) {
                if (!annots[a]) {
                    fn++;
                }
            }

            eval_string += "\""+image_name+"\",\t" + tp + ",\t" + fp + ",\t"+fn+ ",\t" + IJ.d2s(tp/(float)(tp+fp), 2) + ",\t" + IJ.d2s(tp/(float)(tp+fn),2) +",\t";

        }
        else {

            eval_string += "\""+image_name+"\",\t" + Float.NaN + ",\t" + Float.NaN + ",\t"+Float.NaN+ ",\t" + Float.NaN + ",\t" + Float.NaN+",\t";

        }

        return ov;

    }

	public Overlay exportEndpointDetection()
	{

        /* create yellow Overlay with ENDs */
		// create detection image
		int w = critpoint_det[0].length;
		int h = critpoint_det[0][0].length;

		byte[] t = new byte[w*h];
		for (int x=0; x<w; x++) {
			for (int y=0; y<h; y++) {
				float curr_on = critpoint_det[0][x][y];
				if (curr_on>=output_membership_th) {
					t[y*w+x] = (byte) 255;
				}
			}
		}

		// take detections (binary image), find connected regions, and extract out the overlay with the detections
		ByteProcessor bp = new ByteProcessor(w, h, t);
		Find_Connected_Regions conn_reg = new Find_Connected_Regions(new ImagePlus("END,th="+IJ.d2s(output_membership_th,2), bp), true);
		conn_reg.run("");
		//conn_reg.showLabels().show();
		//FileSaver fs = new FileSaver(conn_reg.showLabels());
		//fs.saveAsTiff(image_dir+"conn_ends.tif");

		Overlay ov = formPointOverlay(conn_reg.getConnectedRegions(), min_size_end, max_size_end, 1, 1, 0); // add yellow intensity based on average score

        /* before returning Overlay, check if there is a ground truth to compare the results */

        // compare dt with ground truth if it exists and print the results in a log file
        ArrayList<float[]> dt = formDetectionCentroids(conn_reg.getConnectedRegions(), min_size_end);// xy stored

		if (new File(image_dir+image_gndtth_endpoints).exists()) {

//            IJ.log("evaluating ENDPOINT detection using " + image_dir+image_gndtth_endpoints);
			int tp = 0, fp = 0, fn = 0;

			ReadSWC reader = new ReadSWC(image_dir+image_gndtth_endpoints);

			boolean[] annots = new boolean[reader.nodes.size()];

			for (int a=0; a<dt.size(); a++) {

				boolean found = false;

				for (int b=0; b<reader.nodes.size(); b++) {
					float ax = dt.get(a)[0];
					float ay = dt.get(a)[1];
					float bx = reader.nodes.get(b)[reader.XCOORD];
					float by = reader.nodes.get(b)[reader.YCOORD];
					float dist = (float) Math.sqrt(Math.pow(bx-ax,2)+Math.pow(by-ay,2));
					if (dist<=eval_scale*D[D.length-1]) {// compare the distance
						found = true;
						tp++;
						annots[b] = true;
					}
				}

				if (!found) {
					 fp++;
				}

			}

			for (int a=0; a<annots.length; a++) {
				if (!annots[a]) {
					 fn++;
				}
			}

			eval_string += tp + ",\t" + fp + ",\t"+fn+ ",\t" + IJ.d2s(tp/(float)(tp+fp),2) + ",\t" + IJ.d2s(tp/(float)(tp+fn),2);

		}
		else {
            eval_string += Float.NaN + ",\t" + Float.NaN + ",\t"+ Float.NaN + ",\t" + Float.NaN + ",\t" + Float.NaN;
		}

		return ov;

	}

	private static ArrayList<float[]> formDetectionCentroids(ArrayList<ArrayList<int[]>> regs, int minSize)
	{

		ArrayList<float[]> lst = new ArrayList<float[]>();

		for (int i=0; i<regs.size(); i++) {
			if (regs.get(i).size()>minSize) {

				float Cx=0, Cy=0;//, R= (float) Math.sqrt((float)regs.get(i).size()/Math.PI);
//				R = (R<1)? 1 : R ;

				for (int aa=0; aa<regs.get(i).size(); aa++) {
					Cx += regs.get(i).get(aa)[1];
					Cy += regs.get(i).get(aa)[0];
				}

				Cx /= regs.get(i).size();
				Cy /= regs.get(i).size();

				lst.add(new float[]{Cx, Cy});

			}
		}

		return lst;

	}

	private static Overlay formPointOverlay(ArrayList<ArrayList<int[]>> regs, int minSize, int maxSize, float r_col, float g_col, float b_col)
	{

		Overlay detections = new Overlay();

		for (int i=0; i<regs.size(); i++) {
			if (regs.get(i).size()>=minSize && regs.get(i).size()<maxSize) {

				float Cx=0, Cy=0, R = (float) Math.sqrt((float)regs.get(i).size()/Math.PI);
//				R = (R<1)? 1 : 1.5f*R ; // because just 1*R is too small for visualizaiton

				for (int aa=0; aa<regs.get(i).size(); aa++) {
					Cx += regs.get(i).get(aa)[1];
					Cy += regs.get(i).get(aa)[0];
				}

				Cx /= regs.get(i).size();
				Cy /= regs.get(i).size();

				// maybe nice to visualize wrt diameters set and not the size of the connected region

                R = 0.5f * max_D;//D[D.length-1]; // D is static
				OvalRoi ovroi = new OvalRoi(Cx-R+.5, Cy-R+.5, 2*R, 2*R);
				ovroi.setStrokeWidth(3);
				ovroi.setFillColor(new Color(r_col, g_col, b_col, 1.0f));
//				if (critpoint_type==0)
//					ovroi.setFillColor(Color.YELLOW);
//				else if (critpoint_type==1)
//					ovroi.setFillColor(Color.RED);

//				int firstY = regs.get(i).get(0)[0];
//				int firstX = regs.get(i).get(0)[1];
//                if (score2[firstY*W+firstX]==(byte)1) ovroi.setStrokeColor(Color.YELLOW);
//                else ovroi.setStrokeColor(Color.RED);
				detections.add(ovroi);

			}
		}

		return detections;

	}


//    public void run(ImageProcessor imageProcessor) {
//        /*
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
//			//Overlay overlay_with_detections = SimpleDetector2D.drawDetections();
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
//		IJ.selectWindow(cnv.getImage().getWindow().getTitle());
//        IJ.selectWindow(cnv.getImage().getTitle());
//    }

    public void mouseClicked(MouseEvent e)
    {

//        mouseMoved(e);

        int clickX = cnv.offScreenX(e.getX());
        int clickY = cnv.offScreenY(e.getY());

		/*
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

        pfl_is = Delineator2D.getFuzzyAgg(clickX, clickY);
        pfl_im.setStack(pfl_is);
        pfl_im.show();


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