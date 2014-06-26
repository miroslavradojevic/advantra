package detection2d;

import aux.ClusterDirections;
import aux.ReadSWC;
import aux.Stat;
import conn.Find_Connected_Regions;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Line;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.io.FileSaver;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;

import java.awt.*;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by miroslav on 16-5-14.
 */
public class Detector2D {

	String 		image_dir;
	String		image_name;
	float[][] 	inimg_xy;               	// store input image as an array
	String 		image_gndtth_endpoints;		// ground truth swc file with critical points to evaluate, same name as input image
	String 		image_gndtth_bifurcations;
	String 		image_gndtth_crosssections;

	// parameters
	float       	s;
	float 			sigma_ratio; 			// sigma = sigma_ratio * D
	static float[]  D;
	float       	ncc_high;
	float       	ncc_low;

    float 			likelihood_high;
	float 			likelihood_low;

    float           smoothness_high;            // smoothness_high is actually lower than smoothness_low
    float           smoothness_low;             // they are automatically calculated fro mthe distribution of smoothness values

    float			output_sigma;

	float 			output_membership_th;       // based on k and output_sigma



	int         	M = 1;
	float       	minCos = -.5f;
	float			k = 0.8f;                   // hardcoded, defines the stricktness of the threshold for out membership function, k higher => means more strictness
	String          eval_string = "";
	public boolean 	save_midresults = false;
	String midresults_dir = "";
    public boolean  auto_smoothness = false;

	String		output_dir_name; 		    	// parameter coded output folder name
	String 		output_log_name;
	PrintWriter logWriter = null;

    Fuzzy2D sample_fls = null; // will be used for user interaction (to simulate detection for features extracted at single locations)

	int         CPU_NR;

	// output
	static float[][] critpoint_det; 				// store 2d grid with detection scores (size of the image itself)
	ArrayList<CritpointRegion> detected_regions; 	// list with the detections
	static int[][]   region_map;					// pixel at (x,y) assigned to one region from the list by storing its index

    public Detector2D(
							 ImagePlus 	ip_load,
							 float 		_s,
							 float[] 	_D,
							 float 		_sigma_ratio,
							 float 		_ncc_high,
							 float 		_ncc_low,
							 float 		_likelihood_high,
							 float 		_likelihood_low,
                             float      _smoothness_high,
                             float      _smoothness_low,
							 float 		_output_sigma
    )
	{

		image_dir = ip_load.getOriginalFileInfo().directory; //  + File.separator  + image_name
		image_name = ip_load.getShortTitle();

		inimg_xy = new float[ip_load.getWidth()][ip_load.getHeight()]; // x~column, y~row
		if (ip_load.getType()== ImagePlus.GRAY8) {
			byte[] read = (byte[]) ip_load.getProcessor().getPixels();
			for (int idx=0; idx<read.length; idx++) {
				inimg_xy[idx%ip_load.getWidth()][idx/ip_load.getWidth()] = (float) (read[idx] & 0xff);
			}
		}
		else if (ip_load.getType()==ImagePlus.GRAY32) {
			float[] read = (float[]) ip_load.getProcessor().getPixels();
			for (int idx=0; idx<read.length; idx++) {
				inimg_xy[idx%ip_load.getWidth()][idx/ip_load.getWidth()] = read[idx];
			}
		}
		else {
			IJ.log("image type not recognized");
			return;
		}

		image_gndtth_endpoints 		= image_name + ".end";
		image_gndtth_bifurcations 	= image_name + ".bif";
		image_gndtth_crosssections 	= image_name + ".crs";

		// detection scale define
		s = _s;
		D = _D;
		sigma_ratio = _sigma_ratio;

		// detection parameters
		ncc_high = _ncc_high;
		ncc_low = _ncc_low;
		likelihood_high = _likelihood_high;
		likelihood_low = _likelihood_low;
        smoothness_high = _smoothness_high;
        smoothness_low = _smoothness_low;
		output_sigma = _output_sigma;

		String		Dlist = ""; // for the out directory
		for (int i=0; i<_D.length; i++) Dlist += IJ.d2s(_D[i], 1) + ((i == _D.length - 1) ? "" : ":");

		output_membership_th = (float) Math.exp(-(0.5f*0.5f)/(2*output_sigma*output_sigma)); // 0.5 is middle of the output range
		output_membership_th = (float) (1 - Math.pow(output_membership_th,k) * (1-output_membership_th)); // h1 = 1 - h^k * (1-h)

		output_dir_name = image_dir+String.format(
														 "det.Dlist.nncH.nccL.lhoodH.lhoodL_%s_%.2f_%.2f_%.2f_%.2f",
														 Dlist,
														 ncc_high,
														 ncc_low,
														 likelihood_high,
														 likelihood_low
														 );
		output_log_name = output_dir_name+ File.separator+"det.csv";
		midresults_dir = image_dir+image_name + "_midresults" + File.separator;

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

		File f1 = new File(midresults_dir);
		if (!f1.exists()) {
			f1.mkdirs();
		}
        else { // clean it, reset with each exec
            File[] files = f1.listFiles();
            for (File file : files)
            {
                if (!file.delete())
                    System.out.println("Failed to delete " + file);
            }
        }

		CPU_NR = Runtime.getRuntime().availableProcessors() + 1;

		critpoint_det = new float[ip_load.getWidth()][ip_load.getHeight()];
		detected_regions = new ArrayList<CritpointRegion>();

        sample_fls = new Fuzzy2D( // default initialization with given smoothness (yet to calculate auto)
                ncc_high, // ncc
                ncc_low,
                likelihood_high, // lhood
                likelihood_low,
                smoothness_high, // smoothness default
                smoothness_low,
                output_sigma  // std output membership functions - defines separation
        );

    }

    public void run()
	{

		long t1, t2;
		t1 = System.currentTimeMillis();

		for (int didx = 0; didx<D.length; didx++) { // loop scales D[]


			System.out.print("\nD = " + D[didx] + " : ");
			Sphere2D sph2d = new Sphere2D(D[didx], s, sigma_ratio);
            if (save_midresults) {
                IJ.saveAs(sph2d.showSampling(), "Tiff", midresults_dir+"sampling_"+D[didx]+".tif");
                IJ.saveAs(sph2d.showWeights(),  "Tiff", midresults_dir+"weights_"+D[didx]+".tif");
//                sph2d.showWeights().show();
            }
			/********************************************************************/
			System.out.print("Masker2D...");
			float new_masker_radius = 1.5f*sph2d.getOuterRadius();   	// important that it is outer radius of the sphere
			float new_masker_percentile = 30;                   		// used to have these two as argument but not necessary
			Masker2D.loadTemplate(
										 inimg_xy,
										 (int)Math.ceil(new_masker_radius),
										 new_masker_radius,
										 new_masker_percentile); //image, margin, check, percentile
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
			if (save_midresults) {
				ImagePlus mask = new ImagePlus("mask", Masker2D.getMask());
				IJ.saveAs(mask, "Tiff", midresults_dir+"mask_"+D[didx]+".tif");
			}
			/********************************************************************/
			System.out.print("Profiler2D...");
			Profiler2D.loadTemplate(
                    sph2d,
                    Masker2D.i2xy,
                    Masker2D.xy2i,
                    inimg_xy);
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
			/********************************************************************/
			System.out.print("ProfileSpread2D...");
			ProfileSpread2D.loadTemplate(
                    Masker2D.i2xy,
                    Masker2D.xy2i,
                    Profiler2D.prof2,
                    inimg_xy.length,
                    inimg_xy[0].length);
			int totalProfileLocations = Profiler2D.prof2.length;
			ProfileSpread2D pv_jobs[] = new ProfileSpread2D[CPU_NR];
			for (int i = 0; i < pv_jobs.length; i++) {
				pv_jobs[i] = new ProfileSpread2D(i*totalProfileLocations/CPU_NR,  (i+1)*totalProfileLocations/CPU_NR);
				pv_jobs[i].start();
			}
			for (int i = 0; i < pv_jobs.length; i++) {
				try {
					pv_jobs[i].join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			ProfileSpread2D.threshold();
			if (save_midresults) {
				ImagePlus testip = new ImagePlus("", ProfileSpread2D.getMask());
				IJ.saveAs(testip, "Tiff", midresults_dir+"mask_profile_"+D[didx]+".tif");
			}
			System.out.print(" " + IJ.d2s((ProfileSpread2D.getNrCritpointCandidates() * 100f) / (inimg_xy.length * inimg_xy[0].length), 2) + " % candidates...");
			/********************************************************************/
			System.out.print("PeakExtractor2D...");
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
			/********************************************************************/
			System.out.print("Delineator2D...");
			Delineator2D.loadTemplate(
											 Masker2D.i2xy,
											 Masker2D.xy2i,
											 ProfileSpread2D.profile_diverse,
											 PeakExtractor2D.peaks_i,
											 PeakExtractor2D.peaks_w,
											 inimg_xy,
											 sph2d,	//D[didx],
											 M,
											 minCos
			);
			int totalDelineationComponents = PeakExtractor2D.peaks_i.length;
			Delineator2D dl_jobs[] = new Delineator2D[CPU_NR];
			for (int i = 0; i < dl_jobs.length; i++) {
				dl_jobs[i] = new Delineator2D(i*totalDelineationComponents/CPU_NR, (i+1)*totalDelineationComponents/CPU_NR);
				dl_jobs[i].start();
			}
			for (int i = 0; i < dl_jobs.length; i++) {
				try {
					dl_jobs[i].join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			if (save_midresults) {  // TODO
				//new ImagePlus("", Delineator2D.getSmoothnessDistribution(64)).show();
			}

            if (auto_smoothness) {
                int percentile = 90;
                float sensitivity = 0.5f;
                smoothness_high = Delineator2D.getSmoothnessPercentile(percentile);
                smoothness_low = sensitivity * smoothness_high;
                System.out.print("new smoothness_high: " + smoothness_high);
                System.out.print("new smoothness_low:  " + smoothness_low);
                System.out.println("new fls...");
                sample_fls = new Fuzzy2D(// redefine flas used later for simulating, for visualizations (should be the same as the one in FuzzyDetector2D.run())
                        ncc_high, // ncc
                        ncc_low,
                        likelihood_high, // lhood
                        likelihood_low,
                        smoothness_high, // smoothness redefined
                        smoothness_low,
                        output_sigma  // std output membership functions - defines separation
                );
            }
            /********************************************************************/
			System.out.print("Ncc2D...");
			Ncc2D.loadTemplate(
									  Delineator2D.xy2,
									  Delineator2D.vxy2,
									  Delineator2D.L,
									  Delineator2D.dim,
									  Delineator2D.samplingStep,
									  inimg_xy,
									  sigma_ratio
			);
			int totalNccExtractionComponents = Delineator2D.xy2.length;
			Ncc2D ncc_jobs[] = new Ncc2D[CPU_NR];
			for (int i = 0; i < ncc_jobs.length; i++) {
				ncc_jobs[i] = new Ncc2D(i*totalNccExtractionComponents/CPU_NR, (i+1)*totalNccExtractionComponents/CPU_NR);
				ncc_jobs[i].start();
			}
			for (int i = 0; i < ncc_jobs.length; i++) {
				try {
					ncc_jobs[i].join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
            /********************************************************************/
            System.out.print("FuzzyDetector2D...");
            //float[] ab = Stat.get_min_max(FuzzyDetector2D.endpoint_score);
            //System.out.println("MN/MX="+ Arrays.toString(ab));

            FuzzyDetector2D.loadTemplate(
                    Ncc2D.scores,
                    PeakExtractor2D.peaks_lhood,
                    Delineator2D.smoothness,

                    ncc_high,        // user
                    ncc_low,         // user

                    likelihood_high,  // user
                    likelihood_low,   // user

                    smoothness_high,  // automatic or manual
                    smoothness_low,   // automatic or manual

                    output_sigma
            );
            int totalFuzzyDetectorComponents = PeakExtractor2D.peaks_lhood.length;
            FuzzyDetector2D fd_jobs[] = new FuzzyDetector2D[CPU_NR];
            for (int i = 0; i < fd_jobs.length; i++) {
                fd_jobs[i] = new FuzzyDetector2D(i*totalFuzzyDetectorComponents/CPU_NR, (i+1)*totalFuzzyDetectorComponents/CPU_NR);
                fd_jobs[i].start();
            }
            for (int i = 0; i < fd_jobs.length; i++) {
                try {
                    fd_jobs[i].join();
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }
            /********************************************************************/
			/*System.out.print("FeatureExtractor2D...");
			FeatureExtractor2D.loadTemplate(
												   Masker2D.i2xy,
												   Masker2D.xy2i,
												   PeakExtractor2D.peaks_i,
												   PeakExtractor2D.peaks_w,
												   PeakExtractor2D.peaks_lhood,
												   inimg_xy,
												   M,
												   minCos,
												   D[didx],
												   sigma_ratio,
												   "MEAN",
												   ncc_high,
												   ncc_low,
												   likelihood_high,
												   likelihood_low,
												   output_sigma
			);
			int totalPeakAnalyzeComponents = Masker2D.i2xy.length; // number of locations
			FeatureExtractor2D pa_jobs[] = new FeatureExtractor2D[CPU_NR];
			for (int i = 0; i < pa_jobs.length; i++) {
				pa_jobs[i] = new FeatureExtractor2D(i*totalPeakAnalyzeComponents/CPU_NR, (i+1)*totalPeakAnalyzeComponents/CPU_NR);
				pa_jobs[i].start();
			}
			for (int i = 0; i < pa_jobs.length; i++) {
				try {
					pa_jobs[i].join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			*/
			/********************************************************************/
			System.out.print("AppendRegions...");
			float ang_deg = 20f;  // ms parameter
            int ridiculously_low = 2;
            // with D we expect to capture ~(D/2)^2*pi pixel area
			float 		region_radius 	= D[didx];
            float       target_area = Math.round(Math.pow((region_radius/2),2) * Math.PI);
            int			min_region_size = (int) Math.ceil(0.1*target_area);  // (int) Math.round(0.25f*0.25f*Math.pow(D[didx],2));  // smallest connected region area to be valid critical point
            int			max_region_size = Math.round(5*target_area);    // (int) Math.round(2.0f*2.0f*Math.pow(D[didx],2));  // largest connected region to be valid critical point
            min_region_size = (min_region_size<ridiculously_low)? ridiculously_low : min_region_size;

            spatialScoreDistribution(critpoint_det, FuzzyDetector2D.endpoint_score, Masker2D.i2xy);
            appendCritpointRegions(
																					critpoint_det,
																					min_region_size,
																					max_region_size,
									   												region_radius,	          // the same radius will be assigned at each scale
																					CritpointRegion.RegionType.END,
																					ang_deg,
																					detected_regions
																					);
            System.out.println(" " + detected_regions.size() + " endpoints found ("+min_region_size+"-"+max_region_size+")");
//			FeatureExtractor2D.exportCritpointScores(critpoint_det, 2); // store critpoint_score2 list into 2d array critpoint_det
//			getCritpointRegions(
//																					critpoint_det,
//																					min_region_size,
//																					max_region_size,
//									   												region_radius,
//																					CritpointRegion.RegionType.BIF_CROSS,
//																					ang_deg,
//																					detected_regions
//																					);
			System.out.println(" " + didx + "/" + D.length );

        }

		t2 = System.currentTimeMillis();
		IJ.log("done. " + ((t2 - t1) / 1000f) + "sec.");

	}

	public void doEvaluation() // once the  list is completed useful to compare
	{

		// loop detected_regions to evaluate detection wrt .bif .end or .crs ground truth file at the same location
		// as the image evaluation is done per bifurcations, endpoints, cross-points and junctions (bifurcations and cross-points)
		// string format:
		// file_name
			// tp_BIF fp_BIF fn_BIF
				// tp_END fp_END fn_END
					// tp_CRS fp_CRS fn_CRS
						// tp_JUN fp_JUN fn_JUN

		eval_string = ""; // initialize before appending in following calls

		int tp, fp, fn;
		ReadSWC reader;
		boolean[] annots;  // necessary for fn calculation

		boolean bif_annot_exists = new File(image_dir+image_gndtth_bifurcations).exists();
		boolean end_annot_exists = new File(image_dir+image_gndtth_endpoints).exists();
		boolean crs_annot_exists = new File(image_dir+image_gndtth_crosssections).exists();

		if (bif_annot_exists && end_annot_exists && crs_annot_exists) {

			// BIF
			tp = 0; fp = 0; fn = 0;
			reader = new ReadSWC(image_dir+image_gndtth_bifurcations);
			annots = new boolean[reader.nodes.size()];

			// loop all detected BIF regions (all sorts of critical points are in the same list now)
			for (int a=0; a<detected_regions.size(); a++) {
				if (detected_regions.get(a).type== CritpointRegion.RegionType.BIF) {


					boolean found = false;

					float ax = detected_regions.get(a).centroid[0];
					float ay = detected_regions.get(a).centroid[1];
					float ar = detected_regions.get(a).radius;

					for (int b=0; b<reader.nodes.size(); b++) {

						float bx = reader.nodes.get(b)[reader.XCOORD];
						float by = reader.nodes.get(b)[reader.YCOORD];
						float br = reader.nodes.get(b)[reader.RADIUS];

						if (circlesOverlap(ax, ay, ar, bx, by, br)) {// two roi circles overlap
							found = true;
							tp++;
							annots[b] = true;
						}
					}

					if (!found) fp++;  // detected but was not in the list of annotated ones


				}
			}

			for (int a=0; a<annots.length; a++) if (!annots[a]) fn++;
			eval_string += "\""+image_name+"\",\t" + tp + ",\t" + fp + ",\t"+fn+ ",\t";


			// END
			tp = 0; fp = 0; fn = 0;
			reader = new ReadSWC(image_dir+image_gndtth_endpoints);
			annots = new boolean[reader.nodes.size()];

			for (int a=0; a<detected_regions.size(); a++) {
				if (detected_regions.get(a).type== CritpointRegion.RegionType.END) {

					boolean found = false;

					float ax = detected_regions.get(a).centroid[0];
					float ay = detected_regions.get(a).centroid[1];
					float ar = detected_regions.get(a).radius;

					for (int b=0; b<reader.nodes.size(); b++) {

						float bx = reader.nodes.get(b)[reader.XCOORD];
						float by = reader.nodes.get(b)[reader.YCOORD];
						float br = reader.nodes.get(b)[reader.RADIUS];

						if (circlesOverlap(ax, ay, ar, bx, by, br)) {// two roi circles overlap
							found = true;
							tp++;
							annots[b] = true;
						}
					}

					if (!found) fp++;  // detected but was not in the list of annotated ones

				}
			}

			for (int a=0; a<annots.length; a++) if (!annots[a]) fn++;
			eval_string += tp + ",\t" + fp + ",\t" + fn + ",\t";




			// CRS
			tp = 0; fp = 0; fn = 0;
			reader = new ReadSWC(image_dir+image_gndtth_crosssections);
			annots = new boolean[reader.nodes.size()];

			for (int a=0; a<detected_regions.size(); a++) {
				if (detected_regions.get(a).type== CritpointRegion.RegionType.CROSS) {

					boolean found = false;

					float ax = detected_regions.get(a).centroid[0];
					float ay = detected_regions.get(a).centroid[1];
					float ar = detected_regions.get(a).radius;

					for (int b=0; b<reader.nodes.size(); b++) {

						float bx = reader.nodes.get(b)[reader.XCOORD];
						float by = reader.nodes.get(b)[reader.YCOORD];
						float br = reader.nodes.get(b)[reader.RADIUS];

						if (circlesOverlap(ax, ay, ar, bx, by, br)) {// two roi circles overlap
							found = true;
							tp++;
							annots[b] = true;
						}
					}

					if (!found) fp++;  // detected but was not in the list of annotated ones

				}
			}

			for (int a=0; a<annots.length; a++) if (!annots[a]) fn++;
			eval_string += tp + ",\t" + fp + ",\t" + fn + ",\t";





			// JUN (BIF+CRS)
			tp = 0; fp = 0; fn = 0;
			ReadSWC readerBIF = new ReadSWC(image_dir+image_gndtth_bifurcations);
			ReadSWC readerCRS = new ReadSWC(image_dir+image_gndtth_crosssections);
			annots = new boolean[readerBIF.nodes.size() + readerCRS.nodes.size()];

			for (int a=0; a<detected_regions.size(); a++) {
				if (detected_regions.get(a).type== CritpointRegion.RegionType.BIF || detected_regions.get(a).type== CritpointRegion.RegionType.CROSS) {

					boolean found = false;

					float ax = detected_regions.get(a).centroid[0];
					float ay = detected_regions.get(a).centroid[1];
					float ar = detected_regions.get(a).radius;

					for (int b=0; b<readerBIF.nodes.size(); b++) {

						float bx = readerBIF.nodes.get(b)[reader.XCOORD];
						float by = readerBIF.nodes.get(b)[reader.YCOORD];
						float br = readerBIF.nodes.get(b)[reader.RADIUS];

						if (circlesOverlap(ax, ay, ar, bx, by, br)) {// two roi circles overlap
							found = true;
							tp++;
							annots[b] = true;
						}
					}

					for (int b = readerBIF.nodes.size(); b<readerBIF.nodes.size()+readerCRS.nodes.size(); b++) {

						float bx = readerCRS.nodes.get(b-readerBIF.nodes.size())[reader.XCOORD];
						float by = readerCRS.nodes.get(b-readerBIF.nodes.size())[reader.YCOORD];
						float br = readerCRS.nodes.get(b-readerBIF.nodes.size())[reader.RADIUS];

						if (circlesOverlap(ax, ay, ar, bx, by, br)) {// two roi circles overlap
							found = true;
							tp++;
							annots[b] = true;
						}

					}

					if (!found) fp++;  // detected but was not in the list of annotated ones

				}
			}

			for (int a=0; a<annots.length; a++) if (!annots[a]) fn++;
			eval_string += tp + ",\t" + fp + ",\t" + fn;



		}   // incomplete annotation
		else
			eval_string += "\"" + image_name + "\",\t" +
								   Float.NaN + ",\t" + Float.NaN + ",\t" + Float.NaN + ",\t" +    // bif
								   Float.NaN + ",\t" + Float.NaN + ",\t" + Float.NaN + ",\t" +    // end
								   Float.NaN + ",\t" + Float.NaN + ",\t" + Float.NaN + ",\t" +    // crs
								   Float.NaN + ",\t" + Float.NaN + ",\t" + Float.NaN;             // jun (bif+crs)

		logWriter.println(eval_string);
		IJ.log(""+eval_string);

		logWriter.close();

	}

	private boolean circlesOverlap(float x1, float y1, float r1, float x2, float y2, float r2)
	{
		return Math.pow(x1-x2,2)+Math.pow(y1-y2,2) <= Math.pow(r1+r2,2);
	}

	public Overlay getDetectionOverlay()
	{

		// create Overlay with CritpointRegions
		Overlay ov = new Overlay();

		for (int i=0; i<detected_regions.size(); i++) {

            if (detected_regions.get(i) != null) {

                float cx = detected_regions.get(i).centroid[0];
                float cy = detected_regions.get(i).centroid[1];
                float cr = detected_regions.get(i).radius;
                float sc = detected_regions.get(i).score;
                CritpointRegion.RegionType ctype = detected_regions.get(i).type;

                Color region_color = null;
                switch (ctype) {
                    case BIF:
                        region_color = new Color(1, 0, 0, sc);
                        break;

                    case END:
                        region_color = new Color(1, 1, 0, sc);
                        break;

                    case CROSS:
                        region_color = new Color(0, 1, 0, sc);
                        break;

                    default:
                        IJ.log("non valid critical point");
                        break;
                }

                //add region
                OvalRoi ovroi = new OvalRoi(cx-cr+.5, cy-cr+.5, 2*cr, 2*cr);
                ovroi.setStrokeWidth(1);
//			    ovroi.setFillColor(region_color);
                ov.add(ovroi);

                // add directions
                for (int j = 0; j < detected_regions.get(i).outward_directions.length; j++) {
                    float dx = detected_regions.get(i).outward_directions[j][0];
                    float dy = detected_regions.get(i).outward_directions[j][1];
                    Line l = new Line(cx+.5f, cy+.5f, cx+1.5*cr*dx+.5f, cy+1.5*cr*dy+.5f);
                    l.setStrokeWidth(3);
                    l.setStrokeColor(region_color);
                    l.setFillColor(region_color);
                    ov.add(l);
                }

            }

		}

		return ov;

	}

	// appends to the list of CritpointRegion, using connected components to optimize region segmentation
	private void appendCritpointRegions(
																  float[][] _critpoint_det,
																  int minSize,
																  int maxSize,
																  float det_radius,
																  CritpointRegion.RegionType _choose_type,
																  // mean-shift param (for directionality analysis)
																  float alfa_deg,
																  // destination to append
																  ArrayList<CritpointRegion> region_list
																  )
	{

		// create detection image
		int w = _critpoint_det.length;
		int h = _critpoint_det[0].length;

        // sto
		byte[] t = new byte[w*h];
		for (int x=0; x<w; x++) {
			for (int y=0; y<h; y++) {
				float curr_on = _critpoint_det[x][y];
                // threshold score on critpoint detection
                if (curr_on>=output_membership_th) t[y * w + x] = (byte) 255;
			}
		}

		// take detections (binary image), find connected regions
		ByteProcessor bp = new ByteProcessor(w, h, t);
        ImagePlus imp_thresholded_scores = new ImagePlus("", bp);
		Find_Connected_Regions conn_reg = new Find_Connected_Regions(imp_thresholded_scores, true);
		conn_reg.run("");
		//conn_reg.showLabels().show();

        if (save_midresults) {
			if (_choose_type == CritpointRegion.RegionType.END) {
                IJ.saveAs(imp_thresholded_scores,                                   "Tiff", midresults_dir+"ends_thr_"+det_radius+".tif");
                IJ.saveAs(new ImagePlus("", new FloatProcessor(_critpoint_det)),    "Tiff", midresults_dir+"ends_sco_"+det_radius+".tif");
            }
			if (_choose_type == CritpointRegion.RegionType.BIF_CROSS) {
                IJ.saveAs(imp_thresholded_scores,                                   "Tiff", midresults_dir+"jun_thr_"+det_radius+".tif");
                IJ.saveAs(new ImagePlus("", new FloatProcessor(_critpoint_det)),    "Tiff", midresults_dir+"jun_sco_"+det_radius+".tif");
            }
		}

		if (_choose_type== CritpointRegion.RegionType.END) {// DEBUG: just thresholding and connected components will give a lot of noisy regions
			// export binary image
//			ByteProcessor spots = new ByteProcessor(_critpoint_det.length, _critpoint_det[0].length);
//			System.out.println("Endpoints... what happened? total " + regs.size() + "regions" + minSize + " <> " +maxSize);
//			int cnt = 0;
//			for (int i=0; i<regs.size(); i++) {
//				if (true) { // regs.get(i).size()>=minSize && regs.get(i).size()<=maxSize
//					cnt++;
//					float C=0;        	// score
//					for (int aa=0; aa<regs.get(i).size(); aa++) {  // loop elements of the region
//						int xcoord = regs.get(i).get(aa)[1];
//						int ycoord = regs.get(i).get(aa)[0];
//						spots.set(xcoord, ycoord, (byte)255);
//						C += _critpoint_det[xcoord][ycoord];
//					}
//					C /= regs.get(i).size();   // score (calculate it here regardless of the type, says how much average fuzzy score was confident on the output)
//					//System.out.println("region elements:  " + regs.get(i).size() + " : " + C );
//				}
//			}
		}

		//System.out.println("loop regions.. " + minSize + " till " + maxSize);

        ArrayList<ArrayList<int[]>> regs = conn_reg.getConnectedRegions();
		ArrayList<float[]> vxy = new ArrayList<float[]>(); // will be filled up for each SALIENT branch

		// regs
		for (int i=0; i<regs.size(); i++) {

			if (regs.get(i).size()>=minSize && regs.get(i).size()<=maxSize) {

				float Cx=0, Cy=0; 	// centroid
				float C=0;        	// score
				float Cr; 			// radius

				CritpointRegion.RegionType Ctype;   	// type
				float[][] Cdirections;           		// directions

				for (int aa=0; aa<regs.get(i).size(); aa++) {  // loop elements of the region

					int xcoord = regs.get(i).get(aa)[1];
					int ycoord = regs.get(i).get(aa)[0];

					Cx += regs.get(i).get(aa)[1];
					Cy += regs.get(i).get(aa)[0];
					C += _critpoint_det[xcoord][ycoord];

				}

				Cx /= regs.get(i).size();
				Cy /= regs.get(i).size();   // centroid

				C /= regs.get(i).size();    // score (calculate it here regardless of the type, says how much average fuzzy score was confident on the output)

				Cr = det_radius; 		    // radius is not calculated but fixed wrt to the parameter

				// choose salient regions before determining their subtype and direction

				// second part (type, outward_directions) depends on which type of critical point we deal with
				Ctype = null;   			// values to add for this region

				// vxy list
				vxy.clear();  // list of local directions taken from the region

				for (int aa=0; aa<regs.get(i).size(); aa++) {  // loop elements of the region again

					int xcoord = regs.get(i).get(aa)[1];
					int ycoord = regs.get(i).get(aa)[0];
					int icoord = Profiler2D.xy2i[xcoord][ycoord];

					if (icoord!= -1) {  // peaks_i[icoord] is not null
						// element of the region is in foreground, take its thetas (if they exist)
						for (int peak_idx = 0; peak_idx < PeakExtractor2D.peaks_i[icoord].length; peak_idx++) {

							// iccord is the location index which is eqivalent to branch index
							// check if it exists and if it exists check whether the branch is on
							int curr_peak_i = PeakExtractor2D.peaks_i[icoord][peak_idx];
							boolean curr_peak_on = FuzzyDetector2D.branch_score[icoord][peak_idx]>output_membership_th;

							if (curr_peak_i!=-1 && curr_peak_i!=-2 && curr_peak_on) { // indexes of the spatial locations corresponding to peaks

								int peak_x = PeakExtractor2D.i2xy[curr_peak_i][0]; // PeakExtractor2D stores spatial location of the follow-up points
								int peak_y = PeakExtractor2D.i2xy[curr_peak_i][1];

								float[] unit_vxy = new float[]{peak_x-Cx, peak_y-Cy};
								float norm_vxy = (float) Math.sqrt(Math.pow(unit_vxy[0],2)+Math.pow(unit_vxy[1],2));
								unit_vxy[0] = unit_vxy[0] / norm_vxy;
								unit_vxy[1] = unit_vxy[1] / norm_vxy;

								vxy.add(unit_vxy);

							}

						}

					}

				}

				ArrayList<float[]> direction_clusters_vxy_count = ClusterDirections.run(vxy, alfa_deg);

				if (_choose_type == CritpointRegion.RegionType.END) {
					Ctype = CritpointRegion.RegionType.END; // this is known
				}
				else if (_choose_type == CritpointRegion.RegionType.BIF_CROSS) {
					// decide whether it is bif (3- clusters) or CROSS (4 clusters)
					if (direction_clusters_vxy_count.size()==4 &&  direction_clusters_vxy_count.get(3)[2]/direction_clusters_vxy_count.get(2)[2]>0.8) // if the last one is balanced towards smallest remaining
						Ctype = CritpointRegion.RegionType.CROSS;
					else
						Ctype = CritpointRegion.RegionType.BIF;
				}

				Cdirections = new float[direction_clusters_vxy_count.size()][2];  // will take 2 columns (3rd will be used to determine if it is crossing)
				for (int k=0; k<direction_clusters_vxy_count.size(); k++) {
					Cdirections[k][0] = direction_clusters_vxy_count.get(k)[0];
					Cdirections[k][1] = direction_clusters_vxy_count.get(k)[1];
				}

				switch (Ctype) {
					case END:
						region_list.add(new CritpointRegion(Ctype, Cx, Cy, Cr, C, Cdirections, 1)); // take one from the top only, no need to differentiate
						break;
					case BIF:
						region_list.add(new CritpointRegion(Ctype, Cx, Cy, Cr, C, Cdirections, 3));
						break;
					case CROSS:
						region_list.add(new CritpointRegion(Ctype, Cx, Cy, Cr, C, Cdirections, 4));
						break;
					default:
						break;
				}

			}
		}

	}

	private static void spatialScoreDistribution(float[][] critpoint2d, float[] critpoint1d, int[][] _i2xy) {

		for (int ll=0; ll<critpoint1d.length; ll++) {

			//if (critpoint1d[ll]!=null) {

			int x = _i2xy[ll][0];
			int y = _i2xy[ll][1];

			critpoint2d[x][y] = critpoint1d[ll];

			//}

		}
	}

	public static void print(int atX, int atY)
	{

		int atLoc = Masker2D.xy2i[atX][atY];

		IJ.log(String.format("/**** DETECTOR 2D TOOL (%5d, %5d) [%10d] ****/", atX, atY, atLoc));

		if (atLoc != -1) {  // first level check, intensity masker

			if (Delineator2D.delin2[atLoc]!=null) {  // second level check, profile masker

				String printout = "";

				printout += "\nDELINEATION INDEXES:\n";

				for (int ii=0; ii<Delineator2D.delin2[atLoc].length; ii++) {
					printout += ii+"\t->\t";
					for (int jj=0; jj<Delineator2D.delin2[atLoc][ii].length; jj++) {

						if (Delineator2D.delin2[atLoc][ii][jj]==-1) {
							printout += "BGRD";
						}
						else if (Delineator2D.delin2[atLoc][ii][jj]==-2) {
							printout += "NONE";
						}
						else {
							printout += IJ.d2s(Delineator2D.delin2[atLoc][ii][jj], 2);
						}

						if (jj==Delineator2D.delin2[atLoc][ii].length-1) printout += "\n";
						else printout += ",  ";
					}
				}

				printout += "\nREFINED LOCS:\n";
				if (Delineator2D.xy2[atLoc]!=null) {

					for (int b=0; b<Delineator2D.xy2[atLoc].length; b++) {

						printout += b+"\t->\t";

						if (Delineator2D.xy2[atLoc][b]!=null) {

							if (Delineator2D.xy2[atLoc][b][0] != null) {

								for (int l=0; l<Delineator2D.xy2[atLoc][b][0].length; l++) {

									printout += "("+IJ.d2s(Delineator2D.xy2[atLoc][b][0][l], 2)+", "+IJ.d2s(Delineator2D.xy2[atLoc][b][1][l], 2)+")";

									if (l==Delineator2D.xy2[atLoc][b][0].length-1) printout += "\n";
									else printout += ", ";
								}

							}
							else {

								printout += "NULL\n";

							}




						}
						else {
							printout += "NONE\n";
						}

					}

				}
				else {
					printout += "SKIPPED CALCULATING HERE (THERE WAS A THREAD POINTING TO BGRD)\n";
				}

				printout += "\nREFINED VECS:\n";
				if (Delineator2D.vxy2[atLoc]!=null) {

					for (int b=0; b<Delineator2D.vxy2[atLoc].length; b++) {

						printout += b+"\t->\t";

						if (Delineator2D.vxy2[atLoc][b]!=null) {

							if (Delineator2D.vxy2[atLoc][b][0] != null) {

								for (int l=0; l<Delineator2D.vxy2[atLoc][b][0].length; l++) {

									printout += "("+IJ.d2s(Delineator2D.vxy2[atLoc][b][0][l], 2)+", "+IJ.d2s(Delineator2D.vxy2[atLoc][b][1][l], 2)+")";

									if (l==Delineator2D.vxy2[atLoc][b][0].length-1) printout += "\n";
									else printout += ", ";
								}

							}
							else {

								printout += "NULL\n";

							}


						}
						else {
							printout += "NONE\n";
						}

					}

				}
				else {
					printout += "SKIPPED CALCULATING HERE (THERE WAS A THREAD POINTING TO BGRD)\n";
				}

				printout += "\nLHOOD STHNESS NCC\n";

				if (Delineator2D.smoothness[atLoc]!=null && Ncc2D.scores[atLoc]!=null) {
					for (int b=0; b<4; b++)
						printout += b + " -> " + IJ.d2s(PeakExtractor2D.peaks_lhood[atLoc][b], 1) + " " + IJ.d2s(Delineator2D.smoothness[atLoc][b], 1) + " " + IJ.d2s(Ncc2D.scores[atLoc][b], 1) + "\n";
				}
				else {
					printout += "NONE(profile was flat)\n";
				}

                printout += "\nEND <- NONE -> JUN\n";

                printout += "END->" + FuzzyDetector2D.endpoint_score[atLoc] + " JUN ->" + FuzzyDetector2D.junction_score[atLoc];

//                if (FuzzyDetector2D.endpoint_score[atLoc]!=null && Ncc2D.scores[atLoc]!=null) {
//
//                    for (int b=0; b<4; b++)
//                        printout += b + " -> " + IJ.d2s(PeakExtractor2D.peaks_lhood[atLoc][b], 1) + " " + IJ.d2s(Delineator2D.smoothness[atLoc][b], 1) + " " + IJ.d2s(Ncc2D.scores[atLoc][b], 1) + "\n";
//                }
//                else {
//                    printout += "NONE(profile was flat)\n";
//                }

				IJ.log(printout);

			}
			else {
				IJ.log("not enough variation in profile, profile masker said so");
			}

		}
		else {
			IJ.log("background point, intensity masker said so");
		}

	}

    private ImageStack getFuzzyLogicResponse(float ncc_1, float likelihood_1, float smoothness_1)
    {
        // emulate what the Fuzzy2D instance was giving for this input, use sample_fls initiated in the same way
        sample_fls.clearLog(); // to be empty for the next call
        sample_fls.verbose = true; // enable capturing the results of different stages of fuzzy logic detection
        float[] dummy = new float[3];
        float[] dummy1 = new float[4];
        sample_fls.critpointScore(ncc_1, likelihood_1, smoothness_1, dummy, dummy1);
        ImageStack is_out; // retrieve the cummulative log
        is_out = sample_fls.fls_steps.duplicate();
        sample_fls.clearLog(); // to be empty for the next call
        sample_fls.verbose = false;
        return is_out;
    }

    private ImageStack getFuzzyLogicResponse(float ncc_1, float likelihood_1, float smoothness_1,
                                             float ncc_2, float likelihood_2, float smoothness_2)
    {
        // emulate what the Fuzzy2D instance was giving for this input, use sample_fls initiated in the same way
        sample_fls.clearLog(); // to be empty for the next call
        sample_fls.verbose = true; // enable capturing the results of different stages of fuzzy logic detection
        float[] dummy = new float[3];
        float[] dummy1 = new float[4];
        sample_fls.critpointScore(
                ncc_1, likelihood_1, smoothness_1,
                ncc_2, likelihood_2, smoothness_2,
                dummy, dummy1);
        ImageStack is_out; // retrieve the cummulative log
        is_out = sample_fls.fls_steps.duplicate();
        sample_fls.clearLog(); // to be empty for the next call
        sample_fls.verbose = false;
        return is_out;
    }

    private ImageStack getFuzzyLogicResponse(float ncc_1, float likelihood_1, float smoothness_1,
                                             float ncc_2, float likelihood_2, float smoothness_2,
                                             float ncc_3, float likelihood_3, float smoothness_3)
    {
        // emulate what the Fuzzy2D instance was giving for this input, use sample_fls initiated in the same way
        sample_fls.clearLog(); // to be empty for the next call
        sample_fls.verbose = true; // enable capturing the results of different stages of fuzzy logic detection
        float[] dummy = new float[3];
        float[] dummy1 = new float[4];
        sample_fls.critpointScore(
                ncc_1, likelihood_1, smoothness_1,
                ncc_2, likelihood_2, smoothness_2,
                ncc_3, likelihood_3, smoothness_3,
                dummy, dummy1);
        ImageStack is_out; // retrieve the cummulative log
        is_out = sample_fls.fls_steps.duplicate();
        sample_fls.clearLog(); // to be empty for the next call
        sample_fls.verbose = false;
        return is_out;
    }

    private ImageStack getFuzzyLogicResponse(float ncc_1, float likelihood_1, float smoothness_1,
                                             float ncc_2, float likelihood_2, float smoothness_2,
                                             float ncc_3, float likelihood_3, float smoothness_3,
                                             float ncc_4, float likelihood_4, float smoothness_4)
    {
        // emulate what the Fuzzy2D instance was giving for this input, use sample_fls initiated in the same way
        sample_fls.clearLog(); // to be empty for the next call
        sample_fls.verbose = true; // enable capturing the results of different stages of fuzzy logic detection
        float[] dummy = new float[3];
        float[] dummy1 = new float[4];
        sample_fls.critpointScore(
                ncc_1, likelihood_1, smoothness_1,
                ncc_2, likelihood_2, smoothness_2,
                ncc_3, likelihood_3, smoothness_3,
                ncc_4, likelihood_4, smoothness_4,
                dummy, dummy1);
        ImageStack is_out; // retrieve the cummulative log
        is_out = sample_fls.fls_steps.duplicate();
        sample_fls.clearLog(); // to be empty for the next call
        sample_fls.verbose = false;
        return is_out;
    }

    public ImageStack getFuzzyLogicResponse(int atX, int atY)
    {

        ImageStack is_out = new ImageStack(528, 255);

        int atLoc = Masker2D.xy2i[atX][atY];

        if (atLoc != -1) {

            if (Delineator2D.smoothness[atLoc]!=null && Ncc2D.scores[atLoc]!=null) {

                int cnt = 0; // count branches that are existing
                float ncc_1=0, lhood_1=0, smthness_1=0;
                float ncc_2=0, lhood_2=0, smthness_2=0;
                float ncc_3=0, lhood_3=0, smthness_3=0;
                float ncc_4=0, lhood_4=0, smthness_4=0;

                for (int b=0; b<4; b++) {

                    float curr_ncc      = Ncc2D.scores[atLoc][b];
                    float curr_lhood    = PeakExtractor2D.peaks_lhood[atLoc][b];
                    float curr_smooth   = Delineator2D.smoothness[atLoc][b];

                    if (!Float.isNaN(curr_ncc) && !Float.isNaN(curr_lhood) && !Float.isNaN(curr_smooth)) {

                        cnt++;

                        if (cnt==1) {
                            ncc_1 = curr_ncc;
                            lhood_1 = curr_lhood;
                            smthness_1 = curr_smooth;
                        }
                        else if (cnt==2) {
                            ncc_2 = curr_ncc;
                            lhood_2 = curr_lhood;
                            smthness_2 = curr_smooth;
                        }
                        else if (cnt==3) {
                            ncc_3 = curr_ncc;
                            lhood_3 = curr_lhood;
                            smthness_3 = curr_smooth;
                        }
                        else if (cnt==4) {
                            ncc_4 = curr_ncc;
                            lhood_4 = curr_lhood;
                            smthness_4 = curr_smooth;
                        }


                    }

                }

                if (cnt==1) {
                    is_out = getFuzzyLogicResponse(ncc_1, lhood_1, smthness_1);
                }
                else if (cnt==2) {
                    is_out = getFuzzyLogicResponse(ncc_1, lhood_1, smthness_1, ncc_2, lhood_2, smthness_2);
                }
                else if (cnt==3) {
                    is_out = getFuzzyLogicResponse(ncc_1, lhood_1, smthness_1, ncc_2, lhood_2, smthness_2, ncc_3, lhood_3, smthness_3);
                }
                else if (cnt==4) {
                    is_out = getFuzzyLogicResponse(ncc_1, lhood_1, smthness_1, ncc_2, lhood_2, smthness_2, ncc_3, lhood_3, smthness_3, ncc_4, lhood_4, smthness_4);
                }

            }
            else {
                is_out.addSlice("NOTHING", new Plot("","","").getProcessor());
            }

        }
        else {
            is_out.addSlice("NOTHING", new Plot("","","").getProcessor());
        }

        return is_out;

    }

}