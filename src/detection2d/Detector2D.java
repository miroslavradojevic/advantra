package detection2d;

import aux.ClusterDirections;
import aux.ReadSWC;
import conn.Find_Connected_Regions;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.Line;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.io.FileSaver;
import ij.process.AutoThresholder;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;

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
	float			output_sigma;


	float 		output_membership_th;       	// based on k and output_sigma



	int         	M = 1;
	float       	minCos = -.5f;
	float			k = 4.0f;                   // hardcoded, defines the stricktness of the threshold for out membership function, k higher => means more strictness
	String          eval_string = "";





	String		output_dir_name; 		    	// parameter coded output folder name
	String 		output_log_name;
	PrintWriter logWriter = null;


	int         CPU_NR;



	// output
	static float[][] critpoint_det; 		// store 2d grid with detection scores (size of the image itself)
	ArrayList<CritpointRegion> detected_regions; // list with the detections

    public Detector2D(
							 ImagePlus 	ip_load,
							 float 		_s,
							 float[] 	_D,
							 float 		_sigma_ratio,
							 float 		_ncc_high,
							 float 		_ncc_low,
							 float 		_likelihood_high,
							 float 		_likelihood_low,
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

		critpoint_det = new float[ip_load.getWidth()][ip_load.getHeight()];
		detected_regions = new ArrayList<CritpointRegion>();

    }

    public void run()
	{

		long t1, t2;
		t1 = System.currentTimeMillis();

		for (int didx = 0; didx<D.length; didx++) { // loop scales D[]


			System.out.print("\nD = " + D[didx] + " : ");
			Sphere2D sph2d = new Sphere2D(D[didx], s, sigma_ratio);
			//sph2d.showSampling().show(); sph2d.showWeights().show();
			/********/






			System.out.print("Masker2D...");
			float new_masker_radius = 1.5f*sph2d.getOuterRadius();   	// important that it is outer radius of the sphere
			float new_masker_percentile = 30;                   		// used to have these two as argument but not necessary
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
			ImagePlus mask = new ImagePlus("mask", Masker2D.getMask());
			IJ.saveAs(mask, "Tiff", image_dir+image_name+".mask_values.tif");
			/********/







			System.out.print("Profiler2D...");
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





			System.out.print("ProfileSpread2D...");
			ProfileSpread2D.loadTemplate(Masker2D.i2xy, Masker2D.xy2i, Profiler2D.prof2, inimg_xy.length, inimg_xy[0].length);
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
			ImagePlus testip = new ImagePlus("", ProfileSpread2D.getMask());
			IJ.saveAs(testip, "Tiff", image_dir+image_name+".mask_profile.tif");
			System.out.println("" + IJ.d2s( (ProfileSpread2D.getNrCritpointCandidates()*100f)/(inimg_xy.length * inimg_xy[0].length), 2) + " % candidates");
			/********/











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
			/********/







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
//			new ImagePlus("", Delineator2D.getSmoothnessDistribution(64)).show();
			/********/



			System.out.println("DID IT!!!");
			if (true) return;





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




			System.out.print("Append regions... ");
			// bif-crs clustering parameters
			float ang_deg = 20f;

			// imagine a square scaled wrt current D[idx] used for detection
			int			min_region_size = (int) Math.round(0.25f*0.25f*Math.pow(D[didx],2));  // smallest connected region area to be valid critical point
			int			max_region_size = (int) Math.round(2.0f*2.0f*Math.pow(D[didx],2));  // largest connected region to be valid critical point
			float 		region_radius 	= D[didx];
			float sensitivity = 0.7f;

//			System.out.println("diameter " + region_radius + " min = " + min_region_size + " max = " + max_region_size);

			FeatureExtractor2D.exportCritpointScores(critpoint_det, 0); // store critpoint_score2 list into 2d array critpoint_det

			getCritpointRegions(
																					critpoint_det,
																					min_region_size,
																					max_region_size,
									   												region_radius,	          // the same radius will be assigned at each scale
																					CritpointRegion.RegionType.END,
																					ang_deg,
//																					min_clust_cnt,
																					detected_regions
																					);

			FeatureExtractor2D.exportCritpointScores(critpoint_det, 2); // store critpoint_score2 list into 2d array critpoint_det
			getCritpointRegions(
																					critpoint_det,
																					min_region_size,
																					max_region_size,
									   												region_radius,
																					CritpointRegion.RegionType.BIF_CROSS,
																					ang_deg,
																					detected_regions
																					);

																					*/

			System.out.println("finished scale " + D[didx] + " out of " + D.length );
		}

		t2 = System.currentTimeMillis();
		IJ.log("done. " + ((t2 - t1) / 1000f) + "sec.");

		// loop here the list
//		IJ.log(bif_list.size()+" bifurcations ");
//		IJ.log(end_list.size()+" endpoints ");

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
//			ovroi.setFillColor(region_color);
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

		return ov;

	}

	// appends to the list of CritpointRegion, using connected components to optimize region segmentation
	private void getCritpointRegions(
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

		byte[] t = new byte[w*h];
		for (int x=0; x<w; x++) {
			for (int y=0; y<h; y++) {
				float curr_on = _critpoint_det[x][y];
				if (curr_on>=output_membership_th) { // threshold score on critpoint detection
					t[y*w+x] = (byte) 255;
				}
			}
		}

		// take detections (binary image), find connected regions
		ByteProcessor bp = new ByteProcessor(w, h, t);
		Find_Connected_Regions conn_reg = new Find_Connected_Regions(new ImagePlus("", bp), true);
		conn_reg.run("");

		//conn_reg.showLabels().show();

		FileSaver fs = new FileSaver(conn_reg.showLabels());
		if (_choose_type == CritpointRegion.RegionType.END) fs.saveAsTiff(image_dir+"conn_ends"+det_radius+".tif");
		if (_choose_type == CritpointRegion.RegionType.BIF_CROSS) fs.saveAsTiff(image_dir+"conn_jun"+det_radius+".tif");


		ArrayList<ArrayList<int[]>> regs = conn_reg.getConnectedRegions();

		// just thresholding and connected components will give a lot of noisy regions



		ArrayList<float[]> vxy = new ArrayList<float[]>(); // will be filled up for each SALIENT region

		if (_choose_type== CritpointRegion.RegionType.END) {

			// export binary image
			ByteProcessor spots = new ByteProcessor(_critpoint_det.length, _critpoint_det[0].length);


			System.out.println("Endpoints... what happened? total " + regs.size() + "regions" + minSize + " <> " +maxSize);

			int cnt = 0;
			for (int i=0; i<regs.size(); i++) {

				if (true) { // regs.get(i).size()>=minSize && regs.get(i).size()<=maxSize

					cnt++;
					float C=0;        	// score

					for (int aa=0; aa<regs.get(i).size(); aa++) {  // loop elements of the region

						int xcoord = regs.get(i).get(aa)[1];
						int ycoord = regs.get(i).get(aa)[0];
						spots.set(xcoord, ycoord, (byte)255);
						C += _critpoint_det[xcoord][ycoord];

					}

					C /= regs.get(i).size();   // score (calculate it here regardless of the type, says how much average fuzzy score was confident on the output)

					//System.out.println("region elements:  " + regs.get(i).size() + " : " + C );

				}
			}

			//IJ.saveAs(new ImagePlus("", spots), "Tiff", "/home/miroslav/Desktop/ends.tif");



		}

		//System.out.println("loop regions.. " + minSize + " till " + maxSize);

		// regs
		for (int i=0; i<regs.size(); i++) {

//			if (i==1) System.out.println("tried region 1 " + regs.get(i).size());

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
				Cy /= regs.get(i).size();  // centroid

				C /= regs.get(i).size();   // score (calculate it here regardless of the type, says how much average fuzzy score was confident on the output)

				Cr = det_radius; 				// radius is not calculated but fixed wrt to the parameter

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
							boolean curr_peak_on = FeatureExtractor2D.streamline_score2[icoord][2][peak_idx]>output_membership_th;

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
						region_list.add(new CritpointRegion(Ctype, Cx, Cy, Cr, C, Cdirections, 1)); // take one from the top
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

	public static void print(int atX, int atY)
	{

//		System.out.println(Arrays.toString(peaks_i[atLoc]));

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



				printout += "b \t -> \t\t STH. \t\t LHOOD \t\t NCC.\n";

				if (Delineator2D.smoothness[atLoc]!=null) {

					for (int b=0; b<Delineator2D.smoothness[atLoc].length; b++) {
						printout += (b+1) +"\t->\t\t";
						if (!Float.isNaN(Delineator2D.smoothness[atLoc][b])) {
							printout += "" + IJ.d2s(Delineator2D.smoothness[atLoc][b], 2)+
												" \t "+IJ.d2s(Delineator2D.smoothness[atLoc][b], 2)+
												" \t "+IJ.d2s(Delineator2D.smoothness[atLoc][b], 2)
//												+
//												"  \t \t  OFF " + IJ.d2s(1, 2)
//												+
//												"  \t \t   ON "+  IJ.d2s(1, 2)
												+"\n"
							;
						}
						else {
							printout += "NaN\n";
						}
					}
				}
				else {
					printout += "NONE(profile was flat)\n";
				}


				/*
				if (ncc2[atLoc]!=null) {
					for (int b=0; b<ncc2[atLoc].length; b++) {
						printout += (b+1) +"\t->\t"; //+ IJ.d2s(ratio2[atLoc][ii], 2) + "\n"
						if (!Float.isNaN(ncc2[atLoc][b])) {
							printout += "NCC " + IJ.d2s(ncc2[atLoc][b], 2)+"  \t LHOOD "+IJ.d2s(lhoods2[atLoc][b], 2)+"  \t \t  OFF "
												+IJ.d2s(streamline_score2[atLoc][0][b], 2)+"  \t  ON "+IJ.d2s(streamline_score2[atLoc][2][b], 2)+"\n";
						}
						else {
							printout += "NONE\n";
						}
					}
				}
				else {
					printout += "NONE\n";
				}
				*/

//				printout += "\nBRANCH STRENGTH \n...todo...\n";

			/*printout += "\nCRITPOINT \n";
			if (critpoint_score2[atLoc]!=null) {
				printout += "END "+critpoint_score2[atLoc][0]+", NONE "+critpoint_score2[atLoc][1]+", BIF "+critpoint_score2[atLoc][2];
			}
			else {
				printout += "NONE\n";
			}*/

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

}