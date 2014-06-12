package detection2d;

import aux.ClusterDirections;
import aux.ReadSWC;
import conn.Find_Connected_Regions;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Line;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.process.ByteProcessor;
import sun.nio.cs.StreamDecoder;

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
	static float[]  D;
	int         	L;
	float       	ncc_high;
	float       	ncc_low;
	float 			likelihood_high;
	float 			likelihood_low;
	float			output_sigma;


	float 		output_membership_th;       	// based on k and output_sigma



	int         	M = 1;
	float       	minCos = -.5f;
	float       	sampling_crosswise = .5f;
	float			k = 0.5f;                   // hardcoded, defines the stricktness of the threshold for out membership function, k higher => means more strictness
	String          eval_string = "";



	String		output_dir_name; 		    	// parameter coded output folder name
	String 		output_log_name;
	PrintWriter logWriter = null;


	int         CPU_NR;



	// output
	static float[][] critpoint_det; 		// store 2d grid with detection scores (size of the image itself)
	ArrayList<CritpointRegion> detected_regions; // list with the detections

    public Detector2D(
							 ImagePlus ip_load,
							 float _s,
							 float[] _D,
							 int _L,
							 float _ncc_high,
							 float _ncc_low,
							 float _likelihood_high,
							 float _likelihood_low,
							 float _output_sigma
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

		image_gndtth_endpoints = image_name + ".end";
		image_gndtth_bifurcations = image_name + ".bif";
		image_gndtth_crosssections = image_name + ".crs";

		s = _s;
		D = _D; // just pointer assignement, no memory allocated
		L = _L;
		ncc_high = _ncc_high;
		ncc_low = _ncc_low;
		likelihood_high = _likelihood_high;
		likelihood_low = _likelihood_low;
		output_sigma = _output_sigma;

//		min_D = Float.POSITIVE_INFINITY;
//		max_D = Float.NEGATIVE_INFINITY;
		String		Dlist = ""; // for the out directory
		//			if (_D[i]<min_D) min_D = _D[i];
//			if (_D[i]>max_D) max_D = _D[i];
		for (int i=0; i<_D.length; i++) Dlist += IJ.d2s(_D[i], 1) + ((i == _D.length - 1) ? "" : ":");

//		min_region_size = Math.round(min_D);
//		max_region_size = Math.round(min_D * min_D);


		output_membership_th = (float) Math.exp(-(0.5f*0.5f)/(2*output_sigma*output_sigma)); // 0.5 is middle of the output range
//		System.out.println("border " + output_membership_th);
		output_membership_th = (float) (1 - Math.pow(output_membership_th,k) * (1-output_membership_th)); // h1 = 1 - h^k * (1-h)


//		System.out.println("th was "+output_membership_th);

		output_dir_name = image_dir+String.format(       "det.s.Dlist.M.L.nncH.nccL.lhoodH.lhoodL.outSig_%.1f_%s_%d_%d_%.2f_%.2f_%.2f_%.2f_%.2f",
														 this.s,
														 Dlist,
														 M, L,
														 ncc_high,
														 ncc_low,
														 likelihood_high,
														 likelihood_low,
														 output_sigma);
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

		for (int didx = 0; didx<D.length; didx++) {





			Sphere2D sph2d = new Sphere2D(D[didx], s); //sph2d.showSampling().show(); sph2d.showWeights().show();




			float new_masker_radius = 1.5f*sph2d.getOuterRadius();   // important that it is outer radius of the sphere
			float new_masker_percentile = 30;                   // used to have these two as argument but not necessary
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
			//new ImagePlus("MASK,r="+IJ.d2s(new_masker_radius,2)+",per="+IJ.d2s(new_masker_percentile,0), Masker2D.getMask()).show();
			/********/









//			IJ.log("calculating profiles...");
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
//			IJ.log("extracting peaks...");
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
//			IJ.log("fitting model... extracting features... fuzzy detection...");
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
												   "MEDIAN",
												   ncc_high,
												   ncc_low,
												   likelihood_high,
												   likelihood_low,
												   output_sigma,
												   L,
												   sampling_crosswise
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


			// extract the regions to the list






			/********/
//			IJ.log("append ends...");

			// bif-crs clustering parameters
			float ang_deg = 20f;

			// imagine a square scaled wrt current D[idx] used for detection
			int			min_region_size = (int) Math.round(0.5f*0.5f*Math.pow(D[didx],2));  // smallest connected region area to be valid critical point
			int			max_region_size = (int) Math.round(0.9f*0.9f*Math.pow(D[didx],2));  // largest connected region to be valid critical point
			float 		region_radius 	= D[didx];

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

//			IJ.log("append bifs-crosses...");
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
			ovroi.setFillColor(region_color);
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
		//FileSaver fs = new FileSaver(conn_reg.showLabels());
		//fs.saveAsTiff(image_dir+"conn_ends.tif");
		ArrayList<ArrayList<int[]>> regs = conn_reg.getConnectedRegions();

		ArrayList<float[]> vxy = new ArrayList<float[]>(); // will be filled up for each region

		// regs are converted into cls
		for (int i=0; i<regs.size(); i++) {
			if (regs.get(i).size()>=minSize && regs.get(i).size()<=maxSize) {

				float Cx=0, Cy=0; 	// centroid
				float C=0;        	// score
				float Cr; 			// radius

				CritpointRegion.RegionType Ctype;   // type
				float[][] Cdirections;           // directions

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

				// second part (type, outward_directions) depends on which type of critical point we deal with
				Ctype = null;   			// values to add for this region
				Cdirections = null;

				/*
					CASE ENDPOINT
				 */

				if (_choose_type == CritpointRegion.RegionType.END) { // it is endpoint, no need to further reclassify, there is only one direction

					// loop through region members to extract the directions

					vxy.clear();  // list of local directions taken from the region

					for (int aa=0; aa<regs.get(i).size(); aa++) {  // loop elements of the region again

						int xcoord = regs.get(i).get(aa)[1];
						int ycoord = regs.get(i).get(aa)[0];
						int icoord = Profiler2D.xy2i[xcoord][ycoord];

						if (icoord!= -1) {  // peaks_i[icoord] is not null
							// element of the region is in foreground, take its thetas (if they exist)
							for (int peak_idx = 0; peak_idx < PeakExtractor2D.peaks_i[icoord].length; peak_idx++) {

								// iccord is the location index
								// which is eqivalent to branch index
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



					// use mean-shift to find the main direction out of vxy list
//					int min_cluster_cnt = 3;// to avoid outliers (int) Math.round(0.8*regs.get(i).size());
					int[] cls_counts = ClusterDirections.run(vxy, alfa_deg, Cdirections); // gives out descending sorted clusters

					System.out.println(Cdirections.length+" it was");

//					Cdirections = new float[1][2];
//					Cdirections[0][0] = cls.get(0)[0];
//					Cdirections[0][1] = cls.get(0)[1];

					Ctype = CritpointRegion.RegionType.END; // this is known

				}

				/*
					CASE BIFURCATION OR CROSSPOINT (2D)
				 */

				else if (_choose_type == CritpointRegion.RegionType.BIF_CROSS) {
					// (unlike endpoint detection) it still has to be defined whether it is bifurcation of cross section based on the scores

					// loop through region members to extract whether it is bif or cross and to extract the directions
					vxy.clear();
					for (int aa=0; aa<regs.get(i).size(); aa++) {  // loop elements of the region again, collect the directions

						int xcoord = regs.get(i).get(aa)[1];
						int ycoord = regs.get(i).get(aa)[0];
						int icoord = Profiler2D.xy2i[xcoord][ycoord];

						if (icoord!= -1) {  // peaks_i[icoord] is not null
							// element of the region is in foreground, take its thetas (if they exist) and put them together
							for (int peak_idx = 0; peak_idx < PeakExtractor2D.peaks_i[icoord].length; peak_idx++) {

								// iccord is the location index
								// which is eqivalent to branch index
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




					// use mean-shift to find the main direction out of vxy list
//					int min_cluster_cnt = 3;// to avoid outliers (int) Math.round(0.8*regs.get(i).size());

					int[] cls_counts = ClusterDirections.run(vxy, alfa_deg, Cdirections); // list of size-descending sorted clusters for directions (vx, vy)
					// Cdirections are stored now

//					Cdirections = new float[cls.size()][2];
//					for (int k=0; k<cls.size(); k++) {
//						Cdirections[k][0] = cls.get(k)[0];
//						Cdirections[k][1] = cls.get(k)[1];
//					}

					// decide whether it is bif (3 clusters) or CROSS (4 clusters)
					if (cls_counts.length==4 &&  (float)cls_counts[3]/(float)cls_counts[2]>0.8 ) // if the last one is balanced towards smallest remaining
						Ctype = CritpointRegion.RegionType.CROSS;
					else
						Ctype = CritpointRegion.RegionType.BIF;


				}

				/*
					append it here  and prune directions
				 */



				System.out.println("----");
				System.out.println(Ctype);
				for (int ii = 0; ii < Cdirections.length; ii++) {
					System.out.println(""+ Arrays.toString(Cdirections[ii]));
				}
				System.out.println("----");

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

}