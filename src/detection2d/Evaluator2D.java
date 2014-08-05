package detection2d;

import aux.ReadDET;
import aux.ReadSWC;
import aux.Tools;
import ij.*;
import ij.gui.GenericDialog;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;

import java.awt.*;
import java.io.*;
import java.util.ArrayList;

/**
 * Created by miroslav on 30-6-14.
 * the task is to compare one detection (.tif with Overlay, or .det) and ground truth (.swc)
 * and store the log string in eval.csv n the same folder as detection and print it out on IJ output
 */
public class Evaluator2D implements PlugIn {

	int tp_BIF, fp_BIF, fn_BIF;
	float p_BIF, r_BIF;

	// convention - overlay colors
	Color bif_color = Color.RED;
	int   bif_type = 3;

	Color end_color = Color.YELLOW;
	int 	end_type = 1;

	Color crs_color = Color.GREEN;
	int 	crs_type = 4;

	Color ignore_color = new Color(0, 0, 0, 0.5f); // it's not used, doesn't matter, Annotationer needs it only
	int ignore_type = 7;

	int tp_END, fp_END, fn_END;
	float p_END, r_END;

	int tp_CRS, fp_CRS, fn_CRS;
	float p_CRS, r_CRS;

	// BIF and CRS together = JUN
	int tp_JUN, fp_JUN, fn_JUN;
	float p_JUN, r_JUN;

	Overlay			eval_overlay_jun = new Overlay();  // visualize the detections
	Overlay			eval_overlay_end = new Overlay();

	public void run(String s) { // regular PlugIn run() method

        String			det_file_path;     		// either .tif or .det

		boolean			show = false;

		// load the image with detections through the menu
		String in_folder = Prefs.get("id.folder", System.getProperty("user.home"));
		OpenDialog.setDefaultDirectory(in_folder);
		OpenDialog dc = new OpenDialog("Select detection file");
		in_folder = dc.getDirectory();
		det_file_path = dc.getPath();
		if (det_file_path==null) return;
		Prefs.set("id.folder", in_folder);

		String output_dir_name = in_folder;
		String output_log_name = output_dir_name + "eval.csv"; // output (file append) will be stored in the same folder as the detection file

		String 			gndtth_path;			// ground truth swc file with critical points to evaluate, same name as input image
		String			gndtth_tag;

		if (Macro.getOptions()==null) {

			GenericDialog gd = new GenericDialog("GROUND TRUTH?");
			gd.addStringField("gndtth_path (.swc)", 	new File(output_dir_name).getParent(), 60);
			gd.addStringField("gndtth_tag", 	"LABEL", 50);
			gd.addCheckbox("show",  false);
			gd.showDialog();
			if (gd.wasCanceled()) return;
			gndtth_path	= gd.getNextString();
			gndtth_tag	= gd.getNextString();
			show = gd.getNextBoolean();

		}
		else {
			gndtth_path = Macro.getValue(Macro.getOptions(), "gndtth_path", "ground_truth_path");
			gndtth_tag 	= Macro.getValue(Macro.getOptions(), "gndtth_tag", 	"ground_truth_tag");
		}

		/* read ground truth */

		File f_gndtth = new File(gndtth_path);

		if (!f_gndtth.exists()) {
			System.out.println("file does not exist");
			return;
		}

		if (!Tools.getFileExtension(gndtth_path).equals("swc")) {
			System.out.println("file needs to be .swc");
			return;
		}

		// structure to be filled up reading from the ground truth swc: list of <x,y,r>
		ArrayList<float[]> gndtth_juns = new ArrayList<float[]>();
		ArrayList<float[]> gndtth_bifs = new ArrayList<float[]>();
		ArrayList<float[]> gndtth_ends = new ArrayList<float[]>();
		ArrayList<float[]> gndtth_crss = new ArrayList<float[]>();
		ArrayList<float[]> gndtth_ignr = new ArrayList<float[]>();

        ReadSWC reader = new ReadSWC(gndtth_path);

		for (int b=0; b<reader.nodes.size(); b++) {

			int typ = Math.round(reader.nodes.get(b)[reader.TYPE]);
			float bx = reader.nodes.get(b)[reader.XCOORD];
			float by = reader.nodes.get(b)[reader.YCOORD];
			float br = reader.nodes.get(b)[reader.RADIUS];

			if (typ==end_type) {
				gndtth_ends.add(new float[]{bx, by, br});
			}
			else if (typ==bif_type) {
				gndtth_bifs.add(new float[]{bx, by, br});
				gndtth_juns.add(new float[]{bx, by, br});
			}
			else if (typ==crs_type) {
				gndtth_crss.add(new float[]{bx, by, br});
				gndtth_juns.add(new float[]{bx, by, br});
			}
			else if (typ==ignore_type) {
				gndtth_ignr.add(new float[]{bx, by, br});
			}

		}

		System.out.println(
								  "gndtth loaded -> " + (gndtth_juns.size() + gndtth_ends.size()) + " elements: " +
										  gndtth_juns.size() + " JUNS, " +
										  gndtth_ends.size() + " ENDS, " +
										  gndtth_bifs.size() + " BIFS, " +
										  gndtth_crss.size() + " CRSS, " +
										  gndtth_ignr.size() + " IGNR."
		);

		/* read the detections (use gndtth ignore regions to exclude those that we don't wish to evaluate) */
		// lists to be filled up by reading from the .det csv file <x,y,r>
		ArrayList<float[]> det_bifs = new ArrayList<float[]>();
		ArrayList<float[]> det_juns = new ArrayList<float[]>();
		ArrayList<float[]> det_ends = new ArrayList<float[]>();
		ArrayList<float[]> det_crss = new ArrayList<float[]>();

		// append to the list in different scenarios
		if (Tools.getFileExtension(det_file_path).equals("det")) {

			// read csv file (.det) at specified location
			ReadDET det_reader = new ReadDET(det_file_path);

            // print loaded values
            //System.out.println("****************************");
            //det_reader.print();
            //System.out.println("****************************");

			for (int i = 0; i < det_reader.x.size(); i++) { // loop read detections

				float reg_x = det_reader.x.get(i);
				float reg_y = det_reader.y.get(i);
				float reg_r = det_reader.r.get(i);

				// before adding it - check if it is not in some of the ignore regions of the ground truth
				boolean add_it = true;

				for (int j = 0; j < gndtth_ignr.size(); j++) {

					float ignr_x = gndtth_ignr.get(j)[0];
					float ignr_y = gndtth_ignr.get(j)[1];
					float ignr_r = gndtth_ignr.get(j)[2];

					if (Math.pow(reg_x - ignr_x, 2) + Math.pow(reg_y - ignr_y, 2) <= Math.pow(reg_r + ignr_r, 2))
						add_it = false;
				}

				if (!add_it || reg_r<=0) continue; // it overlapped with ignore region

				String region_type = det_reader.t.get(i);

				if (region_type.equals("BIF")) {
					det_bifs.add(new float[]{reg_x, reg_y, reg_r});
					det_juns.add(new float[]{reg_x, reg_y, reg_r});
				} else if (region_type.equals("END")) {
					det_ends.add(new float[]{reg_x, reg_y, reg_r});
				} else if (region_type.equals("CROSS")) {
					det_crss.add(new float[]{reg_x, reg_y, reg_r});
					det_juns.add(new float[]{reg_x, reg_y, reg_r});
				}

			}

			System.out.println("tif detection loaded -> " +
									   det_juns.size() + " JUNS, " +
									   det_ends.size() + " ENDS, " +
									   det_bifs.size() + " BIFS, " +
									   det_crss.size() + " CRSS. "
			);

		}
		else {
			System.out.println("Detection file extension not recognized.");
			return;
		}

		/* evaluation */

		doBifEvaluation(det_bifs, gndtth_bifs);
		doJunEvaluation(det_juns, gndtth_juns);
		doEndEvaluation(det_ends, gndtth_ends);
		doCrsEvaluation(det_crss, gndtth_crss);

		/* store output */

		PrintWriter 	logWriter = null;
		String legend = String.format("%10s,%10s,%6s,%6s,%6s,%6s,%6s,%6s,%6s,%6s,%6s,%6s,%6s,%6s,%6s,%6s,%6s,%6s,%6s,%6s,%6s,%6s",
											 "NAME", "ANNOT_TAG",
											 "TP_JUN", "FP_JUN", "FN_JUN", "P_JUN", "R_JUN",
											 "TP_BIF", "FP_BIF", "FN_BIF", "P_BIF", "R_BIF",
											 "TP_CRS", "FP_CRS", "FN_CRS", "P_CRS", "R_CRS",
											 "TP_END", "FP_END", "FN_END", "P_END", "R_END"
											 );

		File f = new File(output_log_name);
		if (!f.exists()) {
			try {
				logWriter = new PrintWriter(output_log_name);
				logWriter.println(legend);
				logWriter.close();
			} catch (FileNotFoundException ex) {}
		}
		// if it exists already in the folder, just prepare to append on the existing file
		try {
			logWriter = new PrintWriter(new BufferedWriter(new FileWriter(output_log_name, true)));
		} catch (IOException e) {}

		String eval = String.format("%10s,%10s,%6d,%6d,%6d,%6.2f,%6.2f,%6d,%6d,%6d,%6.2f,%6.2f,%6d,%6d,%6d,%6.2f,%6.2f,%6d,%6d,%6d,%6.2f,%6.2f",
										   Tools.getFileName(det_file_path),
										   gndtth_tag,
										   tp_JUN, fp_JUN, fn_JUN, p_JUN, r_JUN,
										   tp_BIF, fp_BIF, fn_BIF, p_BIF, r_BIF,
										   tp_CRS, fp_CRS, fn_CRS, p_CRS, r_CRS,
										   tp_END, fp_END, fn_END, p_END, r_END
										   );


		logWriter.println(eval);
		logWriter.close();

		System.out.println(legend);
		System.out.println(eval);

		if (show) { // show from GenericDialog, this is disabled in headless (because it makes sense only with graphic user interface)

			IJ.open();
			ImagePlus eval_image_jun = IJ.getImage();
			ImagePlus eval_image_end = eval_image_jun.duplicate();

			eval_image_jun.setTitle("*JUN* EVALUATION");
			eval_image_jun.setOverlay(eval_overlay_jun);
			eval_image_jun.show();

			eval_image_end.setTitle("*END* EVALUATION");
			eval_image_end.setOverlay(eval_overlay_end);
			eval_image_end.show();

		}

	}

	public void doBifEvaluation(ArrayList<float[]> _det_bif, ArrayList<float[]> _gndtth_bifs)
	{

		tp_BIF = 0; fp_BIF = 0; fn_BIF = 0;

		boolean[] annots = new boolean[_gndtth_bifs.size()];  // necessary for fn calculation

		// loop all detected BIF regions
		for (int a=0; a<_det_bif.size(); a++) {

				boolean found = false;

				float ax = _det_bif.get(a)[0];
				float ay = _det_bif.get(a)[1];
				float ar = _det_bif.get(a)[2];

				for (int b=0; b<_gndtth_bifs.size(); b++) {

					float bx = _gndtth_bifs.get(b)[0];
					float by = _gndtth_bifs.get(b)[1];
					float br = _gndtth_bifs.get(b)[2];

					if (circlesOverlap(ax, ay, ar, bx, by, br)) {// two roi circles overlap
						found = true;
						tp_BIF++;
						annots[b] = true;
					}
				}

				if (!found) fp_BIF++;  // detected but was not in the list of annotated ones

		}

		for (int a=0; a<annots.length; a++) if (!annots[a]) fn_BIF++;

		p_BIF = (tp_BIF+fp_BIF>0)? tp_BIF/(float)(tp_BIF+fp_BIF) : 0 ;
		r_BIF = (tp_BIF+fn_BIF>0)? tp_BIF/(float)(tp_BIF+fn_BIF) : 0 ;

	}

	public void doJunEvaluation(ArrayList<float[]> _det_jun, ArrayList<float[]> _gndtth_juns)
	{

		tp_JUN = 0; fp_JUN = 0; fn_JUN = 0;
		boolean[] annots = new boolean[_gndtth_juns.size()];  // necessary for fn calculation

		// loop all detected JUN regions
		for (int a=0; a<_det_jun.size(); a++) {

			boolean found = false;

			float ax = _det_jun.get(a)[0];
			float ay = _det_jun.get(a)[1];
			float ar = _det_jun.get(a)[2];

			for (int b=0; b<_gndtth_juns.size(); b++) {

				float bx = _gndtth_juns.get(b)[0];
				float by = _gndtth_juns.get(b)[1];
				float br = _gndtth_juns.get(b)[2];

				if (circlesOverlap(ax, ay, ar, bx, by, br, 1f)) {// two roi circles overlap, margin=1 due to the rounding
					found = true;
					tp_JUN++;

					// add tproi jun
					OvalRoi tproi = new OvalRoi(ax-ar+.5f, ay-ar+.5f, 2*ar, 2*ar);
					tproi.setFillColor(Color.RED);
					eval_overlay_jun.add(tproi);

					annots[b] = true;
				}
			}

			if (!found) {   // detected but was not in the list of annotated ones
				fp_JUN++;

				// add fproi jun
				OvalRoi fproi = new OvalRoi(ax-ar+.5f, ay-ar+.5f, 2*ar, 2*ar);
				fproi.setStrokeColor(Color.CYAN);
				eval_overlay_jun.add(fproi);
			}

		}

		for (int a=0; a<annots.length; a++)
			if (!annots[a]) {

				// this is actually only for visualization, otherwise bx,by,br are not necessary for the fn_JUN number itself
				float bx = _gndtth_juns.get(a)[0];
				float by = _gndtth_juns.get(a)[1];
				float br = _gndtth_juns.get(a)[2];

				fn_JUN++;

				// add fnroi jun
				OvalRoi fnroi = new OvalRoi(bx-br+.5f, by-br+.5f, 2*br, 2*br);
				fnroi.setFillColor(Color.GREEN);
				eval_overlay_jun.add(fnroi);

			}


		p_JUN = (tp_JUN+fp_JUN>0)? tp_JUN/(float)(tp_JUN+fp_JUN) : 0 ;
		r_JUN = (tp_JUN+fn_JUN>0)? tp_JUN/(float)(tp_JUN+fn_JUN) : 0 ;

	}

	public void doEndEvaluation(ArrayList<float[]> _det_end, ArrayList<float[]> _gndtth_ends)
	{

		tp_END = 0; fp_END = 0; fn_END = 0;
		boolean[] annots = new boolean[_gndtth_ends.size()];  // necessary for fn calculation
		// loop all detected END regions
		for (int a=0; a<_det_end.size(); a++) {

			boolean found = false;

			float ax = _det_end.get(a)[0];
			float ay = _det_end.get(a)[1];
			float ar = _det_end.get(a)[2];

			for (int b=0; b<_gndtth_ends.size(); b++) {

				float bx = _gndtth_ends.get(b)[0];
				float by = _gndtth_ends.get(b)[1];
				float br = _gndtth_ends.get(b)[2];

				if (circlesOverlap(ax, ay, ar, bx, by, br)) {// two roi circles overlap
					found = true;
					tp_END++;

					// add tproi end
					OvalRoi tproi = new OvalRoi(ax-ar+.5f, ay-ar+.5f, 2*ar, 2*ar);
					tproi.setFillColor(Color.RED);
					eval_overlay_end.add(tproi);

					annots[b] = true;
				}
			}

			if (!found) {
				fp_END++;  // detected but was not in the list of annotated ones

				// add fproi end
				OvalRoi fproi = new OvalRoi(ax-ar+.5f, ay-ar+.5f, 2*ar, 2*ar);
				fproi.setStrokeColor(Color.CYAN);
				eval_overlay_end.add(fproi);

			}

		}

		for (int a=0; a<annots.length; a++)
			if (!annots[a]){

			// this is actually only for visualization, otherwise bx,by,br are not necessary for the fn_END number itself
			float bx = _gndtth_ends.get(a)[0];
			float by = _gndtth_ends.get(a)[1];
			float br = _gndtth_ends.get(a)[2];

			fn_END++;

			// add fnroi end
			OvalRoi fnroi = new OvalRoi(bx-br+.5f, by-br+.5f, 2*br, 2*br);
			fnroi.setFillColor(Color.GREEN);
			eval_overlay_end.add(fnroi);
		}

		p_END = (tp_END+fp_END>0)? tp_END/(float)(tp_END+fp_END) : 0 ;
		r_END = (tp_END+fn_END>0)? tp_END/(float)(tp_END+fn_END) : 0 ;

	}

	public void doCrsEvaluation(ArrayList<float[]> _det_crs, ArrayList<float[]> _gndtth_crss)
	{

		tp_CRS = 0; fp_CRS = 0; fn_CRS = 0;
		boolean[] annots = new boolean[_gndtth_crss.size()];  // necessary for fn calculation
		// loop all detected CRS regions
		for (int a=0; a<_det_crs.size(); a++) {

			boolean found = false;

			float ax = _det_crs.get(a)[0];
			float ay = _det_crs.get(a)[1];
			float ar = _det_crs.get(a)[2];

			for (int b=0; b<_gndtth_crss.size(); b++) {

				float bx = _gndtth_crss.get(b)[0];
				float by = _gndtth_crss.get(b)[1];
				float br = _gndtth_crss.get(b)[2];

				if (circlesOverlap(ax, ay, ar, bx, by, br)) {// two roi circles overlap
					found = true;
					tp_CRS++;
					annots[b] = true;
				}
			}

			if (!found) fp_CRS++;  // detected but was not in the list of annotated ones

		}

		for (int a=0; a<annots.length; a++) if (!annots[a]) fn_CRS++;

		p_CRS = (tp_CRS+fp_CRS>0)? tp_CRS/(float)(tp_CRS+fp_CRS) : 0 ;
		r_CRS = (tp_CRS+fn_CRS>0)? tp_CRS/(float)(tp_CRS+fn_CRS) : 0 ;

	}

	private boolean circlesOverlap(float x1, float y1, float r1, float x2, float y2, float r2)
	{
		return Math.pow(x1-x2,2)+Math.pow(y1-y2,2) <= Math.pow(r1+r2,2);
	}

	private boolean circlesOverlap(float x1, float y1, float r1, float x2, float y2, float r2, float margin)
	{
		return Math.pow(x1-x2,2)+Math.pow(y1-y2,2) <= Math.pow(r1+r2+margin,2);
	}

//	private String getFileExtension(String file_path)
//	{
//		String extension = "";
//
//		int i = file_path.lastIndexOf('.');
//		if (i >= 0) {
//			extension = file_path.substring(i+1);
//		}
//
//		return extension;
//	}

}
