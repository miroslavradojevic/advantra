package detection2d;

import aux.ReadSWC;
import ij.ImagePlus;
import ij.Macro;
import ij.Prefs;
import ij.gui.GenericDialog;
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

	int tp_END, fp_END, fn_END;
	float p_END, r_END;

	int tp_CRS, fp_CRS, fn_CRS;
	float p_CRS, r_CRS;

	// BIF and CRS together = JUN
	int tp_JUN, fp_JUN, fn_JUN;
	float p_JUN, r_JUN;

	public void run(String s) { // regular PlugIn run() method

		String			det_file_path;     		// either .tif or .det
		String			det_file_name;          //

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
			gd.addStringField("gndtth_path (.swc)", 	new File(output_dir_name).getParent(), 150);
			gd.addStringField("gndtth_tag", 	"ground_truth_tag", 50);
			gd.showDialog();
			if (gd.wasCanceled()) return;
			gndtth_path	= gd.getNextString();
			gndtth_tag	= gd.getNextString();

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

		if (!getFileExtension(gndtth_path).equals("swc")) {
			System.out.println("file needs to be .swc");
			return;
		}

		ArrayList<float[]> gndtth_juns = new ArrayList<float[]>();
		ArrayList<float[]> gndtth_bifs = new ArrayList<float[]>();
		ArrayList<float[]> gndtth_ends = new ArrayList<float[]>();
		ArrayList<float[]> gndtth_crss = new ArrayList<float[]>();

		ReadSWC reader = new ReadSWC(gndtth_path);

		for (int b=0; b<reader.nodes.size(); b++) {

			int typ = Math.round(reader.nodes.get(b)[reader.TYPE]);
			float bx = reader.nodes.get(b)[reader.XCOORD];
			float by = reader.nodes.get(b)[reader.YCOORD];
			float br = reader.nodes.get(b)[reader.RADIUS];

			if (typ==1) {
				gndtth_ends.add(new float[]{bx, by, br});
			}
			else if (typ==3) {
				gndtth_bifs.add(new float[]{bx, by, br});
				gndtth_juns.add(new float[]{bx, by, br});
			}
			else if (typ==4) {
				gndtth_crss.add(new float[]{bx, by, br});
				gndtth_juns.add(new float[]{bx, by, br});
			}


		}

		System.out.println(
								  "gndtth loaded -> " + reader.nodes.size() + " elements: " +
										  gndtth_juns.size() + " JUNS, " +
										  gndtth_bifs.size() + " BIFS, " +
										  gndtth_crss.size() + " CRSS, " +
										  gndtth_ends.size() + " ENDS."
		);


		/* read the detections */

		ArrayList<float[]> det_bifs = new ArrayList<float[]>();
		ArrayList<float[]> det_juns = new ArrayList<float[]>();
		ArrayList<float[]> det_ends = new ArrayList<float[]>();
		ArrayList<float[]> det_crss = new ArrayList<float[]>();

		det_file_name = new File(det_file_path).getName();
		String detection_file_ext = getFileExtension(det_file_path);

		// append to the list in different scenarios
		if (detection_file_ext.equals("tif")) {

			ImagePlus ip_load = new ImagePlus(det_file_path);
			if(ip_load==null) {
				System.out.println("loaded image was NULL");
				return;
			}

			Overlay ov_det = ip_load.getOverlay();

			for (int i = 0; i <ov_det.size(); i++) {

				if (ov_det.get(i).getTypeAsString().equals("Oval"))  {

					float x = (float) ov_det.get(i).getFloatBounds().getX();
					float y = (float) ov_det.get(i).getFloatBounds().getY();
					float w = (float) ov_det.get(i).getFloatBounds().getWidth();
					float h = (float) ov_det.get(i).getFloatBounds().getHeight();
					float r = (float) (Math.max(w, h)/2);

					float xc = x+r/1-.0f;
					float yc = y+r/1-.0f;

					Color obj_color = ov_det.get(i).getFillColor();

					if (obj_color==null) {
						System.out.println("found null color in overlay");
						continue;
					}

					if (obj_color.equals(Color.RED)) {
						det_bifs.add(new float[]{xc, yc, r});
						det_juns.add(new float[]{xc, yc, r});
					}
					else if (obj_color.equals(Color.YELLOW)) {
						det_ends.add(new float[]{xc, yc, r});
					}
					else if (obj_color.equals(Color.GREEN)) {
						det_crss.add(new float[]{xc, yc, r});
						det_juns.add(new float[]{xc, yc, r});
					}

				}

			}

			System.out.println("tif loaded -> " + ov_det.size() + " elements: " +
									   det_juns.size() + " JUNS, " +
									   det_bifs.size() + " BIFS, " +
									   det_crss.size() + " CRSS, " +
									   det_ends.size() + " ENDS."
			);

		}
		else if (detection_file_ext.equals("det")) {





		}
		else {
			System.out.println("Detection file type not recognized.");
			return;
		}


		/* evaluation */

		doBifEvaluation(det_bifs, gndtth_bifs);
		doJunEvaluation(det_juns, gndtth_juns);
		doEndEvaluation(det_ends, gndtth_ends);
		doCrsEvaluation(det_crss, gndtth_crss);

		/* store output */

		PrintWriter 	logWriter = null;
		String legend = "NAME, ANNOT_TAG, " +
								"TP_JUN, FP_JUN, FN_JUN, P_JUN, R_JUN, " +
								"TP_BIF, FP_BIF, FN_BIF, P_BIF, R_BIF, " +
								"TP_CRS, FP_CRS, FN_CRS, P_CRS, R_CRS, " +
								"TP_END, FP_END, FN_END, P_END, R_END";

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

		String eval = 		    "\""+det_file_name + "\", "
							  + "\""+gndtth_tag    + "\", "
							  + tp_JUN + ", " + fp_JUN + ", " + fn_JUN + ", " + p_JUN + ", " + r_JUN + ", "
							  + tp_BIF + ", " + fp_BIF + ", " + fn_BIF + ", " + p_BIF + ", " + r_BIF + ", "
							  + tp_CRS + ", " + fp_CRS + ", " + fn_CRS + ", " + p_CRS + ", " + r_CRS + ", "
							  + tp_END + ", " + fp_END + ", " + fn_END + ", " + p_END + ", " + r_END;


		System.out.println(legend);
		System.out.println(eval);

		logWriter.println(eval);
		logWriter.close();

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

				if (circlesOverlap(ax, ay, ar, bx, by, br)) {// two roi circles overlap
					found = true;
					tp_JUN++;
					annots[b] = true;
				}
			}

			if (!found) fp_JUN++;  // detected but was not in the list of annotated ones

		}

		for (int a=0; a<annots.length; a++) if (!annots[a]) fn_JUN++;

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
					annots[b] = true;
				}
			}

			if (!found) fp_END++;  // detected but was not in the list of annotated ones

		}

		for (int a=0; a<annots.length; a++) if (!annots[a]) fn_END++;

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

	private String getFileExtension(String file_path)
	{
		String extension = "";

		int i = file_path.lastIndexOf('.');
		if (i >= 0) {
			extension = file_path.substring(i+1);
		}

		return extension;
	}


}
