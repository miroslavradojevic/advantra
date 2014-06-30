package detection2d;

import aux.ReadSWC;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.Overlay;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;

import java.awt.*;
import java.io.*;
import java.util.ArrayList;

/**
 * Created by miroslav on 30-6-14.
 */
public class Evaluator2D implements PlugIn {

	String			image_path;
	String			image_name;
	String 			gndtth_path;		// ground truth swc file with critical points to evaluate, same name as input image
	String          eval_string = "";
	String 		    output_log_name;
	String		    output_dir_name;
	PrintWriter 	logWriter = null;

	int tp_BIF, fp_BIF, fn_BIF;
	float p_BIF, r_BIF;

	int tp_END, fp_END, fn_END;
	float p_END, r_END;

	int tp_JUN, fp_JUN, fn_JUN;
	float p_JUN, r_JUN;

	int tp_CRS, fp_CRS, fn_CRS;
	float p_CRS, r_CRS;

	ReadSWC reader;
	Overlay ov_det;

	public void run(String s) {

		// load the image with detections through the menu
		String in_folder = Prefs.get("id.folder", System.getProperty("user.home"));
		OpenDialog.setDefaultDirectory(in_folder);
		OpenDialog dc = new OpenDialog("Select detection file");
		in_folder = dc.getDirectory();
		image_path = dc.getPath();
		if (image_path==null) return;
		Prefs.set("id.folder", in_folder);
		output_dir_name = in_folder;//ip_load.getOriginalFileInfo().directory;

		// load swc with annotations
		dc = new OpenDialog("Select annotation swc");
		gndtth_path = dc.getPath();
		if (gndtth_path==null) return;



		output_log_name = output_dir_name + "det.csv"; // + File.separator
		System.out.println(output_log_name);

		File f = new File(output_log_name);
		if (!f.exists()) {
			try {
				logWriter = new PrintWriter(output_log_name);
				logWriter.print("name," +
										"TP_JUN, FP_JUN,\tFN_JUN,\tP_JUN,\tR_JUN,\t" +
										"TP_END, FP_END,\tFN_END,\tP_END,\tR_END,\t" +
										"TP_BIF, FP_BIF,\tFN_BIF,\tP_BIF,\tR_BIF,\t" +
										"TP_CRS, FP_CRS,\tFN_CRS,\tP_CRS,\tR_CRS\n");
				logWriter.close();
			} catch (FileNotFoundException ex) {}
		}
		// if it exists already in the folder, just prepare to append on the existing file
		try {
			logWriter = new PrintWriter(new BufferedWriter(new FileWriter(output_log_name, true)));
		} catch (IOException e) {}

		ImagePlus ip_load = new ImagePlus(image_path);
		if(ip_load==null) return;
		image_name = ip_load.getShortTitle();

		// read overlay, OvalRois are marking the detected regions
		ov_det = ip_load.getOverlay();

//		for (int i = 0; i < ov_det.size(); i++) {
//			System.out.println("" + ov_det.get(i).getFillColor() + "" + ov_det.get(i).getFillColor().equals(Color.RED));
//		}

		reader = new ReadSWC(gndtth_path);

		doBifEvaluation();
		doJunEvaluation();
		doEndEvaluation();

		String eval = "\""+image_name+"\", " 	+ tp_JUN + ", " + fp_JUN + ", " + fn_JUN + ", " + p_JUN + ", " + r_JUN + ", "
							  					+ tp_END + ", " + fp_END + ", " + fn_END + ", " + p_END + ", " + r_END + ", "
							  					+ tp_BIF + ", " + fp_BIF + ", " + fn_BIF + ", " + p_BIF + ", " + r_BIF + ", ";
		System.out.println(eval);
		logWriter.println(eval);
		logWriter.close();

	}

	public void doBifEvaluation()
	{

		ArrayList<float[]> gndtth_bifs = new ArrayList<float[]>();
		for (int b=0; b<reader.nodes.size(); b++) {

			int typ = Math.round(reader.nodes.get(b)[reader.TYPE]);

			if (typ==3) {
				float bx = reader.nodes.get(b)[reader.XCOORD];
				float by = reader.nodes.get(b)[reader.YCOORD];
				float br = reader.nodes.get(b)[reader.RADIUS];
				gndtth_bifs.add(new float[]{bx, by, br});
			}


		}

		ArrayList<float[]> det_bifs = new ArrayList<float[]>();

		for (int i = 0; i <ov_det.size(); i++) {

			if (ov_det.get(i).getTypeAsString().equals("Oval") && ov_det.get(i).getFillColor().equals(Color.RED))  {

				float x = (float) ov_det.get(i).getFloatBounds().getX();
				float y = (float) ov_det.get(i).getFloatBounds().getY();
				float w = (float) ov_det.get(i).getFloatBounds().getWidth();
				float h = (float) ov_det.get(i).getFloatBounds().getHeight();
				float r = (float) (Math.max(w, h)/2);

				float xc = x+r/1-.0f;
				float yc = y+r/1-.0f;

				det_bifs.add(new float[]{xc, yc, r});
			}

		}


		tp_BIF = 0; fp_BIF = 0; fn_BIF = 0;
		boolean[] annots = new boolean[gndtth_bifs.size()];  // necessary for fn calculation
		// loop all detected BIF regions (all sorts of critical points are in the same list now)
		for (int a=0; a<det_bifs.size(); a++) {

				boolean found = false;

				float ax = det_bifs.get(a)[0];
				float ay = det_bifs.get(a)[1];
				float ar = det_bifs.get(a)[2];

				for (int b=0; b<gndtth_bifs.size(); b++) {

					float bx = gndtth_bifs.get(b)[0];
					float by = gndtth_bifs.get(b)[1];
					float br = gndtth_bifs.get(b)[2];

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

	public void doJunEvaluation()
	{

		ArrayList<float[]> gndtth_juns = new ArrayList<float[]>();
		for (int b=0; b<reader.nodes.size(); b++) {

			int typ = Math.round(reader.nodes.get(b)[reader.TYPE]);

			if (typ==3 || typ==4) {
				float bx = reader.nodes.get(b)[reader.XCOORD];
				float by = reader.nodes.get(b)[reader.YCOORD];
				float br = reader.nodes.get(b)[reader.RADIUS];
				gndtth_juns.add(new float[]{bx, by, br});
			}


		}

		ArrayList<float[]> det_juns = new ArrayList<float[]>();

		for (int i = 0; i <ov_det.size(); i++) {

			if (ov_det.get(i).getTypeAsString().equals("Oval") && (ov_det.get(i).getFillColor().equals(Color.RED) || ov_det.get(i).getFillColor().equals(Color.GREEN)) )  {

				float x = (float) ov_det.get(i).getFloatBounds().getX();
				float y = (float) ov_det.get(i).getFloatBounds().getY();
				float w = (float) ov_det.get(i).getFloatBounds().getWidth();
				float h = (float) ov_det.get(i).getFloatBounds().getHeight();
				float r = (float) (Math.max(w, h)/2);

				float xc = x+r/1-.0f;
				float yc = y+r/1-.0f;

				det_juns.add(new float[]{xc, yc, r});
			}

		}


		tp_JUN = 0; fp_JUN = 0; fn_JUN = 0;
		boolean[] annots = new boolean[gndtth_juns.size()];  // necessary for fn calculation
		// loop all detected BIF regions (all sorts of critical points are in the same list now)
		for (int a=0; a<det_juns.size(); a++) {

			boolean found = false;

			float ax = det_juns.get(a)[0];
			float ay = det_juns.get(a)[1];
			float ar = det_juns.get(a)[2];

			for (int b=0; b<gndtth_juns.size(); b++) {

				float bx = gndtth_juns.get(b)[0];
				float by = gndtth_juns.get(b)[1];
				float br = gndtth_juns.get(b)[2];

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

	public void doEndEvaluation()
	{

		ArrayList<float[]> gndtth_ends = new ArrayList<float[]>();
		for (int b=0; b<reader.nodes.size(); b++) {

			int typ = Math.round(reader.nodes.get(b)[reader.TYPE]);

			if (typ==1) {
				float bx = reader.nodes.get(b)[reader.XCOORD];
				float by = reader.nodes.get(b)[reader.YCOORD];
				float br = reader.nodes.get(b)[reader.RADIUS];
				gndtth_ends.add(new float[]{bx, by, br});
			}

		}

		ArrayList<float[]> det_ends = new ArrayList<float[]>();

		for (int i = 0; i <ov_det.size(); i++) {

			if (ov_det.get(i).getTypeAsString().equals("Oval") && ov_det.get(i).getFillColor().equals(Color.YELLOW))  {

				float x = (float) ov_det.get(i).getFloatBounds().getX();
				float y = (float) ov_det.get(i).getFloatBounds().getY();
				float w = (float) ov_det.get(i).getFloatBounds().getWidth();
				float h = (float) ov_det.get(i).getFloatBounds().getHeight();
				float r = (float) (Math.max(w, h)/2);

				float xc = x+r/1-.0f;
				float yc = y+r/1-.0f;

				det_ends.add(new float[]{xc, yc, r});
			}

		}
		System.out.println(det_ends.size() + "ends detected");


		tp_END = 0; fp_END = 0; fn_END = 0;
		boolean[] annots = new boolean[gndtth_ends.size()];  // necessary for fn calculation
		// loop all detected BIF regions (all sorts of critical points are in the same list now)
		for (int a=0; a<det_ends.size(); a++) {

			boolean found = false;

			float ax = det_ends.get(a)[0];
			float ay = det_ends.get(a)[1];
			float ar = det_ends.get(a)[2];

			for (int b=0; b<gndtth_ends.size(); b++) {

				float bx = gndtth_ends.get(b)[0];
				float by = gndtth_ends.get(b)[1];
				float br = gndtth_ends.get(b)[2];

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

	private boolean circlesOverlap(float x1, float y1, float r1, float x2, float y2, float r2)
	{
		return Math.pow(x1-x2,2)+Math.pow(y1-y2,2) <= Math.pow(r1+r2,2);
	}

}
