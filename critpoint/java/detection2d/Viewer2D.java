package detection2d;

import aux.ReadDET;
import aux.Tools;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.*;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;

import java.awt.*;
import java.io.File;

/**
 * Created by miroslav on 26-7-14.
 */
public class Viewer2D implements PlugIn {

	public void run(String s) {

		IJ.open();
		ImagePlus img = IJ.getImage();

		float r = -1;

		String det_path = img.getOriginalFileInfo().directory;
		boolean show_dirs = false;

		GenericDialog gd = new GenericDialog("Select Detection");
		gd.addStringField("detection (name.det) path", det_path, 50);
		gd.addNumericField("regionradius", r, 1);
		gd.addCheckbox("include_directions", true);
		gd.showDialog();

		if (gd.wasCanceled()) return;

		det_path = gd.getNextString();
		r = (float) gd.getNextNumber();
		show_dirs = gd.getNextBoolean();

		/* read detection */
		File f_det = new File(det_path);

		if (!f_det.exists()) {
			System.out.println("file does not exist");
			return;
		}

		if (!Tools.getFileExtension(det_path).equals("det")) {
			System.out.println("file needs to be .det");
			return;
		}

		// read csv file (.det) at specified location
		ReadDET det_reader = new ReadDET(det_path);

		det_reader.print();

		// loop through read values and extract an Overlay for the viewer

		Overlay ov_read = new Overlay();

		for (int i = 0; i < det_reader.x.size(); i++) { // loop read detections

			float dx,dy;

			float x = det_reader.x.get(i);
			float y = det_reader.y.get(i);
			if (r==-1) r = det_reader.r.get(i); // will change the size of the detection overlay circles

			boolean body_added = false;

			Color clr=null;

			if (det_reader.t.get(i).equals("END")) {
//				clr = Color.YELLOW;
				clr = new Color(1, 1, 0, 0.6f);
				body_added = true;
			}
			else if (det_reader.t.get(i).equals("BIF")) {
//				clr = Color.RED;
				clr = new Color(1, 0, 0, 0.6f);
				body_added = true;
			}
			else if (det_reader.t.get(i).equals("CROSS")){
//				clr = Color.GREEN;
				//clr = new Color(0, 1, 0, 0.6f);
				clr = new Color(1, 0, 0, 0.6f);
				body_added = true;
			}

			OvalRoi regroi = new OvalRoi(x-r, y-r, 2*r, 2*r);
			regroi.setStrokeWidth(1);
			regroi.setStrokeColor(clr);
			regroi.setFillColor(clr);
//			if (det_reader.t.get(i).equals("END"))
				ov_read.add(regroi);

			if (show_dirs && body_added) { // add directions if they are !NaN

				float scale_direction = 2f; // just for visualization
				int nr_directions = det_reader.v.get(i).length;
				for (int j = 0; j < nr_directions; j++) {

					dx = det_reader.v.get(i)[j][0];
					dy = det_reader.v.get(i)[j][1];

                    if (!Float.isNaN(dx) && !Float.isNaN(dy)) {
                        Line l = new Line(x+r*dx, y+r*dy, x+scale_direction*r*dx, y+scale_direction*r*dy);
                        l.setStrokeWidth(1f);
                        l.setStrokeColor(clr);
                        l.setFillColor(clr);
                        ov_read.add(l);
                    }

				}


			}

		}

		System.out.println("adding an overlay...");

		img.setOverlay(ov_read);
		img.show();

	}
}
