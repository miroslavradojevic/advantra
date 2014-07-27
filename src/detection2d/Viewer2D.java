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

		String det_path = img.getOriginalFileInfo().directory;
		boolean show_dirs = false;

		GenericDialog gd = new GenericDialog("Select Detection");
		gd.addStringField("detection.zip path", det_path, 50);
		gd.addCheckbox("include_directions", true);
		gd.showDialog();

		if (gd.wasCanceled()) return;

		det_path = gd.getNextString();
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
			float r = det_reader.r.get(i);

			boolean body_added = false;
			Color clr=null;

			if (det_reader.t.get(i).equals("END")) {
				clr = Color.YELLOW;
				body_added = true;
			}
			else if (det_reader.t.get(i).equals("BIF")) {
				clr = Color.RED;
				body_added = true;
			}
			else if (det_reader.t.get(i).equals("CROSS")){
				clr = Color.GREEN;
				body_added = true;
			}

			OvalRoi regroi = new OvalRoi(x-r, y-r, 2*r, 2*r);
			regroi.setStrokeWidth(1);
			regroi.setStrokeColor(clr);
			regroi.setFillColor(clr);
			ov_read.add(regroi);

			if (show_dirs && body_added) { // add directions

				float scale_direction = 1.3f; // just for visualization
				int nr_directions = det_reader.v.get(i).length;
				for (int j = 0; j < nr_directions; j++) {

					dx = det_reader.v.get(i)[j][0];
					dy = det_reader.v.get(i)[j][1];

					Line l = new Line(x, y, x+scale_direction*r*dx, y+scale_direction*r*dy);
					l.setStrokeWidth(r/2f);
					l.setStrokeColor(clr);
					l.setFillColor(clr);
					ov_read.add(l);

				}


			}

		}

		System.out.println("adding an overlay...");

		img.setOverlay(ov_read);
		img.show();

	}
}
