package detection;

import aux.Tools;
import ij.ImagePlus;
import ij.gui.ImageCanvas;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 9/11/13
 * Time: 1:58 PM
 */
public class IntensityAlongLine implements PlugInFilter {

	ImagePlus 		imp;
	ImageCanvas 	cnv;

	public int setup(String s, ImagePlus imagePlus) {

		if(imagePlus==null) return DONE;
		imp = Tools.convertToFloatImage(imagePlus);

		cnv = imagePlus.getCanvas();

		return DOES_8G+DOES_32+NO_CHANGES;

	}

	public void run(ImageProcessor imageProcessor) {

		int Y = imp.getHeight();
		int X = imp.getWidth();

		float ox1 = 20;//0.2f * X;
		float oy1 = 0;//0.2f * Y;

		float ox2 = 10;//0.5f * X;
		float oy2 = 20;//0.5f * Y;

		Overlay my_Ovl = pointsAlongLine(ox1, oy1, ox2, oy2, imp.getProcessor());

		cnv.setOverlay(my_Ovl);

	}

	private static Overlay pointsAlongLine(float ox1, float oy1, float ox2, float oy2, ImageProcessor ipSource)
	{
		float 	dl 	= 0.75f;
		float 	l 	= (float) Math.sqrt(Math.pow(ox2 - ox1, 2) + Math.pow(oy2 - oy1, 2));

		float 	avg = 0;
		int 	cnt = 0;

		int w = ipSource.getWidth();
		int h = ipSource.getHeight();

		Overlay ov = new Overlay();

		if (ox1<0 || oy1<0)
			return ov;


		for (float x = ox1, y = oy1, ll = 0; ll <= l; x+=dl*(ox2-ox1)/l, y+=dl*(oy2-oy1)/l, ll+=dl) { // x <= ox2 && y <= oy2

			if (x>=w-1 || y>=h-1)
				break;

			PointRoi p = new PointRoi(x+.5, y+.5);
			ov.add(p);
		}
		return ov;
	}

}