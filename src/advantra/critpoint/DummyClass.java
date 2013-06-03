package advantra.critpoint;

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 5/31/13
 * Time: 10:12 AM
 * To change this template use File | Settings | File Templates.
 */
public class DummyClass implements PlugIn {


	public void run(String s) {

		ImagePlus img = convertToFloatImage(IJ.openImage());
		img.show();

		int d = 15;//img.getHeight();
		Conf3 c3 = new Conf3(d);
		new ImagePlus("all", c3.plotAll()).show();
		new ImagePlus("one", c3.plotOne(0)).show();

		ImageProcessor ip2 = new FloatProcessor(img.getWidth(), img.getHeight());
		// compare with calculation per point
		for (int x = 0; x < img.getWidth(); x++) {
			for (int y =0; y < img.getHeight(); y++) {

				float val = c3.scoreAtPos(x, y, 0, img.getProcessor());
				ip2.setf(x, y, val);

			}
		}
//
//		ImagePlus      im2 = new ImagePlus("score.pp."+0, ip2);
//		im2.show();


	}

	private ImagePlus convertToFloatImage(
												 ImagePlus inim
	)
	{

		int W = inim.getWidth();
		int H = inim.getHeight();

		ImageProcessor ip = new FloatProcessor(W, H);
		for (int i = 0; i < H*W; i++ ) {
			ip.setf(i, inim.getProcessor().getPixelValue(i%W, i/W));
		}

		ImagePlus outim = new ImagePlus(inim.getTitle()+"ToFloat", ip);
		return outim;

	}

}
