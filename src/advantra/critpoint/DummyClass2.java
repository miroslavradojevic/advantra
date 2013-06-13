package advantra.critpoint;

import ij.ImagePlus;
import ij.plugin.PlugIn;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 6/7/13
 * Time: 5:22 PM
 */
public class DummyClass2 implements PlugIn {

    public void run(String s) {

		ImagePlus inputimg = new ImagePlus("/home/miroslav/d2.tif");
		inputimg.show();

        CircularConfiguration3 ccf3 = new CircularConfiguration3(10);
		//new ImagePlus("", ccf3.plotKernels()).show();
        new ImagePlus("kernel.0", ccf3.plotKernel(0)).show();

		new ImagePlus("", ccf3.scoreAllRot(0, inputimg.getProcessor())).show();

    }

}
