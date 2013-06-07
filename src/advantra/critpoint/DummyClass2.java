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

        CircularConfiguration3 ccf3 = new CircularConfiguration3(10);
        new ImagePlus("", ccf3.plotKernel(2).getProcessor(1)).show();

    }

}
