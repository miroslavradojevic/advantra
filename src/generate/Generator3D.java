package generate;

import aux.AnalyzeSWC;
import ij.ImagePlus;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 10/23/13
 * Time: 4:46 PM
 */
public class Generator3D {
    // will be used to generate 3D image stacks

    private static int margin = 10;

    public static ImagePlus runFromSWC(String swcPath, float k, float snr) {

        // read swc
        AnalyzeSWC aSWC = new AnalyzeSWC(swcPath);

        return new ImagePlus();

    }

}
