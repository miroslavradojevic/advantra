import ij.IJ;
import ij.plugin.PlugIn;

import java.util.Arrays;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 7/4/13
 * Time: 12:42 PM
 */
public class TestMapTracker  implements PlugIn {

    public void run(String s) {
        IJ.log("testing...");

        int[][] map = new int[3][];
        map[0] = new int[]{-1, 0, 1, 0, -1, 2};
        map[1] = new int[]{-1, 1, 0, 2, -1, 0};
        map[2] = new int[]{-1, 2, 2, 1, -1, 1};

        IJ.log("map[0]:" + Arrays.toString(map[0]));
        IJ.log("map[1]:" + Arrays.toString(map[1]));
        IJ.log("map[2]:" + Arrays.toString(map[2]));



        for (int q=0; q<map[0].length; q++) {
            for (int i=0; i<3; i++) {

                IJ.log("track: "+i+" from "+q+", got: "+BifDetect.mapTracker(map, i, q));

            }

            IJ.log("---");
        }

    }

}
