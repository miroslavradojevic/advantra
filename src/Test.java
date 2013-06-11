import ij.IJ;
import ij.ImagePlus;
import ij.plugin.PlugIn;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 6/11/13
 * Time: 2:53 PM
 */
public class Test implements PlugIn {

    public void run(String s) {

        int r = 15, t = 5;

        Conf c = new Conf(r, t);

        System.out.println("\n\n SUMMARY: \n\n");
        System.out.println("R = "+r+", T = "+t);
        System.out.println(""+c.regionIdxMap.size()+" index maps");

        for (int m = 0; m<c.regionIdxMap.size(); m++) {

            System.out.println("\nReg. "+m+" :");
            System.out.println(""+c.names.get(m)+" ");
            for (int i = 0; i<c.regionSize.get(m).length; i++) {
                System.out.print("["+i+" -> "+c.regionSize.get(m)[i]+"],");
            }

        }

//        System.out.println(""+c.angles.size()+" angle");
//        System.out.println(""+c.names.size()+" names");
        System.out.println("\n\n --- \n\n");


        ImagePlus showC = new ImagePlus("("+c.r+","+c.diam+")", c.plot());

        showC.show();
        showC.getCanvas().zoomIn(0, 0);
        showC.getCanvas().zoomIn(0, 0);
        showC.getCanvas().zoomIn(0, 0);
        showC.getCanvas().zoomIn(0, 0);

    }

}