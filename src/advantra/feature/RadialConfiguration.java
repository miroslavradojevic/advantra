package advantra.feature;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 5/10/13
 * Time: 12:27 PM
 * To change this template use File | Settings | File Templates.
 */
public class RadialConfiguration {

    public int          nrRings;

    public double[]     ringRes;      // the width of the ring

    public double[]     ringStart;    // radial distance at which the ring starts

    public RadialConfiguration(
            double[] ring_resol,
            double[] ring_btw_resol
    ){

        nrRings = ring_resol.length;

        ringRes = new double[nrRings];
        for (int i = 0; i < nrRings; i++){
            ringRes[i] = ring_resol[i];
        }

        ringStart = new double[nrRings];

        /*
        assign initial radial distances wrt. given resolution
        and predefined minimal radial distance between
         */

        ringStart[0] = 0;

        for (int i = 1; i < nrRings; i++){
            ringStart[i] = ringStart[i-1] + ringRes[i-1] + ring_btw_resol[i-1];
        }

    }

    public void print(){

         for (int i = 0; i < ringStart.length; i++){
             System.out.println("ring  :  "+ringStart[i]+" < "+ringRes[i]+" > ");
         }

    }


}
