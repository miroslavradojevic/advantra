import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by miroslav on 28-10-14.
 */
public class SwcDistanceComputer extends Thread  {

    private int begN, endN;

    public static ArrayList<float[]> nodes_A = new ArrayList<float[]>();
    public static ArrayList<float[]> nodes_B = new ArrayList<float[]>();

//    public static float[][] dist;       // matrix with euclidean distances
    public static float[] dAB;          // this one will be filled in threaded run
    public static float[] dBA;          // this one is appended after the run

    public SwcDistanceComputer (int n0, int n1) {
        this.begN = n0;
        this.endN = n1;
    }

    public static void load(ArrayList<float[]> _swc_A, ArrayList<float[]> _swc_B, float dst)
    {

        nodes_A.clear();
        nodes_A = _swc_A;    // nodes_A.size()

        nodes_B.clear();
        nodes_B = _swc_B;    // nodes_B.size()

        if (nodes_A.size()<=0 || nodes_B.size()<=0) {
            System.out.println("Error: empty swc file!");
            return;
        }

            // initialize the matrix with distances and outputs
//        dist = new float[nodes_A.size()][nodes_B.size()];
        dAB = new float[nodes_A.size()];
        Arrays.fill(dAB, Float.POSITIVE_INFINITY);
        dBA = new float[nodes_B.size()]; // will be calculated in a separate method
        Arrays.fill(dBA, Float.POSITIVE_INFINITY);
    }

    public void run()
    {

        // what will be calculated in parallel are the distances and the closest from A->B
        // distances are kept so that B->A can be calculated as well
        for (int locA=begN; locA<endN; locA++) { // threading works per nodes in neuron A

            // calculate differences for point in A towards all the points in B

            float curr_min = Float.POSITIVE_INFINITY; // reset upon every row
            for (int locB = 0; locB < nodes_B.size(); locB++) {

                // since it is 2d it takes into account x and y only
                float dx = nodes_A.get(locA)[ReadSWC.XCOORD]-nodes_B.get(locB)[ReadSWC.XCOORD];
                float dy = nodes_A.get(locA)[ReadSWC.YCOORD]-nodes_B.get(locB)[ReadSWC.YCOORD];
                float dz = nodes_A.get(locA)[ReadSWC.ZCOORD]-nodes_B.get(locB)[ReadSWC.ZCOORD];

                float dd = (float) Math.sqrt( Math.pow(dx,2) + Math.pow(dy,2) + Math.pow(dz,2) );

                if (dd<curr_min) curr_min = dd;

                dBA_store(dd,locB);

            }

            dAB[locA] = curr_min;

        }

//        System.out.println(begN + " -- " + endN +" : "+ Arrays.toString(dBA));

    }

    private static synchronized void dBA_store(float val, int idx) {
        if (val<dBA[idx])
            dBA[idx] = val;
    }

    public static float dAB(){ // directed divergence

        // return average of dAB array with mins. (closest dist of node in A towards B)
        float out = 0;
        for (int i = 0; i < dAB.length; i++) {
            out += dAB[i];
        }
        return out;///dAB.length;

    }

    private static void dAB(float _dst, float[] _count_sum) {

        // average of those that are away from each other at least _dst pixels
        int cnt = 0;
        float sum = 0;
        for (int i = 0; i < dAB.length; i++) {
            if (dAB[i]>=_dst) {
                cnt++;
                sum+=dAB[i];
            }
        }

        _count_sum[0] = cnt;
        _count_sum[1] = sum;

    }

    public static float dBA(){ // directed divergence

        // return average of dBA array with mins.
        float out = 0;
        for (int i = 0; i < dBA.length; i++) {
            out += dBA[i];
        }
        return out;///dBA.length;

    }

    private static void dBA(float _dst, float[] _count_sum) { // returns [count, sum]

        // average of those that are away from each other at least _dst pixels
        int cnt = 0;
        float sum = 0;
        for (int i = 0; i < dBA.length; i++) {
            if (dBA[i]>=_dst) {
                cnt++;
                sum+=dBA[i];
            }
        }

        _count_sum[0] = cnt;
        _count_sum[1] = sum;

    }

    public static float SD(){

//        return (dAB.length+dBA.length>0)? ((dAB()+dBA())/(dAB.length+dBA.length)) : 0;
        return (dAB.length+dBA.length>0)? (float) (.5 * (dAB() / dAB.length) + .5 * (dBA() / dBA.length)) : 0; // alternative way to compute (as in vaa3d)

    }

    public static float SSD(float _dst){ // substanial spatial distance

        // average distance of those that are apart from the other neuron at least _dst pixels
        float[] count_sum_AB = new float[2];
        dAB(_dst, count_sum_AB);

        float[] count_sum_BA = new float[2];
        dBA(_dst, count_sum_BA);

//        return (count_sum_AB[0]+count_sum_BA[0]>0)?( (count_sum_AB[1]+count_sum_BA[1]) / (count_sum_AB[0]+count_sum_BA[0]) ):0;
        return (count_sum_AB[0]+count_sum_BA[0]>0)? (float) (.5 * (count_sum_AB[1] / count_sum_AB[0]) + .5 * count_sum_BA[1] / count_sum_BA[0]) : 0 ; // alternative way to compute (as in vaa3d)

    }

    public static float percSSD(float _dst) {

        // percentage of those that are above _dst
        // serves as a robust indicator on how consistent two reconstructions are

        // A->B
        int cntAB = 0;
        for (int i = 0; i < dAB.length; i++) if (dAB[i] >= _dst) cntAB++;

        // B->A
        int cntBA = 0;
        for (int i = 0; i < dBA.length; i++) if (dBA[i] >= _dst) cntBA++;

//        return (cntAB+cntBA)/(float)(dAB.length+dBA.length);
        return (float) (.5*(cntAB/(float)dAB.length)+.5*(cntBA/(float)dBA.length)); // alternative way to compute (as in vaa3d)

    }

}
