package detection2d;

import java.util.Arrays;

/**
 * Created by miroslav on 3/10/14.
 */
public class Fuzzy2DDemo {

    public static void main(String [] args) {
        System.out.println("TEST FUZZY STUFF");

        int Npoints = 101; ///
        float sigma = 0.25f;

        Fuzzy2D f2d = new Fuzzy2D(Npoints, sigma);
        f2d.showFuzzification();
        f2d.showDefuzzification();

        float[] ii = new float[]{1.0f, 0.5f, 0.1f, 0.1f, 0.9f};
        float[] oo = new float[5];
        f2d.critpointScores(ii[0], ii[1], ii[2], ii[3], ii[4], oo);
        System.out.println(Arrays.toString(oo));

    }

}
