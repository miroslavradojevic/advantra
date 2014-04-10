package demos;

import aux.Stat;
import detection2d.FuzzyBranch2D;

import java.util.Arrays;

/**
 * Created by miroslav on 3/25/14.
 * java -cp "$HOME/critpoint/*:$HOME/jarlib/*" demos.Dummy
 */
public class Dummy {
	public static void main(String[] args) {
		System.out.println("Dummy working...");
		float[] ar = new float[]{1, 2, 3, 4, 5};
		System.out.println(Arrays.toString(ar));
		Stat.normalize(ar);
		System.out.println(Arrays.toString(ar));

        FuzzyBranch2D test = new FuzzyBranch2D(20, 1, 0.1f, 0.6f, 0.1f);
        test.showFuzzification();
        test.showDefuzzification();
        float[] out = new float[test.agg.length];
        test.critpointScores(1f, 1f, out);
        test.showAgg();

	}
}
