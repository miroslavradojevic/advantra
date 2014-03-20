package demos;

import detection2d.Fuzzy2D;

import java.util.Arrays;

/**
 * Created by miroslav on 3/10/14.
 * java -cp "$HOME/critpoint/*:$HOME/jarlib/*" demos.Fuzzy
 */
public class Fuzzy {

    public static void main(String [] args) {
        System.out.println("TEST FUZZY STUFF");

        int Npoints = 41;       // nr points to generate out distribution
        float sigma_ON  = 0.2f;    //
        float sigma_OFF = 0.2f;    //
        float mu_ON = 0f;           //
        float mu_OFF = 1f;          //
        System.out.println("mu_ON = "+mu_ON);
        System.out.println("mu_OFF = "+mu_OFF);

        Fuzzy2D f2d = new Fuzzy2D(Npoints, mu_ON, sigma_ON, mu_OFF, sigma_OFF);
        f2d.showFuzzification();
        f2d.showDefuzzification();

        float[] ii, oo;

		/*
        example 1
         */
		ii = new float[]{1.0f, 0.77f, 0.0f, 0.0f};
        oo = new float[5];
        f2d.critpointScores(ii[0], ii[1], ii[2], ii[3], oo);

        for (int kk=0; kk<oo.length; kk++) System.out.print(String.format("%4d\t", kk));
        System.out.println();
        for (int kk=0; kk<oo.length; kk++) System.out.print(String.format("%1.2f\t", oo[kk]));
		System.out.println();

		/*
		example 2
		 */
		System.out.println("example 2");
		ii = new float[]{0.66f, 0.52f};
		oo = new float[5];
		f2d.critpointScores(ii[0], ii[1], oo);

		for (int kk=0; kk<oo.length; kk++) System.out.print(String.format("%4d\t", kk));
		System.out.println();
		for (int kk=0; kk<oo.length; kk++) System.out.print(String.format("%1.2f\t", oo[kk]));
		System.out.println();




	}

}
