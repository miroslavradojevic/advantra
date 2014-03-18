package demos;

import fit.Fitter1D;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;

import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

/**
 * Created by miroslav on 2/27/14.
 * java -cp "$HOME/critpoint/*:$HOME/jarlib/*" demos.Fitter
 */
public class Fitter {

    public static void main(String[] args) {

        // generate profiles that will be fitted
        ArrayList<float[]> test_profiles = new ArrayList<float[]>();

        /*
            test_profiles
         */
        int profile_length = 14;
        float [] profile_custom;
        profile_custom = new float[]{0.06f, 0.06f, 0.08f, 0.14f, 0.21f, 0.36f, 0.69f, 1.00f, 0.74f, 0.44f, 0.23f, 0.09f, 0.04f, 0.00f}; test_profiles.add(profile_custom);
        profile_custom = new float[]{0.59f, 0.82f, 0.99f, 1.00f, 0.87f, 0.64f, 0.52f, 0.43f, 0.30f, 0.17f, 0.08f, 0.02f, 0.00f, 0.01f}; test_profiles.add(profile_custom);
        profile_custom = new float[]{0.28f, 0.30f, 0.36f, 0.49f, 0.63f, 0.84f, 1.00f, 0.96f, 0.72f, 0.40f, 0.16f, 0.07f, 0.01f, 0.00f}; test_profiles.add(profile_custom);
        profile_custom = new float[]{0.52f, 0.79f, 1.00f, 1.00f, 0.98f, 0.83f, 0.71f, 0.49f, 0.28f, 0.13f, 0.03f, 0.01f, 0.00f, 0.01f}; test_profiles.add(profile_custom);
        profile_custom = new float[]{0.40f, 0.71f, 0.53f, 0.60f, 0.67f, 0.68f, 0.60f, 0.34f, 0.29f, 0.48f, 0.39f, 0.00f, 0.46f, 1.00f}; test_profiles.add(profile_custom);
        // add here more samples

        // random profiles
        int nr_random_profiles = 5;
        Random generator = new Random();
        for (int ii=0; ii<nr_random_profiles; ii++) {

            float[] random_vec = new float[profile_length];

            float min_val = Float.MAX_VALUE;
            float max_val = Float.NEGATIVE_INFINITY;
            for (int aa=0; aa<profile_length; aa++) {
                random_vec[aa] = generator.nextFloat();
                if (random_vec[aa]>max_val) max_val = random_vec[aa];
                if (random_vec[aa]<min_val) min_val = random_vec[aa];
            }
            // normalize
            for (int aa=0; aa<profile_length; aa++) {
                random_vec[aa] = (random_vec[aa]-min_val)/(max_val-min_val);
            }

            test_profiles.add(random_vec);
        }

        int nr_gaussian_profiles = 5;
        for (int ii=0; ii<nr_gaussian_profiles; ii++) {
            float[] gauss_vec = new float[profile_length];
            float std   = (generator.nextFloat()-.5f) * profile_length * 0.2f;
            for (int aa=0; aa<profile_length; aa++) {
                gauss_vec[aa] = (float) Math.exp(-(aa-profile_length/2f)*(aa-profile_length/2f)/(2*std*std));
            }
            test_profiles.add(gauss_vec);
        }

        float[] xaxis = new float[profile_length];
        for (int aa=0; aa<profile_length; aa++) xaxis[aa] = aa;

        // print the profiles
        ImageStack test_profiles_img = new ImageStack(528, 255);
        for (int ii=0; ii<test_profiles.size(); ii++) {
            Plot p = new Plot("", "", "", xaxis, test_profiles.get(ii));
            test_profiles_img.addSlice(p.getProcessor());
        }
        new ImagePlus("input profiles", test_profiles_img).show();

        Fitter1D f1d = new Fitter1D(profile_length, true);
        f1d.showTemplates();

		IJ.log("fitting a family of profiles...");


        float[] fitting_ncc, fitting_mse; // fit output variables
        ImageStack fit_ncc      = new ImageStack(528, 255);
        ImageStack fit_mse      = new ImageStack(528, 255);
        ImageStack fit_scores   = new ImageStack(528, 255);

        // store fit scores
        float[] fit_scores_ncc = new float[test_profiles.size()];
        float[] fit_scores_mse = new float[test_profiles.size()];

        for (int prof_idx = 0; prof_idx<test_profiles.size(); prof_idx++) {

            fitting_ncc = f1d.fit(test_profiles.get(prof_idx), "NCC");
            Plot p_ncc = new Plot("ncc", "idx", "value", xaxis, test_profiles.get(prof_idx)); p_ncc.draw();

            fit_scores_ncc[prof_idx] = fitting_ncc[1];
            int idx_fit = Math.round(fitting_ncc[0]);

            p_ncc.setColor(Color.RED);
            p_ncc.setLineWidth(3);
            p_ncc.addPoints(xaxis, f1d.templates.get(idx_fit), Plot.LINE);
            p_ncc.draw();
            fit_ncc.addSlice("ncc: "+f1d.template_legend.get(idx_fit), p_ncc.getProcessor());

            fitting_mse = f1d.fit(test_profiles.get(prof_idx), "MSE");
            Plot p_mse = new Plot("mse", "idx", "value", xaxis, test_profiles.get(prof_idx)); p_mse.draw();

            fit_scores_mse[prof_idx] = fitting_mse[1];
            idx_fit = Math.round(fitting_mse[0]);

            p_mse.setColor(Color.GREEN);
            p_mse.setLineWidth(3);
            p_mse.addPoints(xaxis, f1d.templates.get(idx_fit), Plot.LINE);
            p_mse.draw();
            fit_mse.addSlice("mse: "+f1d.template_legend.get(idx_fit), p_mse.getProcessor());

        }

        // plot all scores: fit_scores_ncc, fit_scores_mse
        float[] xaxis1 = new float[test_profiles.size()];
        for (int aa=0; aa<xaxis1.length; aa++) xaxis1[aa] = aa;

        float all_max = Float.NEGATIVE_INFINITY;
        float all_min = Float.MAX_VALUE;
        for (int aa=0; aa<fit_scores_ncc.length; aa++) {

            if (fit_scores_ncc[aa]>all_max) all_max=fit_scores_ncc[aa];
            if (fit_scores_mse[aa]>all_max) all_max=fit_scores_mse[aa];

            if (fit_scores_ncc[aa]<all_min) all_min=fit_scores_ncc[aa];
            if (fit_scores_mse[aa]<all_min) all_min=fit_scores_mse[aa];

        }

        Plot p_all_scores = new Plot("", "", "");
        p_all_scores.setLimits(0, xaxis1.length, all_min, all_max);
        p_all_scores.addPoints(xaxis1, fit_scores_ncc, Plot.LINE);
        p_all_scores.addPoints(xaxis1, fit_scores_mse, Plot.LINE);
        p_all_scores.show();

//        new ImagePlus("", p_all_scores.getProcessor()).show();

        new ImagePlus("ncc fit", fit_ncc).show();
        new ImagePlus("mse fit", fit_mse).show();

//        System.exit(0);

    }

}