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
//import java.util.Random;

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
            for (int aa=0; aa<profile_length; aa++) {
                random_vec[aa] = generator.nextFloat();
            }
        }

        int nr_gaussian_profiles = 5;
        for (int ii=0; ii<nr_gaussian_profiles; ii++) {
            float[] gauss_vec = new float[profile_length];
            float std   = (generator.nextFloat()-.5f) * profile_length * 0.2f;
            float shift = (generator.nextFloat()-.5f) * profile_length * 0.1f;
            for (int aa=0; aa<profile_length; aa++) {
                gauss_vec[aa] = (float) Math.exp(-(aa-profile_length/2f+shift)*(aa-profile_length/2f+shift)/(2*std*std));
            }
        }



        Fitter1D f1d = new Fitter1D(profile_length, true);
        f1d.showTemplates();

		IJ.log("fitting a family of profiles...");


        // loop here


        // fit
//        float[] fitting_ncc, fitting_ssd;
//        fitting_ncc = f1d.fit(random_vec, "NCC");
//        IJ.log("fit NCC, idx = "+IJ.d2s(fitting_ncc[0], 0)+", NCC = "+IJ.d2s(fitting_ncc[1],2));
//
//        fitting_ssd = f1d.fit(random_vec, "SSD");
//        IJ.log("fit SSD, idx = "+IJ.d2s(fitting_ssd[0], 0)+", SSD = "+IJ.d2s(fitting_ssd[1],2));
//
//        // show fit
//
//        ImageStack ncc_fit = new ImageStack(528, 255);
//        ImageStack mse_fit = new ImageStack(528, 255);
//
//        float[] xaxis = new float[profile_length];
//        for (int aa=0; aa<profile_length; aa++) xaxis[aa] = aa;
//        Plot p_ncc = new Plot("ncc", "idx", "value", xaxis, random_vec); p_ncc.draw();
//        Plot p_ssd = new Plot("ssd", "idx", "value", xaxis, random_vec); p_ssd.draw();
//
//        int idx_fit = Math.round(fitting_ncc[0]);
//        p_ncc.setColor(Color.RED);
//        p_ncc.setLineWidth(3);
//        p_ncc.addPoints(xaxis, f1d.templates.get(idx_fit), Plot.LINE);
//        p_ncc.draw();
//        is_fit.addSlice("ncc: "+f1d.template_legend.get(idx_fit), p_ncc.getProcessor());
//
//        idx_fit = Math.round(fitting_ssd[0]);
//        p_ssd.setColor(Color.GREEN);
//        p_ssd.addPoints(xaxis, f1d.templates.get(idx_fit), Plot.LINE);
//        p_ssd.draw();
//        is_fit.addSlice("ssd: "+f1d.template_legend.get(idx_fit), p_ssd.getProcessor());
//
//        new ImagePlus("fit", is_fit).show();

    }

}