package demos;

import fit.Fitter1D;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;

import java.awt.*;
//import java.util.Random;

/**
 * Created by miroslav on 2/27/14.
 * java -cp "$HOME/critpoint/*:$HOME/jarlib/*" demos.Fitter
 */
public class Fitter {

    public static void main(String[] args) {

        int profile_length = 121;
        Fitter1D f1d = new Fitter1D(profile_length, true);
        f1d.showTemplates();

		IJ.log("fitting a vector...");
        //Random generator = new Random();
        float[] random_vec = new float[profile_length];
        for (int aa=0; aa<profile_length; aa++) {
            random_vec[aa] = (float) Math.exp(-(aa-profile_length/2f)*(aa-profile_length/2f)/(2*5*5));//generator.nextFloat();
        }

        // fit
        float[] fitting_ncc, fitting_ssd;
        fitting_ncc = f1d.fit(random_vec, "NCC");
        IJ.log("fit NCC, idx = "+IJ.d2s(fitting_ncc[0], 0)+", NCC = "+IJ.d2s(fitting_ncc[1],2));

        fitting_ssd = f1d.fit(random_vec, "SSD");
        IJ.log("fit SSD, idx = "+IJ.d2s(fitting_ssd[0], 0)+", SSD = "+IJ.d2s(fitting_ssd[1],2));

        // show fit

        ImageStack is_fit = new ImageStack(528, 255);

        float[] xaxis = new float[profile_length];
        for (int aa=0; aa<profile_length; aa++) xaxis[aa] = aa;
        Plot p_ncc = new Plot("ncc", "idx", "value", xaxis, random_vec); p_ncc.draw();
        Plot p_ssd = new Plot("ssd", "idx", "value", xaxis, random_vec); p_ssd.draw();

        int idx_fit = Math.round(fitting_ncc[0]);
        p_ncc.setColor(Color.RED);
        p_ncc.setLineWidth(3);
        p_ncc.addPoints(xaxis, f1d.templates.get(idx_fit), Plot.LINE);
        p_ncc.draw();
        is_fit.addSlice("ncc: "+f1d.template_legend.get(idx_fit), p_ncc.getProcessor());

        idx_fit = Math.round(fitting_ssd[0]);
        p_ssd.setColor(Color.GREEN);
        p_ssd.addPoints(xaxis, f1d.templates.get(idx_fit), Plot.LINE);
        p_ssd.draw();
        is_fit.addSlice("ssd: "+f1d.template_legend.get(idx_fit), p_ssd.getProcessor());

        new ImagePlus("fit", is_fit).show();

    }

}