package fit;

import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;
import ij.process.ImageProcessor;

import java.util.ArrayList;

/**
 * Created by miroslav on 2/7/14.
 * Fitting gaussian 1d profile given the input 1d profile vector
 * class contains the library of 1d profiles that can be expected with neuron from projecting along 2d patches
 * (if calculated averages along verticals of the patch)
 * all the profiles are the expected gaussian bells typical for the neuron, normalized 0-1
 * it's the shape that's being matched and quantification of that fit through ssd, ncc, etc...
 *
 */
public class Fitter1D {

    int vector_len;

    // template profiles
    ArrayList<float[]> templates;

    // middle
    float middle_idx, start_sigma, end_sigma, start_width, end_width, d_sigma = 1f, d_width = 2f;

    public Fitter1D(int _vector_len){

        vector_len  = _vector_len;
        start_sigma = vector_len*0.25f;
        end_sigma   = vector_len*0.95f;
        d_sigma     = vector_len*0.05f;
        start_width = 0;
        end_width   = vector_len*0.5f;
        d_width     = vector_len*0.1f;

        middle_idx = (vector_len-1) / 2f;

        templates = new ArrayList<float[]>();

        // width loop, sigma loop
        for(float wd = start_width; wd<=end_width; wd+=d_width) {

            for (float sigma = start_sigma; sigma < end_sigma; sigma+=d_sigma) {

                float[] templates_element = new float[vector_len];

//                float min_val = Float.MAX_VALUE, max_val = Float.NEGATIVE_INFINITY;

                for (int i = 0; i < vector_len; i++) {

                    if (i<middle_idx-wd/2) {
                        double d = middle_idx - wd/2 - i;
                        templates_element[i] = (float) Math.exp(-Math.pow(d, 2)/(2*sigma));
                    }
                    else if (i>=middle_idx-wd/2 && i<middle_idx+wd/2) {
                        float d = i - middle_idx;
                        templates_element[i] = 1;
                        //d = (d<0)? -d : d ;
                        //templates_element[i] = (float) Math.exp(-Math.pow(d, 2)/(2*sigma));
                    }
                    else {
                        float d = i - (middle_idx + wd/2);
                        templates_element[i] = (float) Math.exp(-Math.pow(d, 2)/(2*sigma));
                    }

//                    if (templates_element[i] > max_val) {
//                        max_val = templates_element[i];
//                    }
//                    if (templates_element[i] < min_val) {
//                        min_val = templates_element[i];
//                    }

                }

                templates.add(templates_element);

            }

        }

    }

    public float[] ssd(float[] profile) { // returns array[2] []

        float[] out = new float[]{Float.NaN, Float.MAX_VALUE}; // best score (index, ssd score)

        for (int i=0; i<templates.size(); i++) {

            float   curr_ssd_score  = 0;

			// profile has to have the same length
			for (int j =0; j<templates.get(i)[j]; j++) curr_ssd_score += Math.pow(profile[j] - templates.get(i)[j], 2);

			if (curr_ssd_score < out[1]) {
				out[1] = curr_ssd_score;
				out[0] = i;
			}

        }

        return out;

    }

	public float[] ncc() {

		float[] out = new float[]{Float.NaN, Float.MAX_VALUE};

		return out;

	}

    public void showTemplates() {

        ImageStack is_out = new ImageStack(528, 255);

        float[] xaxis = new float[vector_len];
        for (int aa=0; aa<vector_len; aa++) xaxis[aa] = aa;

        for (int aa= 0; aa<templates.size(); aa++) { // plot each template

            Plot p = new Plot("", "", "", xaxis, templates.get(aa));
            is_out.addSlice(p.getProcessor());

        }

        ImagePlus img_out = new ImagePlus("templates", is_out);
        img_out.show();

    }

	public float[] getTemplate(int template_index) {

		return templates.get(template_index);

	}

	public Plot getTemplatePlot(int template_index) {

		float[] xaxis = new float[vector_len];
		for (int aa=0; aa<vector_len; aa++) xaxis[aa] = aa;

		Plot p = new Plot("", "", "", xaxis, templates.get(template_index));
		return p;

	}

}
