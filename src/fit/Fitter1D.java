package fit;

import aux.Stat;
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
    ArrayList<float[]> 	templates;
	ArrayList<Float> 	templates_mean;
	ArrayList<float[]>	templates_mean_subtr;
	ArrayList<Float> 	templates_mean_subtr_sum;
	ArrayList<String>	template_legend;

    // middle
    float middle_idx,
			start_sigma, d_sigma, end_sigma, // slope
			start_width, d_width, end_width, // basic width
			start_shift, d_shift, end_shift; // sift from middle

    public Fitter1D(int _vector_len, boolean verbose){

		// parameters
        vector_len  = _vector_len;

		middle_idx = (vector_len-1) / 2f;
		// sigma
		start_sigma = vector_len*0.25f;
		d_sigma     = vector_len*0.05f;
		end_sigma   = vector_len*0.95f;
		// width
		start_width = 0;
		d_width     = vector_len*0.1f;
		end_width   = vector_len*0.0f;
		//shift
		start_shift	= -vector_len*0.2f;
		d_shift		= vector_len*0.05f;
		end_shift	= -start_shift;

		// variables will be stored in:
		// templates                                  	(vector)
		// templates_mean (used for NCC calculation)  	(scalar)
		// templates - templates_mean                 	(vector)
		// sum(templates - templates_mean)				(scalar)
		// templates_legend 							(string)
        templates 			= new ArrayList<float[]>();
		templates_mean 		= new ArrayList<Float>();
		templates_mean_subtr= new ArrayList<float[]>();
		templates_mean_subtr_sum = new ArrayList<Float>();
		template_legend = new ArrayList<String>();

        // width loop, sigma loop
        for(float width = start_width; width<=end_width; width+=d_width) {

            for (float sigma = start_sigma; sigma < end_sigma; sigma+=d_sigma) {

				for (float shift = start_shift; shift < end_shift; shift+=d_shift) {

					float[] templates_element = new float[vector_len];

					float boundary_1 = middle_idx + shift - width/2;
					float boundary_2 = middle_idx + shift + width/2;

					for (int i = 0; i < vector_len; i++) {

						if (i < boundary_1) {
							double d = boundary_1 - i;
							templates_element[i] = (float) Math.exp(-Math.pow(d, 2)/(2*sigma*sigma));
						}
						else if (i >= boundary_1 && i < boundary_2) {
							templates_element[i] = 1;
						}
						else {
							float d = i - boundary_2;
							templates_element[i] = (float) Math.exp(-Math.pow(d, 2)/(2*sigma*sigma));
						}

					}

					templates.add(templates_element);
					// calculate mean
					templates_mean.add(Stat.average(templates_element));
					// subtract and store in the same array


				}



            }

        }

		if (verbose) {
			System.out.println("created " + templates.size() + " template profiles.");
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

	public float[] ncc(float[] profile) {

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
