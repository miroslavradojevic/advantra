package fit;

import aux.Stat;
import ij.IJ;
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
    public ArrayList<float[]> 	templates;
	public ArrayList<Float> 	templates_mean;
	public ArrayList<float[]>	templates_mean_subtr;
	public ArrayList<Float> 	templates_mean_subtr_sumsqr;
	public ArrayList<String>	template_legend;

    // middle
    float middle_idx,
			start_sigma, d_sigma, end_sigma, // slope
			start_width, d_width, end_width, // basic width
			start_shift, d_shift, end_shift; // sift from middle

    static float low_boundary = 0.01f;          // bottom limit at the boundary

    public Fitter1D(int _vector_len, boolean verbose){

		// parameters
        vector_len  = _vector_len;

		middle_idx = (vector_len-1) / 2f;
		// sigma
		start_sigma = vector_len*0.02f;
		d_sigma     = vector_len*0.02f;
		end_sigma   = vector_len*0.16f;
		// width
		start_width = 0;
		d_width     = vector_len*0.05f;
		end_width   = vector_len*0.35f;
		//shift
		start_shift	= 0;//-vector_len*0.15f;
		d_shift		= vector_len*0.05f;
		end_shift	= 0;//-start_shift;

		// variables will be stored in:
		// templates                                  	(vector)
		// templates_mean (used for NCC calculation)  	(scalar)
		// templates - templates_mean                 	(vector)
		// sum(templates - templates_mean)				(scalar)
		// templates_legend 							(string)
        templates 			= new ArrayList<float[]>();
		templates_mean 		= new ArrayList<Float>();
		templates_mean_subtr= new ArrayList<float[]>();
		templates_mean_subtr_sumsqr = new ArrayList<Float>();
		template_legend = new ArrayList<String>();

        // width loop, sigma loop
        for(float width = start_width; width<=end_width; width+=d_width) {

            for (float sigma = start_sigma; sigma <= end_sigma; sigma+=d_sigma) {

				for (float shift = start_shift; shift <= end_shift; shift+=d_shift) {

					float[] templates_element = new float[vector_len];

					float boundary_1 = Math.round(middle_idx + shift - width/2);
					float boundary_2 = Math.round(middle_idx + shift + width/2);

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

                    // check boundary elements are low enough
                    if (templates_element[0] < low_boundary && templates_element[vector_len-1] < low_boundary) {

                        templates.add(templates_element.clone());

                        // calculate mean
                        float mn = Stat.average(templates_element);
                        templates_mean.add(mn);

                        // subtract and store in the same array, add the sum as well
                        float sum = 0;
                        for (int aa = 0; aa < vector_len; aa++) {
                            templates_element[aa] = templates_element[aa] - mn;
                            sum += templates_element[aa] * templates_element[aa];
                        }

                        templates_mean_subtr.add(templates_element.clone());
                        templates_mean_subtr_sumsqr.add(sum);
                        template_legend.add("wdt="+ IJ.d2s(width, 2)+","+"sig="+IJ.d2s(sigma, 2)+","+"shf="+IJ.d2s(shift,2));

                    }

				}

            }

        }

		if (verbose) {
			System.out.println("created " + templates.size() + " template profiles.");
		}

    }

    public float[] fit(float[] profile, String mode) {

        // returns the fitting result: (index of the profile, fitting score)
        float[] out = new float[2];
        out[0] = Float.NaN;
        if (mode.equalsIgnoreCase("MSE"))
        {
            out[1] = Float.POSITIVE_INFINITY; // looking for smallest
        }
        else if (mode.equalsIgnoreCase("NCC")) {
            out[1] = Float.POSITIVE_INFINITY; // looking for highest
        }
        else {
            out[1] = Float.NaN;
            return out;
        }

        // loop the templates
        for (int i=0; i<templates.size(); i++) {

            // calculate score
            float curr_score;

            if (mode.equalsIgnoreCase("MSE")) {
                curr_score = mse(profile, templates.get(i));
                if (curr_score < out[1]) {out[0]=i; out[1]=curr_score;}
            }
            else if (mode.equalsIgnoreCase("NCC")) {
                curr_score = ncc(profile, templates_mean_subtr.get(i), templates_mean_subtr_sumsqr.get(i));
                if (curr_score < out[1]) {out[0]=i; out[1]=curr_score;}
            }

        }

        return out;

    }

    /*
        mean squared error
     */
    private float mse(float[] f, float[] t) { // both f[] and t[] are normalized 0-1
        float sc = 0;
        int cnt = 0;
        for (int aa=0; aa<vector_len; aa++) {
            sc += Math.pow(f[aa] - t[aa], 2);  // score computation
            cnt++;
        }
        return sc/cnt;
    }

    /*
        normalized cross correlation
     */
    private float ncc(float[] f, float[] t_tM, float sumsqr_t_tM){

        float f_mean = 0;
        for (int i=0; i < f.length; i++) f_mean += f[i];
        f_mean = f_mean / (float)f.length;

        float sc = 0;
        float f_sub_f_mean_sumsqr = 0;

        for (int aa=0; aa<vector_len; aa++) {
            sc += (f[aa]-f_mean) * t_tM[aa]; // important that input f is vector_len length
            f_sub_f_mean_sumsqr += (f[aa]-f_mean) * (f[aa]-f_mean);
        }

        float ncc_val;

        if (f_sub_f_mean_sumsqr>Float.MIN_VALUE) {
            ncc_val = sc / (float) Math.sqrt(f_sub_f_mean_sumsqr * sumsqr_t_tM);
            return 1- ncc_val;
        }
        else {
            ncc_val = 0;
            return 1 - ncc_val;
        }

    }

    /*
        normalized sum of absolute differences
        normalized wrt cumulative sum of the template signal
        0 - absolute fit
        1 - discrepancy is comparable to the amount of signal
     */
    /*
    private float nsad(float[] f, float[] t) {
        float sc = 0;
        float norm = 0;
        for (int aa=0; aa<vector_len; aa++) {sc+=Math.abs(f[aa]-t[aa]); norm+=Math.abs(t[aa]);}
        return sc / norm;
    }
    */

    /*
        normalized sum of squared differences
        normalized wrt cumulative sum of the template signal squared
        0 - absolute fit
        1 - discrepancy comparable to the amount of signal (higher differences weighted more)
     */
    /*
    private float nssd(float[] f, float[] t) {
        float sc = 0;
        float norm = 0;
        for (int aa=0; aa<vector_len; aa++) {sc+=(f[aa]-t[aa])*(f[aa]-t[aa]); norm+=t[aa]*t[aa];}
        return sc / norm;
    }
    */

    public void showTemplates() {

        ImageStack is_out = new ImageStack(528, 255);

        float[] xaxis = new float[vector_len];
        for (int aa=0; aa<vector_len; aa++) xaxis[aa] = aa;

        for (int aa= 0; aa<templates.size(); aa++) { // plot each template
            Plot p = new Plot("", "idx", "value", xaxis, templates.get(aa));
            is_out.addSlice(template_legend.get(aa), p.getProcessor());
        }

        ImagePlus img_out = new ImagePlus("templates", is_out);
        img_out.show();

    }

	public float[] getTemplate(int template_index) {   // for visualization - plotting

		return templates.get(template_index);

	}

	public Plot getTemplatePlot(int template_index) {

		float[] xaxis = new float[vector_len];
		for (int aa=0; aa<vector_len; aa++) xaxis[aa] = aa;

		Plot p = new Plot(template_legend.get(template_index), "idx", "value", xaxis, templates.get(template_index));
		return p;

	}

}