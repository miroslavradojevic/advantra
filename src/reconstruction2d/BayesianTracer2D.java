package reconstruction2d;

import detection2d.SemiCircle;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;

/**
 * Created by miroslav on 13-2-15.
 */
public class BayesianTracer2D {

    int Niterations;
    int Nstates;

    float[][][]     xt;//      = new float[Ni+1][_Ns][4];    // bayesian filtering states (x,y, vx, vy, r)
    float[][]       wt;//      = new float[Ni+1][_Ns];       // bayesian filtering weights

    // auxiliary arrays for sequential bayesian filtering
    float[][]   trans_xy;
    float[]     pties;
    float[]     lhoods;

    public int          neuron_radius;
    public float       step;
    public float       sigma_deg;

    float[][]   tt;                                          // templates to be matched when filtering 2*Rcnt+1
    float[]     tta;
    float[]     sigmas;

    SemiCircle scirc;

    private static float        cross_profile_step = 0.5f;
    private static float        sg_min = 1.0f;
    private static float        sigma_step = 0.5f;

    public BayesianTracer2D(int _Ni, int _Ns, float _step, float _neuron_radius, float _sigma_deg){

        Niterations = _Ni;
        Nstates = _Ns;

        step = _step;
        scirc = new SemiCircle(step);
        int total_preds = _Ns * scirc.NN;
        System.out.println(total_preds + "  total predictions");

        xt = new float[_Ni][_Ns][5];
        wt = new float[_Ni][_Ns];

        // auxiliary variables assigned with values at each iteration
        trans_xy    = new float[total_preds][4];
        pties       = new float[total_preds];
        lhoods      = new float[total_preds];

        sigma_deg = _sigma_deg;

        int neuron_radius = Math.round(_neuron_radius);
        neuron_radius = (neuron_radius<1)? 1 : neuron_radius;

        // sigmas define

        int cnt = 0;
        for (float sg = sg_min; sg <= neuron_radius; sg+=sigma_step) {
            cnt++;
        }
        sigmas = new float[cnt];

        cnt = 0;
        for (float sg = sg_min; sg <= neuron_radius; sg+=sigma_step) {
            sigmas[cnt++] = sg;
        }

        // templates define
        int ns2 = (int) Math.ceil(sigmas[sigmas.length - 1] / cross_profile_step);
        int ns = ns2*2+1;

        tt = new float[sigmas.length][ns];
        tta = new float[sigmas.length];

        for (int i = 0; i < sigmas.length; i++) {

            float ag = 0;

            for (int j = 0; j < ns; j++) {
                tt[i][j] = (float) ((1/(Math.sqrt(2*Math.PI)*sigmas[i])) * Math.exp(-Math.pow(j-ns2,2)/(2*Math.pow(sigmas[i],2))));
                ag += tt[i][j];
            }

            tta[i] = ag/ns;
        }



    }
    
    public ImagePlus getTemplates(){

        // show templates tt used for measurement

        int L = tt[0].length; // length is the same for all sigmas
        int L2  = L/2;

        ImageStack is = new ImageStack(528, 255);

        float[] yy = new float[L];
        float[] xx = new float[L];

        for (int i = 0; i < sigmas.length; i++) {


            for (int j = -L2; j <= L2; j++) {
                yy[j+L2] = tt[i][j+L2];
                xx[j+L2] = (j+L2) * cross_profile_step;
            }

            Plot p = new Plot("", "cs[pix]", "t", xx, yy, Plot.LINE);

            is.addSlice("sig="+ IJ.d2s(sigmas[i],1), p.getProcessor());

        }

        ImagePlus im_templates = new ImagePlus("neuron radius = "+IJ.d2s(neuron_radius,1), is);

        return im_templates;

    }

    private void getLikelihoods(float _x, float _y, float[][] _inimg_xy) {
        // use tt templates

    }

    public void reset() {

        // reset xt, wt values
        for (int i = 0; i < xt.length; i++) {
            for (int j = 0; j < xt[i].length; j++) {
                for (int k = 0; k < xt[i][j].length; k++) {
                    xt[i][j][k] = 0;
                }
            }
        }

        for (int i = 0; i < wt.length; i++) {
            for (int j = 0; j < wt[i].length; j++) {
                wt[i][j] = 0;
            }
        }

    }

    public void trace (float x, float y, float vx) {

        // before the trace, the values are initialized with zeros if there were any from
        // earlier function call
        reset();

        int iter = 0;

        while (iter < Niterations) {

            if (iter == 0) {
                bayesian_iteration();
            }
            else {
                bayesian_iteration(at_x, at_y, img_xy,iter,semi_circ_step, sigma_deg, xt,wt,pt,ft,mt);
            }




            iter++;
        }

    }



}
