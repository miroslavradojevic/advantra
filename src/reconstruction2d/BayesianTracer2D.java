package reconstruction2d;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;

/**
 * Created by miroslav on 13-2-15.
 */
public class BayesianTracer2D {

    float[][][]     xt;//      = new float[Ni+1][_Ns][4];    // bayesian filtering states (x,y, vx, vy, r)
    float[][]       wt;//      = new float[Ni+1][_Ns];       // bayesian filtering weights

    float       neuron_radius;
    float       sigma_deg;

    float[][]   tt;                                          // templates to be matched when filtering 2*Rcnt+1
    float[]     tta;
    float[]     sigmas;

    private static float        cross_profile_step = 0.5f;
    private static float        sigma_step = 0.5f;

    public BayesianTracer2D(int _Ni, int _Ns, float _neuron_radius, float _sigma_deg){

        xt = new float[_Ni][_Ns][5];
        wt = new float[_Ni][_Ns];

        neuron_radius = _neuron_radius;

        sigma_deg = _sigma_deg;

        int _neuron_radius_round = Math.round(neuron_radius);
        _neuron_radius_round = (_neuron_radius_round<1)? 1 : _neuron_radius_round;


        // sigmas define
        float sg_min = 0.5f;
        int cnt = 0;
        for (float sg = sg_min; sg <= _neuron_radius_round; sg+=sigma_step) {
            cnt++;
        }
        System.out.println(cnt);
        sigmas = new float[cnt];

        cnt = 0;
        for (float sg = sg_min; sg <= _neuron_radius_round; sg+=sigma_step) {
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
    
    public void showTemplates(){

        // show templates tt

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

        new ImagePlus("neuron radius = "+IJ.d2s(neuron_radius,1), is).show();

    }

    private void getLikelihoods(float _x, float _y, float[][] _inimg_xy) {
        // use tt templates

    }

    public void reset(int _Ni, int _Ns) {
    }

    public void trace () {

        // before the trace, the values are initialized with zeros if there were any from
        // earlier function call


    }

//    public


}
