package tracing2d;

import java.util.ArrayList;

/**
 * Created by miroslav on 31-8-14.
 */
public class Sphere2D {

    private static float    arcRes 	        = 3.0f;
    public static float 	arcNbhood       = 2*arcRes;
    private static float    samplingStep    = 0.9f;

    public static float     T_HALF 	        = 0.50f;

    private float 	radius;
    private float   neuronDiameter;
    private float 	sigma_ratio;
    public int     N;

    private int 	limR, limT;

    private static float TWO_PI = (float) (2 * Math.PI);
    private static float ONE_PI = (float) (1 * Math.PI);

    // variables for bayesian tracking
    public static ArrayList<Float>          theta = new ArrayList<Float>(); 	        // list of elements (theta) covering the circle
    public static ArrayList<float[]>        locsXY  = new ArrayList<float[]>();      // list of follow up location offsets (for prediction transistion)
    public static ArrayList<float[]>        vxy     = new ArrayList<float[]>();          // unit vector directions corresponding to thetas


    public Sphere2D(float neuronDiam, float scale)
    {
        this.radius 	= scale * neuronDiam;
        this.neuronDiameter = neuronDiam;

        this.N 	= (int) Math.ceil( ( (2 * Math.PI * radius) / arcRes) );    // N will influence theta list size, and offstXY list size
        this.limT = (int) Math.ceil(T_HALF*neuronDiameter/samplingStep);    // transversal sampling limits
        this.limR = 2 * limT + 1; // how many to take radially with given sampling step

        /*
            theta
         */
        theta.clear();
        for (int i=0; i<N; i++) {
            theta.add( i * ( (float)(2*Math.PI) / N ) );
        }

        /*
            prediction locations locXY, and the directions vxy (for prior)
         */
        locsXY.clear();
        vxy.clear();
        for (int ii = 0; ii < theta.size(); ii++) {

            float curr_theta = theta.get(ii);

            int i = 0;
            int k = - (int) Math.round(limR * (1/2f));// here it scales the step forward

            float px = i * samplingStep;
            float py = k * samplingStep;

            float[][] locPerDirection = new float[1][2];

            locPerDirection[0][0] = px;
            locPerDirection[0][1] = py;

            transY(radius, locPerDirection);
            rotZ(curr_theta, locPerDirection);
            locsXY.add(locPerDirection[0]);

            float vx = 0;
            float vy = 0;

            float[][] vecPerDirection = new float[1][2];

            vecPerDirection[0][0] = vx;
            vecPerDirection[0][1] = vy;

            transY(1, vecPerDirection);
            rotZ(curr_theta, vecPerDirection);
            // normalize vxy (just in case)
            double vnorm = Math.sqrt( vecPerDirection[0][0]*vecPerDirection[0][0] + vecPerDirection[0][1]*vecPerDirection[0][1] );
            vecPerDirection[0][0] /= vnorm;
            vecPerDirection[0][1] /= vnorm;

            vxy.add(vecPerDirection[0]);

        }



    }

    private void rotZ(float ang, float[][] coords)
    {
        for (int i=0; i<coords.length; i++) {
            float x_temp = coords[i][0];
            float y_temp = coords[i][1];
            coords[i][0] = x_temp * (float) Math.cos(ang) - y_temp * (float) Math.sin(ang);
            coords[i][1] = x_temp * (float) Math.sin(ang) + y_temp * (float) Math.cos(ang);
        }
    }

    private void transY(float ty, float[][] coords)
    {
        for (int i=0; i<coords.length; i++){
            coords[i][1] += ty;
        }
    }

    public void getPriors(float[] _vxy, float _deg_std, float[] _store_priors)
    {

        // normalize the input direction
        double _vnorm = Math.sqrt( _vxy[0]*_vxy[0] + _vxy[1]*_vxy[1] );
        _vxy[0] /= _vnorm;
        _vxy[1] /= _vnorm;

        // take the reference direction and assign the probability to each of the transition points
        float sum_priors = 0;
        for (int i = 0; i < vxy.size(); i++) {

            float dot_prod = _vxy[0]*vxy.get(i)[0] + _vxy[1]*vxy.get(i)[1];
            dot_prod = (dot_prod>1)? 1 : dot_prod;
            double ang_rad = Math.acos(dot_prod);
            double ang_deg = ang_rad*(180f/Math.PI);

            // form gaussian weighted prior
            _store_priors[i] = (float) (  Math.exp( -((ang_deg-0)*(ang_deg-0)) / (2*_deg_std*_deg_std) )); // (float) (1 / (_deg_std*Math.sqrt(2*Math.PI)) ) *
            sum_priors += _store_priors[i];

        }

        for (int i = 0; i < vxy.size(); i++) {
            _store_priors[i] /= sum_priors;
        }

    }

}
