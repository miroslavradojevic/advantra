import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.FloatProcessor;

import java.util.ArrayList;

/**
 * Created by miroslav on 10-3-15.
 */
public class BayesianTracer3D {

    int                 Niterations;
    int                 Nstates;
    int                 Nsemicirc;

    public float[]              gcsstd;
    public int                  neuron_radius;
    float[][]                   tt;             // templates to be matched
    float[]                     tta;

    float[] img_vals;                           // reservoir for image values that will be sampled along template cross profile to calculate zncc()

    private static float        gcsstd_step = 0.5f;
    private static float        gcsstd_min  = 0.5f;

    // spatial prediction for one state sample - this is what one state will produce: Nsemicirc*gcsstd.length
    float                       step;
    float[][]                   p;              // 2d locations px,py
    float[][]                   v;              // 2d direction vx, vy

    float[]                     pred_w;         // prior weight assigned to predicted samples wrt. the direction divergence and radius change, Nsemicirc*gcsstd.length

    private static float    solidAngleRes = 1f;    // sampling resolution along arc and also when sampling the template

    public float[][][]              xt;             // bayesian filtering states (x,y, vx, vy, gcsstd)
    public float[][]                wt;             // bayesian filtering weights
    public int                      iter_counter;   // number of iterations carried out before stopped
    public float[][]                xc;             // centroid locations
    public float[]                  rc;             // centroid radiuses
    public int                      last_queue_element_checked; // remembers at which iteration the last

    // sequential bayesian filtering
    float[][]           trans_xyz;
    float[]             pties;
    public float[]      lhoods;

    private static float        cross_profile_step = 0.7f;

    private int mode = 1; // 0 - RAW IMAGE VALUES, 1 - ZNCC

    float[][] rot = new float[3][3]; // for rotations
    float[] ex = new float[3];
    float[] ey = new float[3];

    public static int getSemiSphereNrPoints(float _radius) { // count how many solid angles can fit in half sphere area, solid angle area is defined with arcRes
        return (int) Math.round((2*Math.PI*_radius*_radius)/solidAngleRes);
    }

    public BayesianTracer3D(int _Ni, int _Ns, float _step, float _neuron_radius){ //, float _sigma_deg

        Niterations = _Ni;
        Nstates = _Ns;

        xt = new float[Niterations][Nstates][7]; // final trace outputs
        wt = new float[Niterations][Nstates];
        xc = new float[Niterations][3];
        rc = new float[Niterations];
        iter_counter = Integer.MIN_VALUE;
        last_queue_element_checked = Integer.MIN_VALUE;

        neuron_radius = (int) Math.ceil(_neuron_radius);
        neuron_radius = (neuron_radius<2)? 2 : neuron_radius;

        // gcsstd define
        int cnt = 0;
        for (float sg = gcsstd_min; sg <= neuron_radius/2; sg+=gcsstd_step) cnt++;
        gcsstd = new float[cnt];
        cnt = 0;
        for (float sg = gcsstd_min; sg <= neuron_radius/2; sg+=gcsstd_step) gcsstd[cnt++] = sg;

        // tt, tta templates define
        int ns2 = (int) Math.ceil(neuron_radius/cross_profile_step);
        int ns = ns2*2+1;
        // now they are 2d for 3d tracer
        tt = new float[gcsstd.length][ns*ns]; // 2d template of the cross-section
        tta = new float[gcsstd.length];

        for (int i = 0; i < gcsstd.length; i++) {

            float ag = 0;

            for (int j = 0; j < tt[i].length; j++) {
                int u       = j % ns;
                int w       = j / ns;
                float dst   = (float) Math.sqrt(Math.pow(u-ns2,2)+Math.pow(w-ns2,2));
                tt[i][j] = (float) ((1/(Math.sqrt(2*Math.PI)*gcsstd[i])) * Math.exp(-Math.pow(dst*cross_profile_step, 2) / (2 * Math.pow(gcsstd[i], 2))));
                ag += tt[i][j];
            }

            tta[i] = ag/tt[i].length;

        }

        img_vals = new float[tt[0].length]; // take ns*ns values when calculating the zncc

        step = _step; // for the predictions of one particular state
        Nsemicirc = -1; // initially set to -1 but changed as pred() are called

        p       = null;// new float[Nsemicirc][2];
        v       = null;// new float[Nsemicirc][2];
        pred_w  = null;// new float[Nsemicirc*gcsstd.length];

        // auxiliary variables assigned with values at each recursion iteration
        trans_xyz    = null;// new float[1*Nsemicirc*gcsstd.length][4]; // becomes Nstates*getSemiSphereNrPoints()*gcsstd.length
        pties       = null;// new float[1*Nsemicirc*gcsstd.length];
        lhoods      = null;// new float[1*Nsemicirc*gcsstd.length];

    }

    // first iteration, afterwards becomes getSemiSphereNrPoints(), set so that there is enough to fill Nstates up
    public void predict(float px, float py, float pz,
                        float vx, float vy, float vz,
                        float gcsstd_pv,
                        float zDist,
                        float angula_diff_std_deg,
                        float gcsstd_diff_std_pix)
    {

        // will need a rotation
        float[][] rot = new float[3][3];
        rotation_matrix(0 ,0 , 1,   vx, vy, vz,  rot);
        float[] temp = new float[3];

        // no allocation, use available Nsemicirc value that's been set there already
        // will set the array elements in p, v, gcsstd (predictions), and priors (non-normalized prior weight assigned to each prediction)
        // values are sampled on semi sphere using saaf & kuijlaars method for uniform sampling on the sphere
        double h_k, theta_k, phi_k, phi_k_1 = 0;

        for (int k = 0; k < Nsemicirc; k++) {

            h_k = (double)k/(Nsemicirc-1); // 0 : 1

            theta_k = Math.acos(h_k);

            if(k==0 || k==(Nsemicirc-1)){

                phi_k   = 0;
                phi_k_1 = 0;

            }
            else{

                phi_k = phi_k_1 + 3.6 / ( Math.sqrt(Nsemicirc) * Math.sqrt(1-h_k*h_k));
                phi_k_1 = phi_k;

            }

            // this way top value coordinate is (x,y,z) = [0,0,1]
            p[k][0] = (float) (step * Math.sin(theta_k) * Math.cos(phi_k));
            p[k][1] = (float) (step * Math.sin(theta_k) * Math.sin(phi_k));
            p[k][2] = (float) (step * Math.cos(theta_k));

            // align towards the direction defined by vx, vy, vz, using rot matrix
            rotation_apply(rot, p[k][0], p[k][1], p[k][2], temp);

            // assign it back, calcluate rotated values
            p[k][0] = temp[0];
            p[k][1] = temp[1];
            p[k][2] = temp[2];

            // scale z dimension
            p[k][2] /= zDist;

            float norm = (float) Math.sqrt(Math.pow(p[k][0],2) + Math.pow(p[k][1],2) + Math.pow(p[k][2],2));
            v[k][0] = p[k][0] / norm;
            v[k][1] = p[k][1] / norm;
            v[k][2] = p[k][2] / norm;

            // shift predicted locations
            p[k][0] = px + p[k][0];
            p[k][1] = py + p[k][1];
            p[k][2] = pz + p[k][2];

            float ang_diff = (float) Math.acos(v[k][0]*vx+v[k][1]*vy+v[k][2]*vz); // both have norm 1, not necessary to normalize (vx, vy, vz needs to be unit vector)

            for (int j = 0; j < gcsstd.length; j++) {
                pred_w[k * gcsstd.length + j] =
                        (float) (Math.exp(-Math.pow(ang_diff,2)            / (2*Math.pow((angula_diff_std_deg/180f)*Math.PI,2))) *
                                 Math.exp(-Math.pow(gcsstd[j]-gcsstd_pv,2) / (2*Math.pow(gcsstd_diff_std_pix,   2))));
            }

        }

    }

    public void predict(float px, float py, float pz,
                        float vx, float vy, float vz,
                        float gcsstd_pv,
                        float zDist,
                        int _Nsemicirc,
                        float angula_diff_std_deg,
                        float gcsstd_diff_std_pix)
    {

        Nsemicirc = _Nsemicirc;
        p       = new float[Nsemicirc][2];
        v       = new float[Nsemicirc][2];
        pred_w  = new float[Nsemicirc*gcsstd.length];

        predict(px, py, pz, vx, vy, vz, zDist, gcsstd_pv, angula_diff_std_deg, gcsstd_diff_std_pix);

    }

    public ImagePlus getTemplates(){

        // show templates tt used for measurement
        int W = (int) Math.sqrt(tt[0].length); // length is the same for all sigmas
        int H = W;

        ImageStack is = new ImageStack(W, H);

        for (int i = 0; i < tt.length; i++)
            is.addSlice("gcsstd=" + IJ.d2s(gcsstd[i], 1), new FloatProcessor(W, H, tt[i].clone()));

        ImagePlus im_templates = new ImagePlus("neuron radius = "+IJ.d2s(neuron_radius,0), is);

        return im_templates;

    }

    public void reset() {

        // reset output xt, wt values
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

    public int trace(float x,  float y,  float z,
                     float vx, float vy, float vz,
                     float gcsstd,
                     float zDist,
                     float[][][] _inimg_xyz,
                     float _angula_deg,
                     float _gcsstd_pix,
                     ArrayList<int[]>           _tag_map,
                     ArrayList<boolean[][]>     _mask_xy,
                     ArrayList<boolean[][]>     _queue_map)
    {

        int xnode, ynode, znode;
        int W = _inimg_xyz.length;
        int H = _inimg_xyz[0].length;
        int L = _inimg_xyz[0][0].length;

        reset(); // before each trace, xt, wt are initialized with zeros if there were any from

        iter_counter = 0;
        last_queue_element_checked = 0; // because we start from one

        while (iter_counter < Niterations) { // till the limit

            if (iter_counter == 0)
                iter(x, y, z, vx, vy, vz, gcsstd, zDist, _inimg_xyz, _angula_deg, _gcsstd_pix); // it is different the first time
            else
                iter(iter_counter, zDist, _inimg_xyz, _angula_deg, _gcsstd_pix);

            // round them and check
            xnode = Math.round(xc[iter_counter][0]);
            ynode = Math.round(xc[iter_counter][1]);
            znode = Math.round(xc[iter_counter][2]);

            if (_queue_map.get(znode)[xnode][ynode] && iter_counter>1) last_queue_element_checked = iter_counter;

            if (!_mask_xy.get(znode)[xnode][ynode]) return 0;
            if (_tag_map.get(znode)[ynode * W + xnode] != 0) return _tag_map.get(znode)[ynode * W + xnode];

            iter_counter++; // will be used to describe how far we've gone

        }

        return 0;

    }


    private void iter(
            int                     _iter,
            float                   _zDist,
            float[][][]             _inimg_xyz,
            float                   _angula_deg,
            float                   _gcsstd_pix
    )
    {

        // this one makes sense for _iter>0, there is part that is executed if _iter==1
        // take states at _iter-1
        float[][] prev_x = xt[_iter-1];
        float[]   prev_w = wt[_iter-1];

        int L = tt[0].length; // length is the same for all sigmas
        int L2  = L/2;
        int imgW = _inimg_xyz.length;
        int imgH = _inimg_xyz[0].length;
        int imgL = _inimg_xyz[0][0].length;

        int count = 0;
        float sum_lhoods = 0;// used to check validity
        float sum_prior = 0;

        for (int i = 0; i < prev_x.length; i++) { // will loop throughout the states Nstates

            float px        = prev_x[i][0];
            float py        = prev_x[i][1];
            float pz        = prev_x[i][2];

            float vx        = prev_x[i][3];
            float vy        = prev_x[i][4];
            float vz        = prev_x[i][5];

            float gcsstdc   = prev_x[i][6];

            if (i==0 & _iter==1) {

                predict(px, py, pz, vx, vy, vz, gcsstdc, _zDist, getSemiSphereNrPoints(neuron_radius), _angula_deg, _gcsstd_pix);

                trans_xyz   = new float[Nstates*pred_w.length][7];   // to store all the predictions at current iteration
                pties       = new float[Nstates*pred_w.length];      // store the probabilities temporarily
                lhoods      = new float[Nstates*pred_w.length];      // store the likelihoods

            }
            else
                predict(px, py, pz, vx, vy, vz, gcsstdc, _zDist, _angula_deg, _gcsstd_pix);

            // the reminder is the same as in the other version of iter()
            for (int j = 0; j < pred_w.length; j++) { // loop all options, calculate likelihoods

                trans_xyz[count][0] = p[j/gcsstd.length][0]; // px
                trans_xyz[count][1] = p[j/gcsstd.length][1]; // py
                trans_xyz[count][2] = p[j/gcsstd.length][2]; // pz

                trans_xyz[count][3] = v[j/gcsstd.length][0]; // vx
                trans_xyz[count][4] = v[j/gcsstd.length][1]; // vy
                trans_xyz[count][5] = v[j/gcsstd.length][2]; // vz

                trans_xyz[count][6] = gcsstd[j%gcsstd.length]; // gcsstd[]

                lhoods[count] = 0; // will take the highest calculated correlation out of all possible correlations











                // check if this sampling range falls into the image at all - borders are enough to check, before picking the values from the image
                // ex, ey, ez, where ez is aligned with vx0, vy0, vz0
                float bdry = (float) (Math.sqrt(2) * L2 * cross_profile_step);

                if (
                        trans_xyz[count][0]>=0+bdry && trans_xyz[count][0]<=imgW-1-bdry &&
                        trans_xyz[count][1]>=0+bdry && trans_xyz[count][1]<=imgH-1-bdry &&
                        trans_xyz[count][2]>=0+bdry/_zDist && trans_xyz[count][2]<=imgL-1-(bdry/_zDist)
                        )
                    if (mode==0)
                        lhoods[count] = getVal(trans_xyz[count][0], trans_xyz[count][1], trans_xyz[count][2], _inimg_xyz);
                    else if (mode == 1)
                        lhoods[count] = zncc(
                                trans_xyz[count][0], trans_xyz[count][1], trans_xyz[count][2],
                                trans_xyz[count][3], trans_xyz[count][4], trans_xyz[count][5], j % gcsstd.length, _inimg_xyz);
                    else
                        System.exit(-1);






//                // check if this sampling range falls into the image at all - borders are enough to check
//                float p1x = trans_xy[count][0] + L2*cross_profile_step*trans_xy[count][3]; // x+L*max_step*vy
//                float p1y = trans_xy[count][1] - L2*cross_profile_step*trans_xy[count][2]; // y-L*max_step*vx
//                float p2x = trans_xy[count][0] - L2*cross_profile_step*trans_xy[count][3]; // x-L*max_step*vy
//                float p2y = trans_xy[count][1] + L2*cross_profile_step*trans_xy[count][2]; // y+L*max_step*vx
//
//                if ((p1x>=0 && p1x<=imgW-1 && p1y>=0 && p1y<=imgH-1) && (p2x>=0 && p2x<=imgW-1 && p2y>=0 && p2y<=imgH-1)) // both limits fall within the image
//                    if (mode==0)
//                        lhoods[count] = getVal(trans_xy[count][0], trans_xy[count][1], _inimg_xy);
//                    else if (mode==1)
//                        lhoods[count] = zncc(trans_xy[count][0], trans_xy[count][1], trans_xy[count][2], trans_xy[count][3], j % gcsstd.length, _inimg_xy);
//                    else
//                        System.exit(-1);





                sum_lhoods += lhoods[count];

                pties[count] = pred_w[j]; // w stores priors, possible to immediately multiply posterior here unless it is all zeros
                sum_prior += pties[count];

                count++;
            }

        }

        // just in case check - if all the lhoods were zero then stay with the priors
        if (sum_lhoods>Float.MIN_VALUE)
            for (int i = 0; i < pties.length; i++)
                pties[i] = prev_w[i/pred_w.length] * (pties[i]/sum_prior) * lhoods[i];

        // take best Ns (we'll always have enough to take)
        int[] sort_idx = Toolbox.descending(pties);      // ptes will be sorted as a side effect

        int cnt = 0;
        for (int i = 0; i < sort_idx.length; i++) {

            xt[_iter][cnt][0] = trans_xyz[sort_idx[i]][0];
            xt[_iter][cnt][1] = trans_xyz[sort_idx[i]][1];
            xt[_iter][cnt][2] = trans_xyz[sort_idx[i]][2];
            xt[_iter][cnt][3] = trans_xyz[sort_idx[i]][3];
            xt[_iter][cnt][4] = trans_xyz[sort_idx[i]][4];
            xt[_iter][cnt][5] = trans_xyz[sort_idx[i]][5];
            xt[_iter][cnt][6] = trans_xyz[sort_idx[i]][6];

            wt[_iter][cnt]    = pties[i]; // because they are already sorted
            cnt++;

            if (cnt==Nstates) break; // add top Ns points that are in the circle for further tracing

        }

        Toolbox.probability_distribution(wt[_iter]); // _wt will be normalized

        // set xc, rc at the end of the iteration
        xc[_iter][0] = 0;
        xc[_iter][1] = 0;
        xc[_iter][2] = 0;
        rc[_iter]    = 0;

        for (int j = 0; j < xt[_iter].length; j++) {
            xc[_iter][0]    += wt[_iter][j] * xt[_iter][j][0];
            xc[_iter][1]    += wt[_iter][j] * xt[_iter][j][1];
            xc[_iter][2]    += wt[_iter][j] * xt[_iter][j][2];
            rc[_iter]       += wt[_iter][j] * xt[_iter][j][6];
        }

    }

    private void iter(
            float                   x0,
            float                   y0,
            float                   z0,
            float                   vx0,
            float                   vy0,
            float                   vz0,
            float                   gcsstdc0,
            float                   _zDist,
            float[][][]             _inimg_xyz,
            float                   _angula_deg,
            float                   _gcsstd_pix
    )
    {

        // normally there is a list of states, but now there is only one so this step will be skipped
        // auxilliary variable used for rotations at each 3d state location+direction
        int L = (int) Math.sqrt(tt[0].length); // length is the same for all sigmas
        int L2 = L/2;
        int imgW = _inimg_xyz.length;
        int imgH = _inimg_xyz[0].length;
        int imgL = _inimg_xyz[0][0].length;

        predict(x0, y0, z0, vx0, vy0, vz0, gcsstdc0, _zDist, Nstates, _angula_deg, _gcsstd_pix); // predict so that number of elements at the semicirc is set (will do allocation of p,v,pred_w inside)

        // some aux values for bayesian sequential filtering (Nsemicirc is set to Nstates in predict())
        trans_xyz   = new float[1*pred_w.length][7];   // to store all the predictions at current iteration
        pties       = new float[1*pred_w.length];      // store the probabilities temporarily
        lhoods      = new float[1*pred_w.length];      // store the likelihoods

        int count = 0;
        float sum_lhoods = 0;// used to check validity
        float sum_prior = 0;

        for (int j = 0; j < pred_w.length; j++) { // loop all options, calculate likelihoods

            trans_xyz[count][0] = p[j/gcsstd.length][0]; // px
            trans_xyz[count][1] = p[j/gcsstd.length][1]; // py
            trans_xyz[count][2] = p[j/gcsstd.length][2]; // pz

            trans_xyz[count][3] = v[j/gcsstd.length][0]; // vx
            trans_xyz[count][4] = v[j/gcsstd.length][1]; // vy
            trans_xyz[count][5] = v[j/gcsstd.length][2]; // vy

            trans_xyz[count][6] = gcsstd[j%gcsstd.length]; // gcsstd[]

            lhoods[count] = 0; // will take the highest calculated correlation out of all possible correlations

            // check if this sampling range falls into the image at all - borders are enough to check, before picking the values from the image
            // ex, ey, ez, where ez is aligned with vx0, vy0, vz0
            float bdry = (float) (Math.sqrt(2) * L2 * cross_profile_step);

            if (
                    trans_xyz[count][0]>=0+bdry && trans_xyz[count][0]<=imgW-1-bdry &&
                    trans_xyz[count][1]>=0+bdry && trans_xyz[count][1]<=imgH-1-bdry &&
                    trans_xyz[count][2]>=0+bdry/_zDist && trans_xyz[count][2]<=imgL-1-(bdry/_zDist)
            )
                if (mode==0)
                    lhoods[count] = getVal(trans_xyz[count][0], trans_xyz[count][1], trans_xyz[count][2], _inimg_xyz);
                else if (mode == 1)
                    lhoods[count] = zncc(
                            trans_xyz[count][0], trans_xyz[count][1], trans_xyz[count][2],
                            trans_xyz[count][3], trans_xyz[count][4], trans_xyz[count][5], j % gcsstd.length, _inimg_xyz);
                else
                    System.exit(-1);

            sum_lhoods += lhoods[count];

            pties[count] = pred_w[j]; // w stores priors, possible to immediately multiply posterior here unless it is all zeros
            sum_prior += pties[count];

            count++;
        }

        // just in case check - if all the lhoods were zero then stay with the priors
        if (sum_lhoods>Float.MIN_VALUE)
            for (int i = 0; i < pties.length; i++)
                pties[i] = 1 * (pties[i]/sum_prior) * lhoods[i];

        // take best Ns (we'll always have enough to take)
        int[] sort_idx = Toolbox.descending(pties);       // ptes will be sorted as a side effect

        int cnt = 0;
        for (int i = 0; i < sort_idx.length; i++) {

            xt[0][cnt][0] = trans_xyz[sort_idx[i]][0];
            xt[0][cnt][1] = trans_xyz[sort_idx[i]][1];
            xt[0][cnt][2] = trans_xyz[sort_idx[i]][2];
            xt[0][cnt][3] = trans_xyz[sort_idx[i]][3];
            xt[0][cnt][4] = trans_xyz[sort_idx[i]][4];
            xt[0][cnt][5] = trans_xyz[sort_idx[i]][5];
            xt[0][cnt][6] = trans_xyz[sort_idx[i]][6];

            wt[0][cnt]    = pties[i];                   // because they are already sorted
            cnt++;

            if (cnt==Nstates) break;                    // add top Ns points that are in the circle for further tracing

        }

        Toolbox.probability_distribution(wt[0]);           // _wt will be normalized

        // set xc, rc at the end of the iteration
        xc[0][0] = 0;
        xc[0][1] = 0;
        xc[0][2] = 0;
        rc[0]    = 0;

        for (int j = 0; j < xt[0].length; j++) {
            xc[0][0]    += wt[0][j] * xt[0][j][0];
            xc[0][1]    += wt[0][j] * xt[0][j][1];
            xc[0][2]    += wt[0][j] * xt[0][j][2];
            rc[0]       += wt[0][j] * xt[0][j][6];
        }

    }



    private float getVal(float px, float py, float pz, float[][][] _inimg_xyz)
    {
        return Interpolator.interpolateAt(px, py, pz, _inimg_xyz);
    }

    private float zncc(float px, float py, float pz, float vx, float vy, float vz, int gcsstd_idx, float[][][] _inimg_xyz)
    {

        rotation_matrix(0,0,1, vx, vy, vz, rot);
        rotation_apply(rot, 1,0,0, ex);
        rotation_apply(rot, 0,1,0, ey);

        int L = (int) Math.sqrt(tt[0].length); // length is the same for all sigmas
        int L2 = L/2;
        float bdry = L2*cross_profile_step;

        float ag = 0;

        for (int j = 0; j < tt[gcsstd_idx].length; j++) {

            int u       = j % L;
            int w       = j / L;

            float at_x = (u-L2)*ex[0] + (w-L2)*ey[0];
            float at_y = (u-L2)*ex[1] + (w-L2)*ey[1];
            float at_z = (u-L2)*ex[2] + (w-L2)*ey[2];

            img_vals[j] = Interpolator.interpolateAt(at_x, at_y, at_z, _inimg_xyz);
            ag += img_vals[j];

        }

        ag = ag / (float)(img_vals.length);

        // calculate zncc, use template that corresponds to this sigma
        float corra = 0;
        float corrb = 0;
        float corrc = 0;

        for (int profileidx = 0; profileidx <img_vals.length; profileidx++) {
            corra += (img_vals[profileidx]-ag) * (tt[gcsstd_idx][profileidx]-tta[gcsstd_idx]);
            corrb += Math.pow(img_vals[profileidx]-ag, 2);
            corrc += Math.pow(tt[gcsstd_idx][profileidx]-tta[gcsstd_idx], 2);
        }

        return  (float) (corra / Math.sqrt(corrb*corrc));

    }

    private void rotation_matrix(
            float a1, float a2, float a3,
            float b1, float b2, float b3,
            float[][] R
    )
    {

        // v is cross product of (a1, a2, a3) and (b1, b2, b3)
        float v1 = a2*b3 - b2*a3;
        float v2 = -(a1*b3-b1*a3);
        float v3 = a1*b2-b1*a2;

        float tt = (1-(a1*b1+a2*b2+a3*b3))/(v1*v1+v2*v2+v3*v3);

        // from http://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
        R[0][0] = 1 + 0     + tt * (-v3*v3-v2*v2);
        R[0][1] = 0 + (-v3) + tt * (v1*v2);
        R[0][2] = 0 + (v2)  + tt * (-v1*v3);

        R[1][0] = 0 + (v3)  + tt * (v1*v2);
        R[1][1] = 1 + 0     + tt * (-v3*v3-v1*v1);
        R[1][2] = 0 + (-v1) + tt * (v2*v3);

        R[2][0] = 0 + (-v2) + tt * (v1*v3);
        R[2][1] = 0 + (v1)  + tt * (v2*v3);
        R[2][2] = 1 + 0     + tt * (-v2*v2-v1*v1);

    }

    private void rotation_apply(
            float[][] R,
            float v1, float v2, float v3,
            float[] Rvec
    ){

        Rvec[0] = R[0][0]*v1 + R[0][1]*v2 + R[0][2]*v3;
        Rvec[1] = R[1][0]*v1 + R[1][1]*v2 + R[1][2]*v3;
        Rvec[2] = R[2][0]*v1 + R[2][1]*v2 + R[2][2]*v3;

    }

    public ArrayList<float[]> fullSpherePts(int N){

        ArrayList<float[]> points = new ArrayList<float[]>(N);// float[N][3];

        double h_k, theta_k, phi_k, phi_k_1 = 0;

        for (int k = 0; k < N; k++) {

            h_k = -1 + 2 * (double)k/(N-1); // -1 : 1

            theta_k = Math.acos(h_k);

            if(k==0 || k==(N-1)){

                phi_k   = 0;
                phi_k_1 = 0;

            }
            else{

                phi_k = phi_k_1 + 3.6 / ( Math.sqrt(N) * Math.sqrt(1-h_k*h_k));
                phi_k_1 = phi_k;

            }


            float px = (float) (Math.sin(theta_k) * Math.cos(phi_k));
            float py = (float) (Math.sin(theta_k) * Math.sin(phi_k));
            float pz = (float) Math.cos(theta_k);

            points.add(new float[]{px, py, pz});

        }

        return points;

    }

    public ArrayList<float[]> fullSpherePts(float radius, float solid_angle){

        // N is determined so that number of solid angles is fit onto the sphere
        int N = (int) Math.round((4 * Math.PI * radius * radius) / solid_angle);

        ArrayList<float[]> points = new ArrayList<float[]>(N);// float[N][3];

        double h_k, theta_k, phi_k, phi_k_1 = 0;

        for (int k = 0; k < N; k++) {

            h_k = -1 + 2 * (double)k/(N-1); // -1 : 1

            theta_k = Math.acos(h_k);

            if(k==0 || k==(N-1)){

                phi_k   = 0;
                phi_k_1 = 0;

            }
            else{

                phi_k = phi_k_1 + 3.6 / ( Math.sqrt(N) * Math.sqrt(1-h_k*h_k));
                phi_k_1 = phi_k;

            }

            float px = radius * (float) (Math.sin(theta_k) * Math.cos(phi_k));
            float py = radius * (float) (Math.sin(theta_k) * Math.sin(phi_k));
            float pz = radius * (float)  Math.cos(theta_k);

            points.add(new float[]{px, py, pz});

        }

        return points;

    }

}
