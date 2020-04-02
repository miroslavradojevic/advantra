package reconstruction;

import aux.Interpolator;
import aux.Stat;
import aux.Tools;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.*;

import java.awt.*;

/**
 * Created by miroslav on 13-2-15.
 */
public class BayesianTracer2D {

    int                         Niterations;    // number of iterations to allocate (max nr. iterations)
    int                         Nstates;        // number of states to maintain
    int                         Nsemicirc;      // number of states capturing semi circle

    public float[]              gcsstd;         // list of gaussian cross section standard dev. (will take discrete set of values, several scales)
    public int                  neuron_radius;
    float[][]                   tt;             // templates to be matched when filtering 2*Rcnt+1
    float[]                     tta;

    float[] img_vals;                           // reservoir for image values that will be sampled along template cross profile to calculate zncc()

    private static float        gcsstd_step = 0.5f;
    private static float        gcsstd_min  = 0.5f;

    // spatial prediction for one state sample - this is what one state will produce: Nsemicirc*gcsstd.length
    float                       step;
    float[][]                   p;              // 2d locations px,py
    float[][]                   v;              // 2d direction vx, vy
    // this will run for each of the Nstates, Niterations times, each time picking the best Nstates to maintain

    float[]              pred_w;         // prior weight assigned to predicted samples wrt. the direction divergence and radius change, Nsemicirc*gcsstd.length

    private static float    arcRes = 1f;    // sampling resolution along arc

    public float[][][]              xt;             // bayesian filtering states (x,y, vx, vy, gcsstd)
    public float[][]                wt;             // bayesian filtering weights
    public int                      iter_counter;   // number of iterations carried out before stopped
    public float[][]                xc;             // centroid locations
    public float[]                  rc;             // centroid radiuses
    public int                      last_queue_element_checked; // remembers at which iteration the last

    // sequential bayesian filtering
    float[][]           trans_xy;
    float[]             pties;
    public float[]      lhoods;

    private static float        cross_profile_step = 0.5f;

    private int mode = 1; // 0 - RAW IMAGE VALUES, 1 - ZNCC

    public static int getSemiCircleNrPoints(float _radius) {
        return (int) Math.ceil((Math.PI*_radius)/arcRes);
    }

    public BayesianTracer2D(int _Ni, int _Ns, float _step, float _neuron_radius){ //, float _sigma_deg

        Niterations = _Ni;
        Nstates = _Ns;

        xt = new float[Niterations][Nstates][5]; // final trace outputs
        wt = new float[Niterations][Nstates];
        xc = new float[Niterations][2];
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

        tt = new float[gcsstd.length][ns];
        tta = new float[gcsstd.length];

        for (int i = 0; i < gcsstd.length; i++) {

            float ag = 0;
            for (int j = 0; j < ns; j++) {
                tt[i][j] = (float) ((1/(Math.sqrt(2*Math.PI)*gcsstd[i])) * Math.exp(-Math.pow((j - ns2) * cross_profile_step, 2) / (2 * Math.pow(gcsstd[i], 2)))); //() *
                ag += tt[i][j];
            }

            tta[i] = ag/ns;
        }

        img_vals = new float[tt[0].length];

        step = _step; // for the predictions of one particular state
        Nsemicirc = -1;

        p       = null;// new float[Nsemicirc][2];
        v       = null;// new float[Nsemicirc][2];
        pred_w  = null;// new float[Nsemicirc*gcsstd.length];

        // auxiliary variables assigned with values at each recursion iteration
        trans_xy    = null;// new float[1*Nsemicirc*gcsstd.length][4]; // becomes Nstates*getArcNrPoints()*gcsstd.length
        pties       = null;// new float[1*Nsemicirc*gcsstd.length];
        lhoods      = null;// new float[1*Nsemicirc*gcsstd.length];

    }

    // first iteration, afterwards becomes getArcNrPoints(), set so that there is enough to fill Nstates up
    public void predict(float px, float py, float vx, float vy, float gcsstd_pv,
                        float angula_diff_std_deg,
                        float gcsstd_diff_std_pix)
    {
        // no allocation, use available Nsemicirc value that's been set there already
        // will set the array elements in p, v, gcsstd (predictions), and priors (non-normalized prior weight assigned to each prediction)
        for (int i = 0; i < Nsemicirc; i++) {

            float alfa = (float) (-Math.PI/2 + Math.PI * ((float)i/(Nsemicirc-1))); // -pi/2 -- +pi/2

            v[i][0] = (float) (vx * Math.cos(alfa) - vy * Math.sin(alfa)); // vx
            v[i][1] = (float) (vx * Math.sin(alfa) + vy * Math.cos(alfa)); // vy

            p[i][0] = px + step * v[i][0]; // px
            p[i][1] = py + step * v[i][1]; // py

            for (int j = 0; j < gcsstd.length; j++) {
                pred_w[i*gcsstd.length+j] =
                        (float) (Math.exp(-Math.pow(alfa,2)                    / (2*Math.pow((angula_diff_std_deg/180f)*Math.PI,2))) *
                                 Math.exp(-Math.pow(gcsstd[j]-gcsstd_pv,2)     / (2*Math.pow(gcsstd_diff_std_pix,   2))));
            }

        }

    }

    public void predict(float px, float py, float vx, float vy, float gcsstd_pv,
                        int _Nsemicirc,
                        float angula_diff_std_deg,
                        float gcsstd_diff_std_pix)
    {

        Nsemicirc = _Nsemicirc;
        p       = new float[Nsemicirc][2];
        v       = new float[Nsemicirc][2];
        pred_w  = new float[Nsemicirc*gcsstd.length];

        predict(px, py, vx, vy, gcsstd_pv, angula_diff_std_deg, gcsstd_diff_std_pix);

    }

    public ImagePlus getTemplates(){

        // show templates tt used for measurement
        int L = tt[0].length; // length is the same for all sigmas
        int L2  = L/2;

        ImageStack is = new ImageStack(528, 255);

        float[] yy = new float[L];
        float[] xx = new float[L];

        for (int i = 0; i < gcsstd.length; i++) {
            for (int j = -L2; j <= L2; j++) {
                yy[j+L2] = tt[i][j+L2];
                xx[j+L2] = (j+L2) * cross_profile_step;
            }
            Plot p = new Plot("", "cs[pix]", "t", xx, yy, Plot.LINE);
            is.addSlice("gcsstd="+ IJ.d2s(gcsstd[i],1), p.getProcessor());
        }

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

    public int trace(float x, float y, float vx, float vy, float gcsstd,
                     float[][] _inimg_xy, float _angula_deg, float _gcsstd_pix,
                     int[] _tag_map,
                     boolean[][] _mask_xy,
                     boolean[][] _queue_map)
    {

//        float[] x_y_r = new float[3];
        int xnode, ynode;
        int W = _inimg_xy.length;

        reset(); // before each trace, the values are initialized with zeros if there were any from

        iter_counter = 0;
        last_queue_element_checked = 0; // because we start from one

        while (iter_counter < Niterations) { // till the limit

            if (iter_counter == 0)
                iter(x, y, vx, vy, gcsstd, _inimg_xy, _angula_deg, _gcsstd_pix); // it is different the first time
            else
                iter(iter_counter, _inimg_xy, _angula_deg, _gcsstd_pix);

            // round them and check
            xnode = Math.round(xc[iter_counter][0]);
            ynode = Math.round(xc[iter_counter][1]);

            if (_queue_map[xnode][ynode] && iter_counter>1) last_queue_element_checked = iter_counter;

            if (!_mask_xy[xnode][ynode]) return 0;
            if (_tag_map[ynode * W + xnode] != 0) return _tag_map[ynode * W + xnode];

            iter_counter++; // will be used to describe how far we've gone

            // check
//            getNode(it, x_y_r);// will return xcentroid, ycentroid and rcentroid at this iteration
//            it++;

        }

        return 0;

    }

    private void iter(
            int                     _iter,
            float[][]               _inimg_xy,
            float                   _angula_deg,
            float                   _gcsstd_pix
    )
    {

        int L = tt[0].length; // length is the same for all sigmas
        int L2  = L/2;
        int imgW = _inimg_xy.length;
        int imgH = _inimg_xy[0].length;

        // this one makes sense for _iter>0
        // there is part that is executed if _iter==1, the rest is

        // take states at _iter-1
        float[][] prev_x = xt[_iter-1];
        float[]   prev_w = wt[_iter-1];

        int count = 0;
        float sum_lhoods = 0;// used to check validity
        float sum_prior = 0;

        for (int i = 0; i < prev_x.length; i++) { // will loop throughout the states Nstates

            float px        = prev_x[i][0];
            float py        = prev_x[i][1];
            float vx        = prev_x[i][2];
            float vy        = prev_x[i][3];
            float gcsstdc   = prev_x[i][4];

            if (i==0 & _iter==1) {

                predict(px, py, vx, vy, gcsstdc, getSemiCircleNrPoints(neuron_radius), _angula_deg, _gcsstd_pix);

                trans_xy   = new float[Nstates*pred_w.length][5];   // to store all the predictions at current iteration
                pties      = new float[Nstates*pred_w.length];      // store the probabilities temporarily
                lhoods     = new float[Nstates*pred_w.length];      // store the likelihoods

            }
            else
                predict(px, py, vx, vy, gcsstdc, _angula_deg, _gcsstd_pix);

            // the reminder is the same as in the other version of iter()
            for (int j = 0; j < pred_w.length; j++) { // loop all options, calculate likelihoods

                trans_xy[count][0] = p[j/gcsstd.length][0]; // px
                trans_xy[count][1] = p[j/gcsstd.length][1]; // py
                trans_xy[count][2] = v[j/gcsstd.length][0]; // vx
                trans_xy[count][3] = v[j/gcsstd.length][1]; // vy
                trans_xy[count][4] = gcsstd[j%gcsstd.length]; // gcsstd[]

                lhoods[count] = 0; // will take the highest calculated correlation out of all possible correlations
                // check if this sampling range falls into the image at all - borders are enough to check
                float p1x = trans_xy[count][0] + L2*cross_profile_step*trans_xy[count][3]; // x+L*max_step*vy
                float p1y = trans_xy[count][1] - L2*cross_profile_step*trans_xy[count][2]; // y-L*max_step*vx
                float p2x = trans_xy[count][0] - L2*cross_profile_step*trans_xy[count][3]; // x-L*max_step*vy
                float p2y = trans_xy[count][1] + L2*cross_profile_step*trans_xy[count][2]; // y+L*max_step*vx

                if ((p1x>=0 && p1x<=imgW-1 && p1y>=0 && p1y<=imgH-1) && (p2x>=0 && p2x<=imgW-1 && p2y>=0 && p2y<=imgH-1)) // both limits fall within the image
                    if (mode==0)
                        lhoods[count] = getVal(trans_xy[count][0], trans_xy[count][1], _inimg_xy);
                    else if (mode==1)
                        lhoods[count] = zncc(trans_xy[count][0], trans_xy[count][1], trans_xy[count][2], trans_xy[count][3], j % gcsstd.length, _inimg_xy);
                    else
                        System.exit(-1);


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
        int[] sort_idx = Tools.descending(pties);      // ptes will be sorted as a side effect

        int cnt = 0;
        for (int i = 0; i < sort_idx.length; i++) {

            xt[_iter][cnt][0] = trans_xy[sort_idx[i]][0];
            xt[_iter][cnt][1] = trans_xy[sort_idx[i]][1];
            xt[_iter][cnt][2] = trans_xy[sort_idx[i]][2];
            xt[_iter][cnt][3] = trans_xy[sort_idx[i]][3];
            xt[_iter][cnt][4] = trans_xy[sort_idx[i]][4];

            wt[_iter][cnt]    = pties[i]; // because they are already sorted
            cnt++;

            if (cnt==Nstates) break; // add top Ns points that are in the circle for further tracing

        }

        Stat.probability_distribution(wt[_iter]); // _wt will be normalized

        // set xc, rc at the end of the iteration
        xc[_iter][0] = 0;
        xc[_iter][1] = 0;
        rc[_iter]    = 0;

        for (int j = 0; j < xt[_iter].length; j++) {
            xc[_iter][0]    += wt[_iter][j] * xt[_iter][j][0];
            xc[_iter][1]    += wt[_iter][j] * xt[_iter][j][1];
            rc[_iter]       += wt[_iter][j] * xt[_iter][j][4];
        }

    }

    private void iter(
            float                   x0,
            float                   y0,
            float                   vx0,
            float                   vy0,
            float                   gcsstdc0,
            float[][]               _inimg_xy,
            float                   _angula_deg,
            float                   _gcsstd_pix
    )
    {

        // take states at _iter-1
//        float[][] prev_x = _xt[_iter-1];
//        float[]   prev_w = _wt[_iter-1];
        // normally there is a list of states, but now there is only one so this step will be skipped
//        for (int i = 0; i < prev_x.length; i++) {
//            float x  = prev_x[i][0]; float y  = prev_x[i][1];
//            float vx = prev_x[i][2]; float vy = prev_x[i][3];
//        }
        int L = tt[0].length; // length is the same for all sigmas
        int L2  = L/2;
        int imgW = _inimg_xy.length;
        int imgH = _inimg_xy[0].length;

        predict(x0, y0, vx0, vy0, gcsstdc0, Nstates, _angula_deg, _gcsstd_pix); // predict so that number of elements at the semicirc is set (will do allocation of p,v,pred_w inside)

        // some aux values for bayesian sequential filtering (Nsemicirc is set to Nstates in predict())
        trans_xy   = new float[1*pred_w.length][5];   // to store all the predictions at current iteration
        pties      = new float[1*pred_w.length];      // store the probabilities temporarily
        lhoods     = new float[1*pred_w.length];      // store the likelihoods

        int count = 0;
        float sum_lhoods = 0;// used to check validity
        float sum_prior = 0;

        for (int j = 0; j < pred_w.length; j++) { // loop all options, calculate likelihoods

            trans_xy[count][0] = p[j/gcsstd.length][0]; // px
            trans_xy[count][1] = p[j/gcsstd.length][1]; // py
            trans_xy[count][2] = v[j/gcsstd.length][0]; // vx
            trans_xy[count][3] = v[j/gcsstd.length][1]; // vy
            trans_xy[count][4] = gcsstd[j%gcsstd.length]; // gcsstd[]

            lhoods[count] = 0; // will take the highest calculated correlation out of all possible correlations
            // check if this sampling range falls into the image at all - borders are enough to check
            float p1x = trans_xy[count][0] + L2*cross_profile_step*trans_xy[count][3]; // x+L*max_step*vy
            float p1y = trans_xy[count][1] - L2*cross_profile_step*trans_xy[count][2]; // y-L*max_step*vx
            float p2x = trans_xy[count][0] - L2*cross_profile_step*trans_xy[count][3]; // x-L*max_step*vy
            float p2y = trans_xy[count][1] + L2*cross_profile_step*trans_xy[count][2]; // y+L*max_step*vx

            if ((p1x>=0 && p1x<=imgW-1 && p1y>=0 && p1y<=imgH-1) && (p2x>=0 && p2x<=imgW-1 && p2y>=0 && p2y<=imgH-1)) // both limits fall within the image
                if (mode==0)
                    lhoods[count] = getVal(trans_xy[count][0], trans_xy[count][1], _inimg_xy);
                else if (mode == 1)
                    lhoods[count] = zncc(trans_xy[count][0], trans_xy[count][1], trans_xy[count][2], trans_xy[count][3], j % gcsstd.length, _inimg_xy);
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
        int[] sort_idx = Tools.descending(pties);       // ptes will be sorted as a side effect

        int cnt = 0;
        for (int i = 0; i < sort_idx.length; i++) {

            xt[0][cnt][0] = trans_xy[sort_idx[i]][0];
            xt[0][cnt][1] = trans_xy[sort_idx[i]][1];
            xt[0][cnt][2] = trans_xy[sort_idx[i]][2];
            xt[0][cnt][3] = trans_xy[sort_idx[i]][3];
            xt[0][cnt][4] = trans_xy[sort_idx[i]][4];

            wt[0][cnt]    = pties[i];                   // because they are already sorted
            cnt++;

            if (cnt==Nstates) break;                    // add top Ns points that are in the circle for further tracing

        }

        Stat.probability_distribution(wt[0]);           // _wt will be normalized

        // set xc, rc at the end of the iteration
        xc[0][0] = 0;
        xc[0][1] = 0;
        rc[0]    = 0;

        for (int j = 0; j < xt[0].length; j++) {
            xc[0][0]    += wt[0][j] * xt[0][j][0];
            xc[0][1]    += wt[0][j] * xt[0][j][1];
            rc[0]       += wt[0][j] * xt[0][j][4];
        }

    }

    private float getVal(float px, float py, float[][] _inimg_xy)
    {
        return Interpolator.interpolateAt(px, py, _inimg_xy);
    }

    private float zncc(float px, float py, float vx, float vy, int gcsstd_idx, float[][] _inimg_xy)
    {

        int L2 = img_vals.length/2;

        img_vals[L2] = Interpolator.interpolateAt(px, py, _inimg_xy);
        float avgI = img_vals[L2];

        for (int shift = 1; shift <= L2; shift++) {

            img_vals[L2+shift] = Interpolator.interpolateAt(
                    px + shift * cross_profile_step * vy,
                    py - shift * cross_profile_step * vx,
                    _inimg_xy
            );

            avgI += img_vals[L2+shift];

            img_vals[L2-shift] = Interpolator.interpolateAt(
                    px - shift * cross_profile_step * vy,
                    py + shift * cross_profile_step * vx,
                    _inimg_xy
            );

            avgI += img_vals[L2-shift];

        }

        avgI = avgI / (float)(img_vals.length);

        // calculate zncc, use template that corresponds to this sigma
        float corra = 0;
        float corrb = 0;
        float corrc = 0;

        for (int profileidx = 0; profileidx <img_vals.length; profileidx++) {
            corra += (img_vals[profileidx]-avgI) * (tt[gcsstd_idx][profileidx]-tta[gcsstd_idx]);
            corrb += Math.pow(img_vals[profileidx]-avgI, 2);
            corrc += Math.pow(tt[gcsstd_idx][profileidx]-tta[gcsstd_idx], 2);
        }

        return  (float) (corra / Math.sqrt(corrb*corrc));

    }

    public Overlay getTrace()
    {

        float[] xline = new float[iter_counter];
        float[] yline = new float[iter_counter];

        Overlay ov = new Overlay();

        for (int i = 0; i < iter_counter; i++) {            // iterations

            float[] wgts = Stat.min_max_normalize1(wt[i]);  // visualizaton purposes

//            float min = Float.POSITIVE_INFINITY;
//            float max = Float.NEGATIVE_INFINITY;

            for (int j = 0; j <xt[i].length; j++) {         // loop all the states

//                xline[i] += wt[i][j] * xt[i][j][0];
//                yline[i] += wt[i][j] * xt[i][j][1];

                float r = wgts[j];
                OvalRoi p = new OvalRoi(xt[i][j][0]-r+.5, xt[i][j][1]-r+.5, 2*r, 2*r);
                p.setFillColor(new Color(1f, 1f, 1f, wgts[j]));
                ov.add(p);

//                if (wt[i][j]>max) max = wt[i][j];
//                if (wt[i][j]<min) min = wt[i][j];

            }

//            System.out.println("min = " + min + " max = " + max);

            xline[i] = xc[i][0]+.5f;
            yline[i] = xc[i][1]+.5f;

        }

        PolygonRoi pr = new PolygonRoi(xline, yline, PolygonRoi.FREELINE);
        pr.setFillColor(Color.YELLOW);
        pr.setStrokeColor(Color.YELLOW);
        pr.setStrokeWidth(0.5);
        ov.add(pr);

        return ov;
    }

}
