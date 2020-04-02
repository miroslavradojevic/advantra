package tracing2d;

import aux.Stat;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.*;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.ArrayList;

/**
 * Created by miroslav on 30-10-14.
 */
public class BayesianTracer2D_Demo implements PlugIn, MouseListener {

    float           _D                   ;    // will define the step in Bayesian tracking scheme - size of the sampling sphere
    String          _scales_list         ;    // scales used when calculating neuriteness
    float           _prior_sigma_deg     ;    // gaussian sigma for the angular prediction
    int             _Nt                  ;    // number of samples to maintain throughout recursive estimation
    int             _MAX_ITER            ;    // limit number of iterations when tracing recursively
    int             _MIN_ITER            ;    // minimum amount of iterations in recursive tracing to elapse before considering joining some other region
    boolean         _USE_ORIGINAL        ;    // use original as input measurement instead of neuriteness

    BayesianTracer2D    btracer;
    ImagePlus curr_img;
    String _image_path;

    Overlay         curr_ov;
    int             count_click;
    int             curr_center_x;
    int             curr_center_y;
    ImageCanvas curr_can;

    ArrayList<float[][][]>  tracer2d_Xt_xy      = new ArrayList<float[][][]>();  // (trace#)[trace_len][elements][xy]
    ArrayList<float[][]>    tracer2d_wt_xy      = new ArrayList<float[][]>();    // (trace#)[trace_len][elements]
    ArrayList<float[][]>    tracer2d_est_xy     = new ArrayList<float[][]>();    // (trace#)[trace_len][xy]

    public void mouseClicked(MouseEvent e)
    {
        int x = e.getX();
        int y = e.getY();
        int offscreenX = curr_can.offScreenX(x);
        int offscreenY = curr_can.offScreenY(y);
        count_click++;



        if (count_click==1) {
            curr_center_x = offscreenX; // only one starting state (expect more in sequential tracing)
            curr_center_y = offscreenY;

            // add center marker
            float r = 1.0f;
            OvalRoi central = new OvalRoi(offscreenX-r+.5, offscreenY-r+.5f, 2*r, 2*r);
            central.setStrokeColor(Color.RED);
            central.setFillColor(Color.RED);
            curr_ov.add(central);

        }
        else if (count_click==2) { /// bayesian tracking iteration


            float[] vxy = new float[]{offscreenX-curr_center_x, offscreenY-curr_center_y};

            int out_label = btracer.run(        // bayesian recursive tracing
                    // to start the trace
                    curr_center_x,
                    curr_center_y,
                    vxy[0],
                    vxy[1],

                    1,
                    new int[curr_img.getWidth()][curr_img.getHeight()]     // to stop the trace

            );

            boolean draw_states = true;

            /*
                draw direction line
             */
            Line ln = new Line(offscreenX+.5f, offscreenY+.5f, curr_center_x+.5f, curr_center_y+.5f);
            ln.setStrokeColor(Color.RED);
            ln.setFillColor(Color.RED);
            if (draw_states) ln.setStrokeWidth(0.3f);
            else ln.setStrokeWidth(1f);
            curr_ov.add(ln);

            float ll = (float) Math.sqrt(Math.pow(curr_center_y-offscreenY,2)+Math.pow(curr_center_x-offscreenX,2));
            float dx = -curr_center_x+offscreenX;
            dx /=ll;
            float dy = -curr_center_y+offscreenY;
            dy/=ll;

            float LL = 0.85f;
            float KK = 0.05f;
            ln = new Line(curr_center_x+.5f+ LL * ll * dx - dy * KK * ll, curr_center_y+.5f+ LL * ll * dy + dx * KK * ll,  offscreenX+.5f , offscreenY+.5f );
            ln.setStrokeColor(Color.RED);
            ln.setFillColor(Color.RED);
            if (draw_states) ln.setStrokeWidth(0.3f);
            else ln.setStrokeWidth(1f);
            curr_ov.add(ln);

            ln = new Line(curr_center_x+.5f+ LL * ll * dx + dy * KK * ll, curr_center_y+.5f+ LL * ll * dy - dx * KK * ll,  offscreenX+.5f , offscreenY+.5f );
            ln.setStrokeColor(Color.RED);
            ln.setFillColor(Color.RED);
            if (draw_states) ln.setStrokeWidth(0.3f);
            else ln.setStrokeWidth(1f);
            curr_ov.add(ln);


            /*
            add the traced variables
             */

            // *** btracer.Xt_xy -> tracer2d_Xt_xy
            float[][][] _Xt_xy = new float[btracer.Xt_xy.size()][][];

            for (int i = 0; i < btracer.Xt_xy.size(); i++) {

                _Xt_xy[i] = new float[btracer.Xt_xy.get(i).length][2];

                for (int j = 0; j < btracer.Xt_xy.get(i).length; j++) {
                    _Xt_xy[i][j][0] = btracer.Xt_xy.get(i)[j][0];
                    _Xt_xy[i][j][1] = btracer.Xt_xy.get(i)[j][1];
                }

            }

            tracer2d_Xt_xy.add(_Xt_xy);


            // *** btracer.w_xy -> tracer2d_wt_xy
            float[][] _wt_xy = new float[btracer.wt_xy.size()][];

            for (int i = 0; i < btracer.wt_xy.size(); i++) {

                _wt_xy[i] = new float[btracer.wt_xy.get(i).length];

                for (int j = 0; j < btracer.wt_xy.get(i).length; j++) {
                    _wt_xy[i][j] = btracer.wt_xy.get(i)[j];
                }

            }

            tracer2d_wt_xy.add(_wt_xy);

            // *** btracer.est_xy -> tracer2d_est_xy
            float[][] _est_xy = new float[btracer.est_xy.size()][];

            for (int i = 0; i < btracer.est_xy.size(); i++) {

                _est_xy[i] = new float[btracer.est_xy.get(i).length];

                for (int j = 0; j < btracer.est_xy.get(i).length; j++) {
                    _est_xy[i][j] = btracer.est_xy.get(i)[j];
                }

            }

            tracer2d_est_xy.add(_est_xy);



            /*
                draw estimated trace
             */

            Overlay ovl = new Overlay();


            if (draw_states) {
                Overlay trace_Xt_xy = getTrace_Xt_xy();
                for (int i = 0; i < trace_Xt_xy.size(); i++) curr_ov.add(trace_Xt_xy.get(i));
            }

            Overlay trace_est_xy = getTrace_est_xy();
            for (int i = 0; i < trace_est_xy.size(); i++) curr_ov.add(trace_est_xy.get(i));

//          setOverlay(ovl);
//          orig.show();
//          curr_can.setOverlay(ovl);


//            float[] est_x = new float[est_xy.size()];
//            float[] est_y = new float[est_xy.size()];
//
//            for (int i = 0; i < est_xy.size(); i++) {
//                est_x[i] = est_xy.get(i)[0];
//                est_y[i] = est_xy.get(i)[1];
//            }
//            PolygonRoi pr = new PolygonRoi(est_x, est_y, PolygonRoi.FREELINE);
//            pr.setFillColor(Color.YELLOW);
//            pr.setStrokeWidth(2);
//            curr_ov.add(pr);

        }
        else if (count_click==3) { // reset tracing

            count_click=0;
            curr_ov.clear();

            // remove the traces that were there earlier
            tracer2d_Xt_xy.clear();
            tracer2d_est_xy.clear();
            tracer2d_wt_xy.clear();

        }

        curr_img.draw();
    }

    public void mousePressed(MouseEvent e) {}
    public void mouseReleased(MouseEvent e) {}
    public void mouseEntered(MouseEvent e) {}
    public void mouseExited(MouseEvent e) {}

    public void run(String s) {

        // load the image
        String in_folder = Prefs.get("id.folder", System.getProperty("user.home"));
        OpenDialog.setDefaultDirectory(in_folder);
        OpenDialog dc = new OpenDialog("Select image file");
        in_folder = dc.getDirectory();
        _image_path = dc.getPath();
        if (_image_path==null) return;
        Prefs.set("id.folder", in_folder);

        curr_img = new ImagePlus(_image_path);


        _D              =   (float) Prefs.get("critpoint.tracing2d.d", 5f);
        _scales_list    =           Prefs.get("critpoint.tracing2d.scales_list", "3,5,7,9");
        _prior_sigma_deg=   (float) Prefs.get("critpoint.tracing2d.prior_sigma_deg", 65f);
        _Nt             =   (int)   Prefs.get("critpoint.tracing2d.Nt", 50);
        _MAX_ITER       =   (int)   Prefs.get("critpoint.tracing2d.max_iter", 200);
        _MIN_ITER       =   (int)   Prefs.get("critpoint.tracing2d.min_iter", 5);
        _USE_ORIGINAL   =           Prefs.get("critpoint.tracing2d.use_original", true);

        GenericDialog gd = new GenericDialog("BayesianTracer2D_Demo");

        gd.addNumericField("D", _D, 2, 10, "[pix.]");
        gd.addStringField("ScalesList", _scales_list, 30);
        gd.addNumericField("AngularPriorSigma", _prior_sigma_deg, 2, 10, "[deg.]");
        gd.addNumericField("Nt",                       _Nt,                 0,  10, "count");
        gd.addNumericField("MAX_ITER",                 _MAX_ITER,           0,  15, "");
        gd.addNumericField("MIN_ITER",                 _MIN_ITER,           0,  15, "");
        gd.addCheckbox("USE_ORIGINAL", _USE_ORIGINAL);

        gd.showDialog();
        if (gd.wasCanceled()) return;

        _D                  = (float) gd.getNextNumber();       Prefs.set("critpoint.tracing2d.d",               _D);
        _scales_list        =         gd.getNextString();       Prefs.set("critpoint.tracing2d.scales_list",     _scales_list);
        _prior_sigma_deg    = (float) gd.getNextNumber();       Prefs.set("critpoint.tracing2d.prior_sigma_deg", _prior_sigma_deg);
        _Nt                 = (int)   gd.getNextNumber();       Prefs.set("critpoint.tracing2d.Nt",              _Nt);
        _MAX_ITER           = (int)   gd.getNextNumber();       Prefs.set("critpoint.tracing2d.max_iter",        _MAX_ITER);
        _MIN_ITER           = (int)   gd.getNextNumber();       Prefs.set("critpoint.tracing2d.min_iter",        _MIN_ITER);
        _USE_ORIGINAL       =         gd.getNextBoolean();      Prefs.set("critpoint.tracing2d.use_original",    _USE_ORIGINAL);

        btracer = new BayesianTracer2D(
                _scales_list,
                curr_img,
                _D,
                _prior_sigma_deg,
                _Nt,
                _MAX_ITER,
                _MIN_ITER,
                _USE_ORIGINAL
        );

        System.out.println("ready!");
        count_click = 0; // reset the counts

        curr_center_x = -1;//Float.NaN;
        curr_center_y = -1;//Float.NaN;

        curr_ov = new Overlay();
        curr_img.show();
        curr_can = curr_img.getCanvas();
        curr_can.addMouseListener(this);
        curr_img.setOverlay(curr_ov);

        IJ.setTool("hand");

        System.out.println("Click!");

    }

    private Overlay getTrace_Xt_xy(){

        Overlay ov = new Overlay();

        // read tracer2d_Xt_xy and output the states

        for (int i = 0; i < tracer2d_Xt_xy.size(); i++) { // i - trace


                for (int j = 0; j < tracer2d_Xt_xy.get(i).length; j++) {  // j - iteration along the trace


                    // weights at this iteration will define the radiuses
                    float[] radiuses = tracer2d_wt_xy.get(i)[j].clone();
                    Stat.normalize(radiuses); // max. radius becomes 1

                    for (int k = 0; k < tracer2d_Xt_xy.get(i)[j].length; k++) { // k - one particle state

                        Color colr = new Color(1f,1f,1f,radiuses[k]);
                        OvalRoi overlay_Xt = new OvalRoi(
                                tracer2d_Xt_xy.get(i)[j][k][0]+.5f - 0.5f*radiuses[k],
                                tracer2d_Xt_xy.get(i)[j][k][1]+.5f - 0.5f*radiuses[k],
                                radiuses[k], radiuses[k]);
                        overlay_Xt.setStrokeColor(colr);
                        overlay_Xt.setFillColor(colr);

                        ov.add(overlay_Xt);

                    }

                }


        }

        return ov;

    }

    private Overlay getTrace_est_xy(){

        Overlay ov = new Overlay();

        for (int i = 0; i < tracer2d_est_xy.size(); i++) { // i - trace

            if (tracer2d_est_xy.get(i)!=null) { // if it is null then -> nothing



                float[] trc_x = new float[tracer2d_est_xy.get(i).length];
                float[] trc_y = new float[tracer2d_est_xy.get(i).length];

                for (int j = 0; j < tracer2d_est_xy.get(i).length; j++) { // j - iteration along the trace

                    trc_x[j] = tracer2d_est_xy.get(i)[j][0]   + .5f;
                    trc_y[j] = tracer2d_est_xy.get(i)[j][1]   + .5f;

                }

                PolygonRoi trc = new PolygonRoi(trc_x, trc_y, PolygonRoi.FREELINE);
                trc.setStrokeWidth(0.4f);
                trc.setStrokeColor(Color.YELLOW);
                trc.setFillColor(Color.YELLOW);
                ov.add(trc);

            }

        }


        return ov;

    }

}
