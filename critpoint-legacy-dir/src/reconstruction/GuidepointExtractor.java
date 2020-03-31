package reconstruction;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;

import java.awt.*;
import java.util.ArrayList;

/**
 *
 * will use Profiler2D/3D (included usage of Zncc2D/3D) and the ForegroundExtractor output
 * to get the guidepoints (points with high correlation)
 * output is the queue with points that have high correlation
 * that would be used to initialize tracing
 *
 * Created by miroslav on 5-3-15.
 */
public class GuidepointExtractor {

    public ArrayList<boolean[][]>   isqueue_loc;    // map of the guidepoints
    public ArrayList<int[]>         queue;          // list of selected points
    public ArrayList<float[]>       i2sigma;        // scale
    public ArrayList<float[][]>     i2dir;          // direction
    public ArrayList<float[]>       i2zncc;         // correlation

    public GuidepointExtractor(){
        isqueue_loc = new ArrayList<boolean[][]>();
        queue   = new ArrayList<int[]>();
        i2sigma = new ArrayList<float[]>();
        i2dir   = new ArrayList<float[][]>();
        i2zncc  = new ArrayList<float[]>();
    }

    public void work2D(
            float[][]           _inimg_xy,
            float               _neuron_radius,
            int[][]             _i2xy,
            ArrayList<double[]>  _soma_list,
            float               _correlation_boundary,
            int                 _mode,
            String              _midresults_dir
    ){

        Zncc2D zncc = new Zncc2D(_neuron_radius);

        Profiler2D.loadTemplate(
                zncc,
                _i2xy,
                _soma_list,
                _inimg_xy);

        int totalProfileComponents = _i2xy.length;

        int CPU_NR = Runtime.getRuntime().availableProcessors() + 1;

        Profiler2D pf_jobs[] = new Profiler2D[CPU_NR];

        for (int i = 0; i < pf_jobs.length; i++) {
            pf_jobs[i] = new Profiler2D(i*totalProfileComponents/CPU_NR,  (i+1)*totalProfileComponents/CPU_NR);
            pf_jobs[i].start();
        }
        for (int i = 0; i < pf_jobs.length; i++) {
            try {
                pf_jobs[i].join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }

        // add outputs of the Profiler2D to the class variables
        int[] obtained_queue = Profiler2D.createQueue(_correlation_boundary, _mode);
        queue.add(obtained_queue);
        isqueue_loc.add(Profiler2D.getQueueLocationMap(obtained_queue));
        i2sigma.add(Profiler2D.i2sigma.clone());
        i2dir.add(Profiler2D.i2vxy.clone());
        i2zncc.add(Profiler2D.i2zncc.clone());

        if (!_midresults_dir.equalsIgnoreCase("")) {
            IJ.saveAs(zncc.getSampling(),       "Tiff",     _midresults_dir + "sampling.tif");
            IJ.saveAs(zncc.getTemplates(),      "Tiff",     _midresults_dir + "templates.tif");
            ArrayList<int[][]> aa = new ArrayList<int[][]>(); aa.add(_i2xy.clone());
            IJ.saveAs(getCorrelations(aa),      "Tiff",     _midresults_dir + "zncc.tif");
            IJ.saveAs(getQueueMap(),            "Tiff",     _midresults_dir + "queueMap.tif");
        }

        Profiler2D.clean();

    }

    public void work3D(
            ArrayList<float[][]>    _inimg_xyz,
            float                   _neuron_radius,
            ArrayList<int[][]>      _i2xy,
            ArrayList<double[]>     _soma_list,
            float                   _correlation_boundary,
            int                     _mode,
            String                  _midresults_dir
    ){

        for (int layIdx = 0; layIdx < _inimg_xyz.size(); layIdx++)
            work2D(_inimg_xyz.get(layIdx), _neuron_radius, _i2xy.get(layIdx), _soma_list, _correlation_boundary, _mode, "");

        if (!_midresults_dir.equalsIgnoreCase("")) {
            Zncc2D zncc = new Zncc2D(_neuron_radius);
            IJ.saveAs(zncc.getSampling(),       "Tiff",     _midresults_dir + "sampling.tif");
            IJ.saveAs(zncc.getTemplates(),      "Tiff",     _midresults_dir + "templates.tif");
            IJ.saveAs(getCorrelations(_i2xy),   "Tiff",     _midresults_dir + "zncc.tif");
            IJ.saveAs(getQueueMap(),            "Tiff",     _midresults_dir + "queueMap.tif");
        }

    }


//    public void work3D(
//            float[][][]         _inimg_xyz,
//            float               _neuron_radius,
////            int                 _Ndirections,
//            float               _zDist,
//            int[][]             _i2xyz,
//            ArrayList<double[]>  _soma_list,
//            float               _correlation_boundary,
//            String              _midresults_dir
//    )
//    {
//        Zncc3D zncc = new Zncc3D(_neuron_radius, _zDist);
//
//
////        zncc.exportFilterSampling(_midresults_dir);
////        zncc.exportTemplates(_midresults_dir);
//
//        Profiler3D.loadTemplate(
//            zncc,
//            _i2xyz,
//            _soma_list,
//            _inimg_xyz
//        );
//
//        int totalProfileComponents = _i2xyz.length;
//
//        int CPU_NR = Runtime.getRuntime().availableProcessors() + 1;
//
//        Profiler3D pf_jobs[] = new Profiler3D[CPU_NR];
//
//        System.out.println("profiler 3d...");
//
//        for (int i = 0; i < pf_jobs.length; i++) {
//            pf_jobs[i] = new Profiler3D(i*totalProfileComponents/CPU_NR,  (i+1)*totalProfileComponents/CPU_NR);
//            pf_jobs[i].start();
//        }
//        for (int i = 0; i < pf_jobs.length; i++) {
//            try {
//                pf_jobs[i].join();
//            } catch (InterruptedException e) {
//                e.printStackTrace();
//            }
//        }
//
//        System.out.println("profiler 3d finished!");
//
//        // add outputs of the Profiler2D to the class variables
//        queue = Profiler3D.createQueue(_correlation_boundary, 0); // expansion type is there
//        isqueue_loc = Profiler3D.getQueueLocationMap(queue);
//
//        i2sigma = Profiler3D.i2sigma.clone();
//        i2dir = Profiler3D.i2vxyz.clone();
//
////        System.out.println("saving in " + _midresults_dir);
//        if (!_midresults_dir.equalsIgnoreCase("")) {
//            zncc.exportFilterSampling(_midresults_dir);
//            zncc.exportTemplates(_midresults_dir);
//            IJ.saveAs(Profiler3D.getScores(),   "Tiff", _midresults_dir + "zncc.tif");
//            Profiler3D.getQueueLocationSwc(queue, _midresults_dir + "queue.swc");
//        }
//
//        Profiler3D.clean();
//
//
//    }

    public ImagePlus getQueueMap() {

        int W = isqueue_loc.get(0).length;
        int H = isqueue_loc.get(0)[0].length;
        ImageStack isqueue = new ImageStack(W, H);

        Overlay ov = new Overlay();

        for (int i = 0; i < isqueue_loc.size(); i++) {

            byte[] arry = new byte[W*H];

            for (int x = 0; x < W; x++) {
                for (int y = 0; y < H; y++) {
                    arry[y*W+x] = (isqueue_loc.get(i)[x][y])?(byte)255:(byte)0;
                    if (isqueue_loc.get(i)[x][y]) {
                        OvalRoi ovr = new OvalRoi(x+.5-1, y+.5-1, 2, 2);
                        ovr.setStrokeColor(Color.RED);
                        ovr.setStrokeWidth(1.5);
//                        ovr.setFillColor(new Color(1,0,0,0.4f));
                        ovr.setPosition(i+1);
                        ov.add(ovr);
                    }
                }
            }

            isqueue.addSlice(new ByteProcessor(W, H, arry));
        }

        ImagePlus iout = new ImagePlus("queue_map", isqueue);
        iout.setOverlay(ov);
        iout.show();
        return iout;
    }

    public ImagePlus getCorrelations(ArrayList<int[][]> i2xy){

        int W = isqueue_loc.get(0).length;
        int H = isqueue_loc.get(0)[0].length;

        ImageStack ipzncc = new ImageStack(W, H);

        for (int i = 0; i < i2xy.size(); i++) {

            float[] arry = new float[W*H];

            for (int j = 0; j < i2xy.get(i).length; j++) {
                int x = i2xy.get(i)[j][0];
                int y = i2xy.get(i)[j][1];
                arry[y*W+x] = i2zncc.get(i)[j];
            }

            ipzncc.addSlice(new FloatProcessor(W, H, arry));

        }

        return new ImagePlus("zncc_map", ipzncc);

    }

}
