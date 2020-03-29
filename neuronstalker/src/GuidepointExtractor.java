import cuda.Profiler2Dcuda;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Line;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;

import java.awt.*;
import java.io.*;
import java.util.ArrayList;

/**
 *
 * will use Profiler2D (included usage of Zncc2D) and the ForegroundExtractor output
 * to get the queue of guidepoints (points with high correlation, above given threshold)
 * output is the queue with foreground location indexes that have high correlation
 * that would be used to initialize tracing
 *
 * Created by miroslav on 5-3-15.
 */
public class GuidepointExtractor {

    public ArrayList<int[][]>       gp2lab;                 // label map of the guidepoint locations
    public ArrayList<int[]>         gp2i;                   // list of guidepoint foreground indexes
    public ArrayList<float[]>       i2sigma;                // scale for foreground locations
    public ArrayList<float[][]>     i2dir;                  // direction for foreground locations
    public ArrayList<float[]>       i2zncc;                 // correlation for foreground locations

    public GuidepointExtractor(){
        gp2lab              = new ArrayList<int[][]>();
        gp2i                = new ArrayList<int[]>();
        i2sigma             = new ArrayList<float[]>();
        i2dir               = new ArrayList<float[][]>();
        i2zncc              = new ArrayList<float[]>();
    }

    public void work2D(
            float[][]           _inimg_xy,
            float               _neuron_radius,
            int[][]             _i2xy,
            ArrayList<double[]> _soma_list,
            float               _correlation_boundary,
            int                 _mode,
            String              _midresults_dir,
            boolean             _doCuda,
            int                 _block_size,
            int                 _max_threads_nr
    )
    {

        work(_inimg_xy, _neuron_radius, _i2xy, _soma_list, _doCuda, _block_size, _max_threads_nr);

        int imgW = _inimg_xy.length;
        int imgH = _inimg_xy[0].length;

        gp2i = new ArrayList<int[]>(1);
        gp2lab = new ArrayList<int[][]>(1);

        ArrayList<int[][]> i2xy = new ArrayList<int[][]>(1);
        i2xy.add(_i2xy.clone());

        getGuidepoint(
                _correlation_boundary,
                _mode,
                _soma_list,
                i2zncc, // it is a list after work()
                i2sigma,
                i2xy,
                imgW,
                imgH,
                Integer.MIN_VALUE,
                Float.MIN_VALUE,
                gp2i,    // out_gp2i = new ArrayList<int[]>(nr_layers);
                gp2lab   // out_gp2lab = new ArrayList<int[][]>(nr_layers);
        );

        if (!_midresults_dir.equalsIgnoreCase("")) {

            Profiler2Dcuda p2Dcuda = new Profiler2Dcuda(_neuron_radius, _i2xy, _inimg_xy, _block_size, _max_threads_nr);
            IJ.saveAs(p2Dcuda.getKernels(),            "Tiff",     _midresults_dir + "kernels_cuda.tif");

            Zncc2D zncc = new Zncc2D(_neuron_radius);
            IJ.saveAs(zncc.getKernels(), "Tiff",        _midresults_dir + "kernels.tif");
            IJ.saveAs(zncc.getSampling(), "Tiff",       _midresults_dir + "sampling.tif");
            IJ.saveAs(zncc.getTemplates(), "Tiff", _midresults_dir + "templates.tif");

            IJ.saveAs(getCorrelations(i2xy, imgW, imgH), "Tiff",    _midresults_dir + "zncc.tif");
            IJ.saveAs(getQueueMap(i2xy), "Tiff",                    _midresults_dir + "queue.tif");
            queueSwc(i2xy,                                          _midresults_dir + "queue.swc");
        }

    }

    public void work3D(
            ArrayList<float[][]>    _inimg_xyz,
            float                   _zDist,
            float                   _neuron_radius,
            ArrayList<int[][]>      _i2xy,
            ArrayList<double[]>     _soma_list,
            float                   _correlation_boundary,
            int                     _mode,
            String                  _midresults_dir,
            boolean                 _doCuda,
            int                     _block_size,
            int                     _max_threads_nr
    )
    {

        for (int layIdx = 0; layIdx < _inimg_xyz.size(); layIdx++)
            work(_inimg_xyz.get(layIdx), _neuron_radius, _i2xy.get(layIdx), _soma_list, _doCuda, _block_size, _max_threads_nr);

        int imgW = _inimg_xyz.get(0).length;
        int imgH = _inimg_xyz.get(0)[0].length;
        int imgL = _inimg_xyz.size();

        int nr_layers = _i2xy.size();

        gp2i = new ArrayList<int[]>(nr_layers);
        gp2lab = new ArrayList<int[][]>(nr_layers);

        getGuidepoint(
                _correlation_boundary,
                _mode,
                _soma_list,
                i2zncc,
                i2sigma,
                _i2xy,
                imgW,
                imgH,
                imgL,
                _zDist,
                gp2i,
                gp2lab
        );

        if (!_midresults_dir.equalsIgnoreCase("")) {

            Profiler2Dcuda p2Dcuda = new Profiler2Dcuda(_neuron_radius, _i2xy.get(0), _inimg_xyz.get(0), _block_size, _max_threads_nr);
            IJ.saveAs(p2Dcuda.getKernels(),            "Tiff",     _midresults_dir + "kernels_cuda.tif");

            Zncc2D zncc = new Zncc2D(_neuron_radius);
            IJ.saveAs(zncc.getKernels(), "Tiff",        _midresults_dir + "kernels.tif");
            IJ.saveAs(zncc.getSampling(), "Tiff",       _midresults_dir + "sampling.tif");
            IJ.saveAs(zncc.getTemplates(), "Tiff", _midresults_dir + "templates.tif");

            IJ.saveAs(getCorrelations(_i2xy, imgW, imgH),   "Tiff",     _midresults_dir + "zncc.tif");
            IJ.saveAs(getQueueMap(_i2xy),                   "Tiff",     _midresults_dir + "queue.tif");
            queueSwc(_i2xy,                                             _midresults_dir + "queue.swc");
        }

    }

    private void work(
            float[][]           _inimg_xy,
            float               _neuron_radius,
            int[][]             _i2xy,
            ArrayList<double[]> _soma_list,
            boolean             _doCuda,
            int                 _block_size,
            int                 _max_threads_nr
    )
    {

        if (_doCuda) {

            Profiler2Dcuda p2Dcuda = new Profiler2Dcuda(_neuron_radius, _i2xy, _inimg_xy, _block_size, _max_threads_nr);

            p2Dcuda.run();

            i2zncc.add(p2Dcuda.i2zncc_cuda.clone());

            for (int i = 0; i < p2Dcuda.i2sigma_cuda.length; i++)
                p2Dcuda.i2sigma_cuda[i] = p2Dcuda.sigmas[(int) Math.round(p2Dcuda.i2sigma_cuda[i])];
            i2sigma.add(p2Dcuda.i2sigma_cuda.clone());


            float[][] i2dir_cuda = new float[p2Dcuda.i2vx_cuda.length][2];
            for (int i = 0; i < p2Dcuda.i2vx_cuda.length; i++) {
                i2dir_cuda[i][0] = p2Dcuda.i2vx_cuda[i];
                i2dir_cuda[i][1] = p2Dcuda.i2vy_cuda[i];
            }
            i2dir.add(i2dir_cuda);

            p2Dcuda.clean();

        }
        else {

            Zncc2D zncc = new Zncc2D(_neuron_radius);

            Profiler2D.loadTemplate(zncc, _i2xy, _soma_list, _inimg_xy);

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

            i2zncc.add(Profiler2D.i2zncc.clone());
            i2sigma.add(Profiler2D.i2sigma.clone());
            i2dir.add(Profiler2D.i2vxy.clone());

            Profiler2D.clean();

        }
    }

    private void getGuidepoint(
            float                   boundary,
            int                     weighting_type,
            ArrayList<double[]>     soma_list,
            ArrayList<float[]>      _i2zncc,
            ArrayList<float[]>      _i2sigma,
            ArrayList<int[][]>      _i2xy,
            int _W,
            int _H,
            int _L,
            float _zDist,
            ArrayList<int[]>        out_gp2i,    // out_gp2i = new ArrayList<int[]>(nr_layers);
            ArrayList<int[][]>      out_gp2lab   // out_gp2lab = new ArrayList<int[][]>(nr_layers);
    )
    {

        ArrayList<Float>    sel_scrs = new ArrayList<Float>();
        ArrayList<Integer>  sel_idxs = new ArrayList<Integer>();
        ArrayList<Integer>  sel_lays = new ArrayList<Integer>();

        int nr_layers = _i2zncc.size();

        for (int i = 0; i < nr_layers; i++) {
            for (int j = 0; j < _i2zncc.get(i).length; j++) {
                if (_i2zncc.get(i)[j]>boundary) {
                    sel_idxs.add(j);
                    sel_scrs.add(_i2zncc.get(i)[j]);
                    sel_lays.add(i);
                }
            }
        }

        if (weighting_type==2){ // don't consider distances to soma
        }
        else if (weighting_type==0 || weighting_type==1) {
            if (soma_list.size() >= 1) {
                // there is 1+ soma's

                float min_dist = Float.POSITIVE_INFINITY;
                float max_dist = Float.NEGATIVE_INFINITY;

                // add the distance to the scores of selected locations
                ArrayList<Float>    sel_scrs_dist_to_soma = new ArrayList<Float>(sel_idxs.size());

                for (int i = 0; i < sel_idxs.size(); i++) {

                    float currx = _i2xy.get(sel_lays.get(i))[sel_idxs.get(i)][0];
                    float curry = _i2xy.get(sel_lays.get(i))[sel_idxs.get(i)][1];
                    float currz = sel_lays.get(i);

                    float min_dist_soma = Float.POSITIVE_INFINITY;
                    for (int j = 0; j < soma_list.size(); j++) {
                        float dist1 = (float) (
                                Math.pow(currx - soma_list.get(j)[0], 2) +
                                        Math.pow(curry - soma_list.get(j)[1], 2) +
                                            Math.pow((currz-soma_list.get(j)[2])*_zDist, 2)
                        );
                        if (dist1 < min_dist_soma) min_dist_soma = dist1;
                    }

                    sel_scrs_dist_to_soma.add(i, min_dist_soma);

                    if (min_dist_soma < min_dist) min_dist = min_dist_soma;
                    if (min_dist_soma > max_dist) max_dist = min_dist_soma;

                }

                // loop through all the points once they're calculated append the zncc score to the min/max normalized distances to soma by now
                for (int i = 0; i < sel_scrs.size(); i++) {

                    float normalized_dist = -1;

                    if (weighting_type == 0) { // those closer to the soma are weighted higher
                        normalized_dist = 1 - (sel_scrs_dist_to_soma.get(i) - min_dist) / (max_dist - min_dist);
                    } else if (weighting_type == 1) { // those that are the most distant from soma are weighted higher
                        normalized_dist = (sel_scrs_dist_to_soma.get(i) - min_dist) / (max_dist - min_dist);
                    }

                    sel_scrs.set(i, normalized_dist * sel_scrs.get(i));

                }

            }
            else {
                // there is no soma, make list without taking into account the distances towards soma
//                for (int i = 0; i < sel_idxs.size(); i++) sel_scrs.add(i, _i2zncc[sel_idxs.get(i)]);
            }
        }


        int[] queue = Toolbox.descendingFloat(sel_scrs);
//        for (int i = 0; i < queue.length; i++)
//            queue[i] = sel_idxs.get(queue[i]);
        // now select subset of the queue elements based on supressing the neighbourhood of the selected points
//        boolean[][] map = new boolean[_W][_H];
//        int[][] gp_map_xy = new int[_W][_H]; // allocate output, initialize with zeros

        for (int i = 0; i < nr_layers; i++) out_gp2lab.add(i, new int[_W][_H]);

        ArrayList<Integer> reduced_lays = new ArrayList<Integer>(); // will be combination (layer-indexFromThatLayer)
        ArrayList<Integer> reduced_idxs = new ArrayList<Integer>();

        int count_them = 0;

        for (int i = 0; i < queue.length; i++) { // start from  those that were largest weight

            int ll      = sel_lays.get(queue[i]);
            int ii      = sel_idxs.get(queue[i]);

            int xloc = _i2xy.get(ll)[ii][0]; // take loc index within that layer
            int yloc = _i2xy.get(ll)[ii][1];
            int zloc = ll;
            int rloc = (int) Math.ceil(_i2sigma.get(ll)[ii]);

            if (nr_layers>1) {
                if (xloc<rloc || xloc>=_W-rloc || yloc<rloc || yloc>=_H-rloc || Math.ceil(zloc*_zDist)<rloc || Math.ceil((_L - zloc)*_zDist)<=rloc)
                    continue;
            }

            if (nr_layers==1) {
                if (xloc<rloc || xloc>=_W-rloc || yloc<rloc || yloc>=_H-rloc)
                    continue;
            }

            boolean overlaps = false;

            for (int j = -rloc; j <= rloc; j++) {
                for (int k = -rloc; k <= rloc; k++) {

                    if (nr_layers==1) {
                        if ((int)Math.round(Math.sqrt(j*j+k*k))<=rloc) {
                            if(out_gp2lab.get(zloc)[xloc+j][yloc+k]>0)
                                overlaps = true;
                        }
                    }
                    else{
                        for (int l = -rloc; l <= rloc; l++) {
                            if ((int)Math.round(Math.sqrt(j*j+k*k+(l*_zDist)*(l*_zDist)))<=rloc) {
                                if(out_gp2lab.get(zloc+l)[xloc+j][yloc+k]>0)
                                    overlaps = true;
                            }
                        }
                    }

                }
            }

            if (!overlaps) {

                reduced_lays.add(ll); // add it to the reduced lists
                reduced_idxs.add(ii);

                count_them++;

                // update the map
                for (int j = -rloc; j <= rloc; j++) {
                    for (int k = -rloc; k <= rloc; k++) {
                        if (nr_layers==1) {
                            if ((int)Math.round(Math.sqrt(j*j+k*k))<=rloc)
                                out_gp2lab.get(zloc)[xloc + j][yloc + k] = count_them;
                        }
                        else {
                            for (int l = -rloc; l <= rloc; l++) {
                                if ((int)Math.round(Math.sqrt(j*j+k*k+(l*_zDist)*(l*_zDist)))<=rloc)
                                    out_gp2lab.get(zloc+l)[xloc + j][yloc + k] = count_them;
                            }
                        }
                    }
                }
            }

        }

        ArrayList<ArrayList<Integer>> out_gp2i_temp = new ArrayList<ArrayList<Integer>>(nr_layers);
        for (int i = 0; i < nr_layers; i++) out_gp2i_temp.add(new ArrayList<Integer>());

        for (int i = 0; i < reduced_lays.size(); i++) {
            int at_lay = reduced_lays.get(i);
            int at_idx = reduced_idxs.get(i);
            out_gp2i_temp.get(at_lay).add(at_idx);
        }

        // now it is necesary to allocate int[] at each layer
        for (int i = 0; i < nr_layers; i++) {
            if (out_gp2i_temp.get(i).size()>0) {
                int[] idxs_to_add = new int[out_gp2i_temp.get(i).size()];
                for (int j = 0; j < out_gp2i_temp.get(i).size(); j++) idxs_to_add[j] = out_gp2i_temp.get(i).get(j);
                out_gp2i.add(i, idxs_to_add);
            }
            else out_gp2i.add(i ,null);
        }

    }

    private void queueSwc(ArrayList<int[][]> _i2xy, String swc_path) {

        PrintWriter logWriter = null;
        try {
            logWriter = new PrintWriter(swc_path);
            logWriter.print("");
            logWriter.close();
        } catch (FileNotFoundException ex) {}

        try {
            logWriter = new PrintWriter(new BufferedWriter(new FileWriter(swc_path, true)));
            logWriter.println("#");
        } catch (IOException e) {}

        
        int cnt = 1;
        for (int i = 0; i < gp2i.size(); i++) {
            if (gp2i.get(i)!=null) {
                for (int j = 0; j < gp2i.get(i).length; j++) {

                    int guidepoint_foreground_index = gp2i.get(i)[j];

                    float x = _i2xy.get(i)[guidepoint_foreground_index][0];
                    float y = _i2xy.get(i)[guidepoint_foreground_index][1];
                    float z = i;

                    float r = i2sigma.get(i)[guidepoint_foreground_index];

                    logWriter.println(String.format(
                            "%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d",
                            cnt, 6, x, y, z, r, -1));
                    cnt++;

                }
            }
        }


        logWriter.close();


    }

    private ImagePlus getQueueMap(ArrayList<int[][]> _i2xy) {

        int W = gp2lab.get(0).length;
        int H = gp2lab.get(0)[0].length;
        ImageStack isqueue = new ImageStack(W, H);

        for (int i = 0; i < gp2lab.size(); i++) {

            int[] arry = new int[W*H];

            for (int x = 0; x < W; x++) {
                for (int y = 0; y < H; y++) {
                    arry[y*W+x] = gp2lab.get(i)[x][y];    //?(byte)255:(byte)0;
                }
            }

            isqueue.addSlice(new FloatProcessor(W, H, arry));
        }

        Overlay ov = new Overlay();
        ov.drawNames(true);

        for (int j = 0; j < gp2i.size(); j++) {

            if (gp2i.get(j)!=null) {

                for (int ii = 0; ii < gp2i.get(j).length; ii++) {

                    int xx = _i2xy.get(j)[gp2i.get(j)[ii]][0];
                    int yy = _i2xy.get(j)[gp2i.get(j)[ii]][1];
                    int rr = (int) Math.ceil(i2sigma.get(j)[gp2i.get(j)[ii]]);

                    PointRoi p = new PointRoi(xx+.5, yy+.5);
                    p.setPosition(j + 1);
                    p.setName(Integer.toString(ii));
                    ov.add(p);

                    OvalRoi ovr = new OvalRoi(xx+.5-rr, yy+.5-rr, 2*rr, 2*rr);
                    ovr.setFillColor(new Color(1, 0, 0, 0.3f));
                    ovr.setPosition(j + 1);
                    ov.add(ovr);

                    Line l = new Line(xx+.5, yy+.5, xx+.5+5*i2dir.get(j)[gp2i.get(j)[ii]][0], yy+.5+5*i2dir.get(j)[gp2i.get(j)[ii]][1]);
                    l.setStrokeColor(Color.RED);
                    l.setPosition(j + 1);
                    ov.add(l);

                    l = new Line(xx+.5, yy+.5, xx+.5-5*i2dir.get(j)[gp2i.get(j)[ii]][0], yy+.5-5*i2dir.get(j)[gp2i.get(j)[ii]][1]);
                    l.setStrokeColor(Color.RED);
                    l.setPosition(j + 1);
                    ov.add(l);

                }

            }

        }

        ImagePlus iout = new ImagePlus("queue_map", isqueue);
        iout.setOverlay(ov);
        return iout;
    }

    public ImagePlus getCorrelations(ArrayList<int[][]> i2xy, int W, int H){

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
