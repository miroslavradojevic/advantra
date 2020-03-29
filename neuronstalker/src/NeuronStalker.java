import ij.*;
import ij.gui.*;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import ij.plugin.Text;
import ij.process.FloatProcessor;
import javafx.scene.paint.*;

import java.awt.*;
import java.awt.Color;
import java.io.*;
import java.util.*;

/**
 * Created by miroslav on 8-3-15.
 */
public class NeuronStalker implements PlugIn {

    String      inimg_path;
    ImagePlus   inimg;
    String      inimg_dir_path;
    String      inimg_name;
    String      output_dir_midresults;

    int W,H,L;

    float       neuron_diameter;
    float       percentile;
    float       correlation_boundary;
    float       std_angle_deg;
    float       std_gcsstd_pix;
    int         Ni;
    int         Ns;
    float       step = 1.5f; // fix this one
    float       zDist;

    int         search_mode;
    boolean     save_midresults = false;
    boolean     include_loose_traces = false;
    boolean     include_loose_somas = false;
    boolean     doCuda;
    int         block_size;
    int         max_threads;

    ArrayList<Node>                 stalker         = new ArrayList<Node>();
    int[]                           gpLab2NodeLab   = new int[1]; // filled up during stalking2d, mapping of the guidepoint labels to node labels in hte stalker

    ArrayList<Trace>                trace_list  = new ArrayList<Trace>();
    Overlay                         stalker_gps = new Overlay();
    int  max_tag_disc = 4; // if the trace reached the tag that was max_tag_disc before then exclude due to over tracing or closed loop

    public void run(String s) {

        /**
         * load the image
         */
        String in_folder = Prefs.get("id.folder", System.getProperty("user.home"));
        OpenDialog.setDefaultDirectory(in_folder);
        OpenDialog dc = new OpenDialog("Select image file");
        in_folder = dc.getDirectory();
        inimg_path = dc.getPath();
        if (inimg_path==null) return;
        Prefs.set("id.folder", in_folder);

        inimg = new ImagePlus(inimg_path);
        if (inimg.getType()!=ImagePlus.GRAY8) return;
        inimg.setTitle("inimg");

        W = inimg.getWidth();
        H = inimg.getHeight();
        L = inimg.getNSlices();

        inimg_dir_path      = inimg.getOriginalFileInfo().directory;
        inimg_name          = inimg.getShortTitle();
        inimg_name          = Toolbox.removeExtension(inimg.getOriginalFileInfo().fileName);


        /**
         * load the parameters
         */

        String[] mode = new String[]{"SOMA_CLOSE", "SOMA_AWAY", "HIGH_CORR"};

        if (Macro.getOptions()==null) {

            neuron_diameter         = (float)   Prefs.get("neuronstalker.neuron_diameter", 8f);
            percentile              = (float)   Prefs.get("neuronstalker.percentile", 85f);
            correlation_boundary    = (float)   Prefs.get("neuronstalker.correlation_boundary", 0.75f);
            std_angle_deg           = (float)   Prefs.get("neuronstalker.std_angle_deg", 70f);
            std_gcsstd_pix          = (float)   Prefs.get("neuronstalker.std_gcsstd_pix", 3f);
            Ni                      = (int)     Prefs.get("neuronstalker.Ni", 100);
            Ns                      = (int)     Prefs.get("neuronstalker.Ns", 50);
            zDist                   = (float)   Prefs.get("neuronstalker.zDist", 2f);
            search_mode             = (int)     Prefs.get("neuronstalker.search_mode", 2);
            save_midresults         = (boolean) Prefs.get("neuronstalker.save_midresults", false);
            include_loose_traces    = (boolean) Prefs.get("neuronstalker.include_loose_traces", false);
            include_loose_somas     = (boolean) Prefs.get("neuronstalker.include_loose_somas", false);
            doCuda                  = (boolean) Prefs.get("neuronstalker.do_cuda", true);
            block_size              = (int)     Prefs.get("neuronstalker.cuda_block_size", 1024);
            max_threads             = (int)     Prefs.get("neuronstalker.cuda_max_threads", 20*1024);
            GenericDialog gd = new GenericDialog("NeuronStalker");

            gd.addNumericField("neuron_diameter",       neuron_diameter, 0);
            gd.addNumericField("percentile",            percentile, 0);
            gd.addNumericField("correlation_boundary",  correlation_boundary, 2);
            gd.addNumericField("std_angle_deg",         std_angle_deg, 0);
            gd.addNumericField("std_gcsstd_pix",        std_gcsstd_pix, 0);
            gd.addNumericField("Ni", Ni, 0);
            gd.addNumericField("Ns", Ns, 0);
            gd.addNumericField("zDist", zDist, 1);
            gd.addChoice("search_mode", mode, mode[search_mode]);
            gd.addCheckbox("save_midresults", save_midresults);
            gd.addCheckbox("include_loose_traces", include_loose_traces);
            gd.addCheckbox("include_loose_somas", include_loose_somas);
            gd.addCheckbox("do_cuda", doCuda);
            gd.addNumericField("cuda_block_size", block_size, 0);
            gd.addNumericField("cuda_max_threads", max_threads, 0);

            gd.showDialog();
            if (gd.wasCanceled()) return;

            neuron_diameter = (float) gd.getNextNumber();       Prefs.set("neuronstalker.neuron_diameter", neuron_diameter);            //System.out.println(neuron_diameter);
            percentile = (float) gd.getNextNumber();            Prefs.set("neuronstalker.percentile", percentile);                      //System.out.println(percentile);
            correlation_boundary = (float) gd.getNextNumber();  Prefs.set("neuronstalker.correlation_boundary", correlation_boundary);  //System.out.println(correlation_boundary);
            std_angle_deg = (float) gd.getNextNumber();         Prefs.set("neuronstalker.std_angle_deg", std_angle_deg);                //System.out.println(std_angle_deg);
            std_gcsstd_pix = (float) gd.getNextNumber();        Prefs.set("neuronstalker.std_gcsstd_pix", std_gcsstd_pix);              //System.out.println(std_gcsstd_pix);
            Ni = (int) gd.getNextNumber();                      Prefs.set("neuronstalker.Ni", Ni);                                      //System.out.println(Ni);
            Ns = (int) gd.getNextNumber();                      Prefs.set("neuronstalker.Ns", Ns);                                      //System.out.println(Ns);
            zDist = (float) gd.getNextNumber();                 Prefs.set("neuronstalker.zDist", zDist);                                //System.out.println(zDist);
            search_mode = gd.getNextChoiceIndex();              Prefs.set("neuronstalker.search_mode", search_mode);                    //System.out.println(search_mode);
            save_midresults = gd.getNextBoolean();              Prefs.set("neuronstalker.save_midresults", save_midresults);            //System.out.println(save_midresults);
            include_loose_traces = gd.getNextBoolean();         Prefs.set("neuronstalker.include_loose_traces", include_loose_traces);  //System.out.println(include_loose_traces);
            include_loose_somas  = gd.getNextBoolean();         Prefs.set("neuronstalker.include_loose_somas", include_loose_somas);    //System.out.println(include_loose_somas);
            doCuda = gd.getNextBoolean();                       Prefs.set("neuronstalker.do_cuda", doCuda);
            block_size = (int) gd.getNextNumber();              Prefs.set("neuronstalker.cuda_block_size", block_size);
            max_threads = (int) gd.getNextNumber();             Prefs.set("neuronstalker.cuda_max_threads", max_threads);
        }
        else {

            neuron_diameter         = Float.valueOf(Macro.getValue(Macro.getOptions(),  "neuron_diameter", String.valueOf(8)));
            percentile              = Float.valueOf(Macro.getValue(Macro.getOptions(),  "percentile", String.valueOf(85)));
            correlation_boundary    = Float.valueOf(Macro.getValue(Macro.getOptions(),  "correlation_boundary", String.valueOf(0.75)));
            std_angle_deg           = Float.valueOf(Macro.getValue(Macro.getOptions(),  "std_angle_deg", String.valueOf(70)));
            std_gcsstd_pix          = Float.valueOf(Macro.getValue(Macro.getOptions(),  "std_gcsstd_pix", String.valueOf(3)));
            Ni                      = Integer.valueOf(Macro.getValue(Macro.getOptions(),"Ni", String.valueOf(100)));
            Ns                      = Integer.valueOf(Macro.getValue(Macro.getOptions(),"Ns", String.valueOf(50)));
            zDist                   = Float.valueOf(Macro.getValue(Macro.getOptions(),  "zDist", String.valueOf(2)));

            String read_mode        = Macro.getValue(Macro.getOptions(),"search_mode", mode[2]);
            if (read_mode.equalsIgnoreCase(mode[0]))        search_mode = 0;
            else if (read_mode.equalsIgnoreCase(mode[1]))   search_mode = 1;
            else if (read_mode.equalsIgnoreCase(mode[2]))   search_mode = 2;
            else return;

            save_midresults         = false;
            include_loose_traces    = Boolean.valueOf(Macro.getValue(Macro.getOptions(), "include_loose_traces", String.valueOf(true)));
            include_loose_somas     = Boolean.valueOf(Macro.getValue(Macro.getOptions(), "include_loose_somas", String.valueOf(true)));
            doCuda                  = Boolean.valueOf(Macro.getValue(Macro.getOptions(), "do_cuda",         String.valueOf(true)));
            block_size              = Integer.valueOf(Macro.getValue(Macro.getOptions(), "cuda_block_size", String.valueOf(1024)));
            max_threads             = Integer.valueOf(Macro.getValue(Macro.getOptions(), "cuda_max_threads",String.valueOf(1024*20)));

        }

        output_dir_midresults = inimg_dir_path + inimg_name + "_midresults" + File.separator;
        File f = new File(output_dir_midresults);
        if (!f.exists() && save_midresults) f.mkdirs();

        long t1 = System.currentTimeMillis();

        /**
         * initialize classes used for reconstruction stages
         */
        SomaExtractor se = new SomaExtractor();
        ForegroundExtractor fe = new ForegroundExtractor();
        GuidepointExtractor gpe = new GuidepointExtractor();

        /**
         * extract soma(s)
         */
        float open_radius = .4f*neuron_diameter;
        System.out.print("_____________________\nextracting soma(s) (R=" + IJ.d2s(open_radius,1) + ")... ");
        se.work(inimg, open_radius, (int) neuron_diameter, (int) (0.25*inimg.getWidth()*inimg.getHeight()), save_midresults?output_dir_midresults:"");
        System.out.println(se.soma_list.size() + " soma(s) found");

        if (L==1) { // 2d

            float[][] inimg_xy = new float[W][H]; 	// x~column, y~row
            byte[] read = (byte[]) inimg.getProcessor().getPixels();
            for (int idx=0; idx<read.length; idx++) inimg_xy[idx % W][idx / W] = (float) (read[idx] & 0xff);

            /**
             * extract foreground
             */
            System.out.print("_____________________\nextracting foreground (R=" + IJ.d2s(neuron_diameter, 0) + ")... ");
            fe.work2D(inimg_xy, (int) Math.ceil(Math.sqrt(2) * neuron_diameter), neuron_diameter, percentile, save_midresults ? output_dir_midresults : "");
            System.out.println(IJ.d2s(fe.i2xy.get(0).length/1000f, 0) + "k foreground locs, " + IJ.d2s((100f * fe.i2xy.get(0).length) / (float) (W * H * L), 1) + "%");

            /**
             * extract guidepoints
             */
            System.out.println("_____________________\nextracting guidepoints");
            gpe.work2D(inimg_xy, .5f * neuron_diameter, fe.i2xy.get(0), se.soma_list, correlation_boundary, search_mode, save_midresults ? output_dir_midresults : "",
                    doCuda, block_size, max_threads);
            System.out.println(IJ.d2s(gpe.gp2i.get(0).length/1000f,2) + "k guidepoints");

            /**
             * stalking in 2d...will build up the list of neuron nodes, "stalker"
             */
            System.out.println("_____________________\nstalking... ");
            stalking2d(se, fe, gpe, inimg_xy);
            System.out.println(stalker.size() + " nodes");

            ImagePlus iii = inimg.duplicate();
            iii.setOverlay(getStalkerNodes());
            if (save_midresults) IJ.saveAs(iii, "Tiff", output_dir_midresults+"nodes.tif");

        }
        else { // 3d

            ArrayList<float[][]> inimg_xyz = new ArrayList<float[][]>();    // x~column, y~row, z~slice
            for (int i = 0; i < L; i++) {                                   // loop through the number of slices
                byte[] read = (byte[]) inimg.getStack().getProcessor(i+1).getPixels();
                float[][] tt = new float[W][H];
                for (int idx=0; idx<read.length; idx++) tt[idx % W][idx / W] = (float) (read[idx] & 0xff);
                inimg_xyz.add(tt);
            }

            fe.work3D(inimg_xyz, (int) Math.ceil(Math.sqrt(2)*neuron_diameter), neuron_diameter, percentile, output_dir_midresults);
            int nr_locs = 0;
            for (int i = 0; i < fe.i2xy.size(); i++) nr_locs += fe.i2xy.get(i).length;
            System.out.println(IJ.d2s(nr_locs/1000f, 0) + "k locations of foreground, " + IJ.d2s((100f*nr_locs) / (float) (W * H * L), 1) + "%");

            gpe.work3D(inimg_xyz, zDist, .5f*neuron_diameter, fe.i2xy, se.soma_list, correlation_boundary, search_mode, output_dir_midresults,
                    doCuda, block_size, max_threads);
            if (true) {System.out.println("OUT3D!!");return;}
            /**
             * stalking in 3d...will build up the list of neuron nodes, "stalker"
             */
            System.out.println("_____________________\nstalking... ");
//            stalking3d(se, fe, gpe, inimg_xyz);
            System.out.println(stalker.size() + " nodes");

        }

        // all the info on digital reconstruction are stored in ArrayList<Node> stalker
        // connections, nodes... stalker is used to generate the tree, first as the list of traces and then as a .swc format

        /**
         *  build tree using breath-first search
         *  SomaExtractor extracted list of somas = individual nodes = individual traces
         *  they are added at the beginning of the list ot Nodes (in stalker)
         *  they are added at the beginning of the list of Traces (in trace_list)
         */
        System.out.print("_____________________\ngenerating tree... ");
        // first group of trace elements of the trace list (list of tree elements) are SOMA's and they are taken as seeds of the tree building scheme,
        // stalker.get(1),...,stalker.get(nr_soma) -> trace(0), trace(1), ..., trace(nr_soma-1)
        // if there are no somas nothing will be appended
        ArrayList<Integer> seed_indexes = new ArrayList<Integer>(se.soma_list.size());
        for (int i = 0; i <se.soma_list.size(); i++) seed_indexes.add(i+1); // indexes in stalker node list

//        System.out.println("before tree_build(), nbrNodes of those nodes that corr to soma:");
//        for (int i = 0; i < se.soma_list.size(); i++) {
//            System.out.println("SOMA " + i  + " neighbours: " + stalker.get(i+1).nbrNode + stalker.get(i+1).nbrNode.size());
//        }
//        System.out.println("))))))");
        
        ArrayList[] discovered = init_discovered_list(stalker);
        trace_list = build_tree(stalker, discovered, seed_indexes, Trace.SOMA, Trace.AXON, Trace.LOOSE); // start with soma as seeds
        System.out.println(trace_list.size() + " traces.");

        /**
         * exporting .swc reconstruction
         **/
        System.out.print("_____________________\nswc export... ");

        File swc_dir = new File(inimg_dir_path +
                String.format("NS.D.perc.corr.ang.gcs.Ni.zDist.mode_%.0f_%.0f_%.2f_%.0f_%.0f_%d_%.1f_%d",
                        neuron_diameter, percentile, correlation_boundary, std_angle_deg, std_gcsstd_pix, Ni, zDist, search_mode));
        if (!swc_dir.exists()) swc_dir.mkdirs();
        String swc_path = swc_dir.getAbsolutePath() + File.separator + inimg_name + ".swc";


        // also each of those seed trace elements can happen that no other trace (and hence node) reached them - so they are isolated seed points (somas)
        // stalker can have that info, store it in boolean list
        ArrayList<Integer> loose_soma_idxs = new ArrayList<Integer>();
        for (int i = 0; i < se.soma_list.size(); i++)
            if (stalker.get(i + 1).nbrNode.size() == 0) loose_soma_idxs.add(i); // because they will be later used to index traces and traces start from zero, while nodes start from 1

        exportSwc(trace_list, loose_soma_idxs, include_loose_traces, include_loose_somas, swc_path); // first se.soma_list.size() elements in the trace_list will correspond to soma elements
        System.out.println(swc_path);

        System.out.println("DONE.");

        long t2 = System.currentTimeMillis();
        System.out.println("elapsed " + IJ.d2s((t2 - t1) / 1000f, 1) + " s.");

        if (Macro.getOptions()==null) IJ.setTool("hand");
//        if (true) {System.out.println("OUT2D!!");return;}

    }

    private int get_undiscovered(ArrayList[] discovered){

        for (int i = 1; i < discovered.length; i++) { // first element in stalker list of nodes is nothing
            for (int j = 0; j < discovered[i].size(); j++) {
                if (!(Boolean) discovered[i].get(j)) {
                    return i;
                }
            }
        }

        return -1; // all are discovered

    }

    private ArrayList[] init_discovered_list(ArrayList<Node> _stalker) {
        // will be used to book-keep discovered traces (graph vertices), all are not discovered at the beginning
        ArrayList[] discovered = new ArrayList[_stalker.size()];
        for (int i = 1; i < _stalker.size(); i++) { // skip the first index as the first element of stalker is dummy
            discovered[i] = new ArrayList<Boolean>(_stalker.get(i).nbrNode.size());
            for (int j = 0; j < _stalker.get(i).nbrNode.size(); j++) {
                discovered[i].add(false);
            }
        }
        return discovered;
    }

    private ArrayList<Trace> build_tree(
            ArrayList<Node> _stalker,
            ArrayList[] discovered,
            ArrayList<Integer> _seed_indexes,
            int seed_trace_type,
            int tree_trace_type,
            int loose_trace_type)
    {

        /**
         *  breadth-first search (BFS) to traverse the tree from extracted node list
         *  http://en.wikipedia.org/wiki/Breadth-first_search
         *
         1  procedure BFS(G,v) is
         2      let Q be a queue
         3      Q.enqueue(v)
         4      label v as discovered
         5      while Q is not empty
         6         v ← Q.dequeue()
         7         for all edges from v to w in G.adjacentEdges(v) do
         8             if w is not labeled as discovered
         9                 Q.enqueue(w)
         10                label w as discovered
         *
         */

        // will essentially convert ArrayList<Node> with all its linking (2 directional connections)
        // into ArrayList<Trace> which contains the list of traces (1 directional connections)
        // BFS needs Queue data structure that will be implemented in class BfsQueue
        // Queue keeps the links between the nodes, link is described as int[] where int[] ~ [curr_node_idx, adjacent_node_idx]
        // knowing just the node index in the queue is not enough, need to know the index of the mother node
        // discovered is arry of lists that will keep the labels of the discovered adjacent node pairs, bookekeeping for the BFS

        // the key reason for keeping two values is that each time we need the mother index and so for each trace element, including the first node of the trace
        boolean is3D = _stalker.get(1).loc.length==3;

        BfsQueue bfsQueue = new BfsQueue();

        ArrayList<Trace> tlist = new ArrayList<Trace>(); // output trace list

        for (int i = 0; i < _seed_indexes.size(); i++) { // add seed indexes to queue, mark the as discovered

            int sidx = _seed_indexes.get(i);

            Trace t = new Trace(seed_trace_type); // init seed trace

            float x_node = _stalker.get(sidx).loc[0];
            float y_node = _stalker.get(sidx).loc[1];
            float z_node = (is3D)?_stalker.get(sidx).loc[2]:0;
            float r_node = _stalker.get(sidx).r;

            t.add(-1, sidx, x_node, y_node, r_node);
            tlist.add(t); // add the seed points as 1 element traces

            // add the neighbors to the queue and label them as discovered
            for (int j = 0; j <_stalker.get(sidx).nbrNode.size(); j++) {
                int next = _stalker.get(sidx).nbrNode.get(j);
                // enqueue(), add to FIFO structure, http://en.wikipedia.org/wiki/Queue_%28abstract_data_type%29
                bfsQueue.enqueue(new int[]{sidx, next});
                discovered[sidx].set(j, true);                                          // set label to discovered in both neighbouting index lists
                discovered[next].set(_stalker.get(next).nbrNode.indexOf(sidx), true);   // index where the background link was found
            }
        }

        System.out.println("BFS expands from " + _seed_indexes.size() + " points. |Q| = " + bfsQueue.size() + " elements.");

        while (bfsQueue.hasItems()) {

//            System.out.println(bfsQueue.size() + " elements in queue.");

            // dequeue(), take from FIFO structure, http://en.wikipedia.org/wiki/Queue_%28abstract_data_type%29
            int [] getLnk = (int[]) bfsQueue.dequeue();

//            System.out.println(bfsQueue.size() + " elements in queue.");

            // next neighbour at the time it was added to the queue becomes current
            int prev = getLnk[0]; //Q.get(Q.size()-1)[0];
            int curr = getLnk[1]; //Q.get(Q.size()-1)[1];

            Trace t = new Trace(tree_trace_type); // init branch trace

            // always add the first node (it exists since this one was stored in the queue)
            t.add(
                    prev, curr,
                    _stalker.get(curr).loc[0], _stalker.get(curr).loc[1], (is3D)?_stalker.get(curr).loc[2]:0,
                    _stalker.get(curr).r
            );

            while(Collections.frequency(discovered[curr], false)==1) { // while the number of undiscovered is 1 just step further // _stalker.get(curr).nbrNode.size()==2 // while it is a 'BODY' point

//                System.out.print("."+curr+"("+prev+"). " + discovered[curr] + "|");
                prev = curr;
                curr = _stalker.get(curr).nbrNode.get(discovered[curr].indexOf(false)); // curr takes the value of the next step
//                System.out.print(" "+curr + "   ");

                t.add(
                        prev, curr,
                        _stalker.get(curr).loc[0], _stalker.get(curr).loc[1], (is3D)?_stalker.get(curr).loc[2]:0,
                        _stalker.get(curr).r
                );

                // mark as discovered the connections curr--prev and prev--curr
                discovered[curr].set(_stalker.get(curr).nbrNode.indexOf(prev), true);
                discovered[prev].set(_stalker.get(prev).nbrNode.indexOf(curr), true);

            }

            tlist.add(t);

            // means that it was not on the neurite anymore, check adjacent traces, now it is either enpoint or junction
            for (int i = 0; i < discovered[curr].size(); i++) {
                boolean isDiscovered =  (Boolean) discovered[curr].get(i);
                if (!isDiscovered) { // if it was not discovered
                    int next = _stalker.get(curr).nbrNode.get(i);

                    bfsQueue.enqueue(new int[]{curr, next});   // enqueue()

                    discovered[curr].set(i, true); // label as discovered
                    discovered[next].set(_stalker.get(next).nbrNode.indexOf(curr), true);

//                    System.out.println(bfsQueue.size() + " elements in queue (reached point with !=2 neighbours).");
                }
            }

        }

        // before returning do bfs from those  that were not discovered, seed each time from an undiscovered one until all are discovered
        // those that had undiscovered links
        int sidx = get_undiscovered(discovered);

        while (sidx != -1) {

            // do a new BFS traverse from the first undiscovered seed till those exist
            Trace t = new Trace(loose_trace_type); // init seed trace

            float x_node = _stalker.get(sidx).loc[0];
            float y_node = _stalker.get(sidx).loc[1];
            float z_node = (is3D)?_stalker.get(sidx).loc[2]:0;
            float r_node = _stalker.get(sidx).r;

            t.add(-1, sidx, x_node, y_node, r_node);

            tlist.add(t);

            // add the neighbors to the queue and label them as discovered
            for (int j = 0; j <_stalker.get(sidx).nbrNode.size(); j++) {

                if ((Boolean) discovered[sidx].get(j))
                    System.out.println("ERROR: isolated point HAD at least one discovered link");

                int next = _stalker.get(sidx).nbrNode.get(j);

                bfsQueue.enqueue(new int[]{sidx, next});
                discovered[sidx].set(j, true);
                discovered[next].set(_stalker.get(next).nbrNode.indexOf(sidx), true);

//                System.out.println(bfsQueue.size() + " elements in queue.");

            }

            System.out.println("BFS expands from isolated node #" + sidx + "  |Q| = " + bfsQueue.size() + " elements.");

            while (bfsQueue.hasItems()) {

                int[] getLnk = (int[]) bfsQueue.dequeue();
                int prev = getLnk[0]; //Q.get(Q.size()-1)[0];
                int curr = getLnk[1]; //Q.get(Q.size()-1)[1];

//                System.out.println(bfsQueue.size() + " elements in queue.");

                t = new Trace(loose_trace_type); // init branch trace

                // always add the first node (it exists since this one was stored in the queue)
                t.add(
                        prev, curr,
                        _stalker.get(curr).loc[0], _stalker.get(curr).loc[1], (is3D)?_stalker.get(curr).loc[2]:0,
                        _stalker.get(curr).r
                );

                while(Collections.frequency(discovered[curr], false)==1) { // _stalker.get(curr).nbrNode.size()==2

                    prev = curr;
                    curr = _stalker.get(curr).nbrNode.get(discovered[curr].indexOf(false)); // curr takes the value of the next step

                    t.add(
                            prev, curr,
                            _stalker.get(curr).loc[0], _stalker.get(curr).loc[1], (is3D)?_stalker.get(curr).loc[2]:0,
                            _stalker.get(curr).r
                    );

                    // mark as discovered the connections curr--prev and prev--curr
                    discovered[curr].set(_stalker.get(curr).nbrNode.indexOf(prev), true);
                    discovered[prev].set(_stalker.get(prev).nbrNode.indexOf(curr), true);

                }

                tlist.add(t);

                for (int i = 0; i < discovered[curr].size(); i++) {
                    boolean isDiscovered =  (Boolean) discovered[curr].get(i);
                    if (!isDiscovered) { // if it was not discovered

                        int next = _stalker.get(curr).nbrNode.get(i);

                        bfsQueue.enqueue(new int[]{curr, next}); // enqueue

                        discovered[curr].set(i, true);  // label as discovered
                        discovered[next].set(_stalker.get(next).nbrNode.indexOf(curr), true);

                    }
                }

            }

            sidx = get_undiscovered(discovered);

        }

        return tlist;

    }

    private ImagePlus getTagMap(int[] tag_map, int w, int h) {

        return new ImagePlus("tag_map", new FloatProcessor(w, h, tag_map.clone()));

    }

    private ImagePlus getTagMap(ArrayList<int[]> tag_map, int w, int h) {

        ImageStack isout = new ImageStack(w, h);

        for (int i = 0; i < tag_map.size(); i++) {
            isout.addSlice(new FloatProcessor(w, h, tag_map.get(i)));
        }

        return new ImagePlus("tag_map", isout);

    }

    private void addToMap2D(int _tag, int[] _tag_map, int _mapW, int _mapH, float atx, float aty, float atr) {

        int rr = (int) Math.ceil(atr);
        int xx_min = (int) Math.floor(atx-rr);
        int xx_max = (int) Math.ceil( atx+rr);
        int yy_min = (int) Math.floor(aty - rr);
        int yy_max = (int) Math.ceil(aty + rr);

        for (int xx = xx_min; xx <= xx_max; xx++) {
            for (int yy = yy_min; yy <= yy_max; yy++) {
                if (xx>=0 && xx<_mapW && yy>=0 && yy<_mapH) {
                    if (  Math.pow(xx-atx,2) + Math.pow(yy-aty,2) <= Math.pow(rr,2) ) {
                        _tag_map[yy * _mapW + xx] = _tag;
                    }
                }
            }
        }

    }

    private void addToMap2D(int _tag, int[] _tag_map, int _mapW, int _mapH, int _gp_lab, int[][] _gp_lab_map) { // 2d

        float curr_x = stalker.get(_tag).loc[0];
        float curr_y = stalker.get(_tag).loc[1];
        int curr_r = (int) Math.ceil(stalker.get(_tag).r);

        int xx_min = (int) Math.floor(curr_x - curr_r);
        int xx_max = (int) Math.ceil(curr_x + curr_r);

        int yy_min = (int) Math.floor(curr_y - curr_r);
        int yy_max = (int) Math.ceil(curr_y + curr_r);

        // this one brings in only new node labels to the map, not the gp labels
        // do not write over guidepoint labels, but allow writing over node labels
        for (int xx = xx_min; xx <= xx_max; xx++) {
            for (int yy = yy_min; yy <= yy_max; yy++) {
                    if (xx>=0 && xx<_mapW && yy>=0 && yy<_mapH) {
                        if (  Math.pow(xx-curr_x,2) + Math.pow(yy-curr_y,2) <= Math.pow(curr_r,2) ) {

                            if (_gp_lab_map[xx][yy]==_gp_lab || _gp_lab_map[xx][yy]==0) { // block writing if guidepoint label is diferent than the one we started with
                                _tag_map[yy * _mapW + xx] = _tag;
                            }

                        }
                    }
            }
        }

    }

    private void addToMapWithLinking2D(int _tag, int[] _tag_map, int _mapW, int _mapH, int _gp_lab, int[][] _gp_lab_map) { // 2d

        float curr_x = stalker.get(_tag).loc[0];
        float curr_y = stalker.get(_tag).loc[1];
        int curr_r = (int) Math.ceil(stalker.get(_tag).r);

        int xx_min = (int) Math.floor(curr_x - curr_r);
        int xx_max = (int) Math.ceil(curr_x + curr_r);

        int yy_min = (int) Math.floor(curr_y - curr_r);
        int yy_max = (int) Math.ceil(curr_y + curr_r);

        // this one brings in only new node labels to the map, not the gp labels
        // do not write over guidepoint labels, but allow writing over node labels
        for (int xx = xx_min; xx <= xx_max; xx++) {
            for (int yy = yy_min; yy <= yy_max; yy++) {
                if (xx>=0 && xx<_mapW && yy>=0 && yy<_mapH) {
                    if (  Math.pow(xx-curr_x,2) + Math.pow(yy-curr_y,2) <= Math.pow(curr_r,2) ) {
                        if (_gp_lab_map[xx][yy]==_gp_lab || _gp_lab_map[xx][yy]==0) {
                            _tag_map[yy * _mapW + xx] = _tag; // owerwrite node tag over guidepoint tag
                        }
                    }
                }
            }
        }

        // linking in between the nodes is tagged
        // part that does linking with  the previous node
        float x2 = stalker.get(_tag - 1).loc[0];
        float y2 = stalker.get(_tag - 1).loc[1];
        int r2 = (int) Math.ceil(stalker.get(_tag - 1).r);

        // curr_x is current node, x2 is previous node
        float vnorm = (float) Math.sqrt(Math.pow(x2 - curr_x, 2) + Math.pow(y2 - curr_y, 2));

        if (vnorm<r2+curr_r) return; // if the distance is necessary to be filled up

        float v1x = (x2 - curr_x) / vnorm;
        float v1y = (y2 - curr_y) / vnorm;
        float v2x = (curr_x - x2) / vnorm;
        float v2y = (curr_y - y2) / vnorm;

        xx_min = (int) Math.floor(Math.min(curr_x - curr_r, x2 - r2));
        xx_max = (int) Math.ceil(Math.max(curr_x + curr_r, x2 + r2));
        yy_min = (int) Math.floor(Math.min(curr_y - curr_r, y2 - r2));
        yy_max = (int) Math.ceil(Math.max(curr_y + curr_r, y2 + r2));

        for (int xx = xx_min; xx <= xx_max; xx++) {
            for (int yy = yy_min; yy <= yy_max; yy++) {
                if (xx >= 0 && xx < _mapW && yy >= 0 && yy < _mapH) {

                    float p1x = xx - curr_x;
                    float p1y = yy - curr_y;

                    float side1 = p1x * v1x + p1y * v1y;

                    float p2x = xx - x2;
                    float p2y = yy - y2;

                    float side2 = p2x * v2x + p2y * v2y;

                    float p2lx = -p1x + side1 * v1x;
                    float p2ly = -p1y + side1 * v1y;

                    float p2l_sqr = (p2lx * p2lx + p2ly * p2ly);

                    if (side1 > 0 && side2 > 0 && p2l_sqr <= Math.pow(Math.min(curr_r, r2), 2)) {
                        if (_tag_map[yy * _mapW + xx]==0) { // if there was untagged location - tag it to connec with previous node
                            _tag_map[yy * _mapW + xx] = _tag;
                        }
                    }
                }
            }
        }






    }

    private void addToMap(int _tag, int[] _tag_map, int _mapW, int _mapH) { // 2d

        float curr_x = stalker.get(_tag).loc[0];
        float curr_y = stalker.get(_tag).loc[1];
        int curr_r = (int) Math.ceil(stalker.get(_tag).r);

        int xx_min = (int) Math.floor(curr_x - curr_r);
        int xx_max = (int) Math.ceil(curr_x + curr_r);

        int yy_min = (int) Math.floor(curr_y - curr_r);
        int yy_max = (int) Math.ceil( curr_y + curr_r);

        for (int xx = xx_min; xx <= xx_max; xx++) {
            for (int yy = yy_min; yy <= yy_max; yy++) {
                if (xx>=0 && xx<_mapW && yy>=0 && yy<_mapH) {
                    if (  Math.pow(xx-curr_x,2) + Math.pow(yy-curr_y,2) <= Math.pow(curr_r,2) ) {
                        if (_tag_map[yy * _mapW + xx] == 0) {// fill new value if it has not been filled already
                                _tag_map[yy * _mapW + xx] = _tag;
                        }
                    }
                }
            }
        }

    }

    private int maxTagJump(ArrayList<Node> trace_check, int[] _tag_map, int _mapW, int _mapH) {

        int max_tag_discrepancy = 0;

        for (int i = 0; i < trace_check.size(); i++) {

            float x1 = trace_check.get(i).loc[0];
            float y1 = trace_check.get(i).loc[1];
            int   r1 = (int) Math.ceil(trace_check.get(i).r);

            int xmin = (int) Math.floor( x1-r1);
            int xmax = (int) Math.floor( x1+r1);

            int ymin = (int) Math.floor( y1-r1);
            int ymax = (int) Math.ceil(  y1+r1);

            for (int xx = xmin; xx <= xmax; xx++) {
                for (int yy = ymin; yy <= ymax; yy++) {
                        if (xx>=0 && xx<_mapW && yy>=0 && yy<_mapH) {
                            if ( Math.pow(xx-x1,2) + Math.pow(yy-y1,2) <= Math.pow(r1,2) ) {
                                if (_tag_map[yy*_mapW+xx] == 0) {
                                    _tag_map[yy*_mapW+xx] = i+1;
                                }
                                else {
                                    int tag_discrepancy = (i+1) - _tag_map[yy*_mapW+xx];
                                    if (tag_discrepancy>max_tag_discrepancy) max_tag_discrepancy = tag_discrepancy;
                                }
                            }
                        }
                }
            }

            if (i==0) continue;

            // part that does linking with  the previous node
            float x2 = trace_check.get(i - 1).loc[0];
            float y2 = trace_check.get(i - 1).loc[1];
            int r2 = (int) Math.ceil(trace_check.get(i - 1).r);

            // x1 is current node, x2 is previous node
            float vnorm = (float) Math.sqrt(Math.pow(x2 - x1, 2) + Math.pow(y2 - y1, 2));

            float v1x = (x2 - x1) / vnorm;
            float v1y = (y2 - y1) / vnorm;
            float v2x = (x1 - x2) / vnorm;
            float v2y = (y1 - y2) / vnorm;


            xmin = (int) Math.floor(Math.min(x1 - r1, x2 - r2));
            xmax = (int) Math.ceil(Math.max(x1 + r1, x2 + r2));
            ymin = (int) Math.floor(Math.min(y1 - r1, y2 - r2));
            ymax = (int) Math.ceil(Math.max(y1 + r1, y2 + r2));

            for (int xx = xmin; xx <= xmax; xx++) {
                for (int yy = ymin; yy <= ymax; yy++) {
                        if (xx >= 0 && xx < _mapW && yy >= 0 && yy < _mapH) {

                            float p1x = xx - x1;
                            float p1y = yy - y1;


                            float side1 = p1x * v1x + p1y * v1y;

                            float p2x = xx - x2;
                            float p2y = yy - y2;

                            float side2 = p2x * v2x + p2y * v2y;

                            float p2lx = -p1x + side1 * v1x;
                            float p2ly = -p1y + side1 * v1y;

                            float p2l = (float) Math.sqrt(p2lx * p2lx + p2ly * p2ly);

                            if (side1 > 0 && side2 > 0 && p2l <= Math.min(r1, r2)) {
                                if (_tag_map[yy * _mapW + xx] == 0) {
                                    _tag_map[yy * _mapW + xx] = i + 1;
                                } else {
                                    int tag_discrepancy = (i + 1) - _tag_map[yy * _mapW + xx];
                                    if (tag_discrepancy > max_tag_discrepancy) max_tag_discrepancy = tag_discrepancy;
                                }
                            }

                        }
                }
            }

        }

        return max_tag_discrepancy;

    }

    private void exportSwc(ArrayList<Trace> trace_list, ArrayList<Integer> loose_soma_idxs, boolean include_loose_traces, boolean include_loose_somas, String swc_path){

        if (trace_list.size()==0) {
            System.out.println("trace list was empty");
            return;
        }

        boolean is3dtrace = trace_list.get(0).locs.get(0).length==3;

        PrintWriter logWriter = null;
        try {
            logWriter = new PrintWriter(swc_path);
            logWriter.print("");
            logWriter.close();
        } catch (FileNotFoundException ex) {}

        try {
            logWriter = new PrintWriter(new BufferedWriter(new FileWriter(swc_path, true)));
            logWriter.println("#NeuronStalker\n#ADVANTRA\n#author: Miro\n# " + trace_list.size() + " traces");
        } catch (IOException e) {}

        exportTrace(trace_list, loose_soma_idxs, logWriter, is3dtrace, include_loose_traces, include_loose_somas);

        logWriter.close();

    }

    private void exportTrace(ArrayList<Trace> trace_list, ArrayList<Integer> loose_soma_idxs, PrintWriter logWriter, boolean is3dtrace,
                             boolean include_loose_traces, boolean include_loose_somas){

        // table mapping Node index into Swc index (which increment as the traces are looped and always concatenate on previous traces)

        int[] NodeId2SwcId = new int[stalker.size()];  NodeId2SwcId[0] = Integer.MIN_VALUE;

        int swc_id = 1;

        for (int i = 0; i < trace_list.size(); i++) {

            Trace tt = trace_list.get(i);

            if (tt.type==Trace.LOOSE && !include_loose_traces) continue; // skip adding elements of this trace

            if (loose_soma_idxs.contains(i) && !include_loose_somas) continue; // skip disconnected somas

            for (int j = 0; j < tt.locs.size(); j++) {

//                logWriter.println(String.format(
//                        "%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d",
//                        tt.curr.get(j), tt.type, tt.locs.get(j)[0], tt.locs.get(j)[1], (is3dtrace)?tt.locs.get(j)[2]:0f, tt.rads.get(j), tt.prev.get(j)));

                int NodeIdx = tt.curr.get(j);
                NodeId2SwcId[NodeIdx] = swc_id;

                logWriter.println(String.format(
                        "%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d",
                        swc_id, tt.type, tt.locs.get(j)[0], tt.locs.get(j)[1], (is3dtrace)?tt.locs.get(j)[2]:0f, tt.rads.get(j),
                        (tt.prev.get(j)==-1)?-1:NodeId2SwcId[tt.prev.get(j)])); // by definition traces are generated so that prev is always something that was added to the swc tree earlier
                swc_id++;

            }

        }

    }

    private void stalkNode(Node soma_node, int[] node2lab) // used for soma,if it does not initialize with any guidepoint
    {
        stalker.add(soma_node);
        addToMap(stalker.size() - 1, node2lab, W, H);
    }

    private void stalkLinker(int node_label) {
        stalker.get(node_label      ).nbrNode.add(stalker.size()-1);
        stalker.get(stalker.size()-1).nbrNode.add(node_label);
    }

    private void gpntLinker(int gpnt_label) {
        stalker.get(stalker.size()-1).nbrGpnt.add(gpnt_label); // this one is yet to be added to node list and get a label, 1 direction only
        // still need to connect from the guidepoint to the current node
    }

    private void stalking2d(
            SomaExtractor           se,
            ForegroundExtractor     fe,
            GuidepointExtractor     gpe,
            float[][]               inimg_xy
    )
    {

        stalker.clear();
        stalker.add(new Node());// corresponding index in the list is initialized with a dummy node

        int nr_gudepoints = 0;
        for (int i = 0; i < gpe.gp2i.size(); i++) nr_gudepoints += gpe.gp2i.get(i).length;
        
        gpLab2NodeLab = new int[nr_gudepoints+1]; // 2d there are as many as there are indexes in the only layer of gpe.gp2i
        gpLab2NodeLab[0] = Integer.MIN_VALUE; // as map labels will be zero for the background, they are i+1

        stalker_gps.clear();

        int[] queue             = gpe.gp2i.get(0); // indexes of the foreground locations in i2xy
        int[][] gp_label_map    = gpe.gp2lab.get(0);
        int[] node2lab          = new int[W*H]; // map of labels assigned to each node during stalking

        // soma addition
        for (int i = 0; i < se.soma_list.size(); i++) {

            float xx = (float) se.soma_list.get(i)[0];
            float yy = (float) se.soma_list.get(i)[1];
            float rr = (float) se.soma_list.get(i)[2];

            stalkNode(new Node(xx, yy, rr), node2lab);

        }

        // bayesian tracing in 2d will give the list of nodes (don't add to the node label map till the trace is finished)
        BayesianTracer2D bt = new BayesianTracer2D(Ni, Ns, step, .5f*neuron_diameter);
        int[] t1_NODE_GDPNT = new int[2];
        boolean t1_reached_NODE = false;
        boolean t1_reached_GDPNT = false;
        boolean t1_completed = false;

        int[] t2_NODE_GDPNT = new int[2];
        boolean t2_reached_NODE = false;
        boolean t2_reached_GDPNT = false;
//        boolean t2_completed = false;

        int[] trace_tag_map = new int[W*H]; // to check the tag discrepancy

            int LOG_EVERY = queue.length / 10;

            for (int qidx = 0; qidx < queue.length; qidx++) {

                if (qidx % LOG_EVERY == 0) System.out.print(IJ.d2s((qidx/(float)LOG_EVERY)*10,0) + "%\t");

                // take the guidepoint to initialize tracng in 2 directions, queue[i] defines the foreground index
                int start_px         =  fe.i2xy.get(0)[    queue[qidx]][0];
                int start_py         =  fe.i2xy.get(0)[    queue[qidx]][1];
                float start_vx       = gpe.i2dir.get(0)[   queue[qidx]][0];
                float start_vy       = gpe.i2dir.get(0)[   queue[qidx]][1];
                float start_gcsstd   = gpe.i2sigma.get(0)[ queue[qidx]];

                int start_GpLab = gp_label_map[start_px][start_py]; // take it's unique label

                if (node2lab[start_py * W + start_px] > 0) {
                    // since it shoud not overlap with other guidepoint region, this can only be covered by soma labels already
                    // no need to trace as it is soma, just update the node label table for later, don't add it to Node list but propagate the labels that were there fo soma component

                    addToMap2D(node2lab[start_py * W + start_px], node2lab, W, H, start_px, start_py, start_gcsstd);

                    gpLab2NodeLab[qidx+1] = node2lab[start_py * W + start_px]; // just in case, but adding to map will disable from referring to this guidepoint region at the first place

                    continue; // skip tracing from this guidepoint, soma will swallow this guidepoint and it;s region and it wont be run
                }
                ArrayList<Node> gptrace = new ArrayList<Node>(); // guidepoint trace will be initialized by guidepoint and go in two directions, with 2 parts concatenated together
                int gpstart_index = -1; // index of the guidepoint within the trace

//                System.out.print("<");

                ////////////////////////////////////////////////
                // trace 1
                ////////////////////////////////////////////////
                bt.trace(
                        start_px, start_py, start_vx, start_vy, start_gcsstd,
                        inimg_xy,
                        std_angle_deg, std_gcsstd_pix,
                        node2lab,
                        fe.mask_xy.get(0),
                        gp_label_map,
                        start_GpLab,
                        t1_NODE_GDPNT
                );

//                System.out.print("-");

                t1_reached_NODE     = t1_NODE_GDPNT[0] != Integer.MIN_VALUE;
                t1_reached_GDPNT    = t1_NODE_GDPNT[1] != Integer.MIN_VALUE;
                t1_completed        = t1_reached_NODE || t1_reached_GDPNT;

                if (t1_reached_NODE) {
                    if (bt.iter_counter==0) {
                        gptrace.add(new Node(start_px, start_py, start_gcsstd));
                        gpstart_index = gptrace.size()-1;
                    }
                    else {
                        int LL = bt.iter_counter-1;
                        for (int j = LL; j >= 0; j--) gptrace.add(new Node(bt.xc[j][0], bt.xc[j][1], bt.rc[j]));

                        gptrace.add(new Node(start_px, start_py, start_gcsstd));
                        gpstart_index = gptrace.size()-1;
                    }
                }
                else if (t1_reached_GDPNT) {
                    if (bt.iter_counter==0) {
                        gptrace.add(new Node(start_px, start_py, start_gcsstd));
                        gpstart_index = gptrace.size()-1;
                    }
                    else {
                        int LL = bt.iter_counter-1;
                        for (int j = LL; j >= 0; j--) gptrace.add(new Node(bt.xc[j][0], bt.xc[j][1], bt.rc[j]));

                        gptrace.add(new Node(start_px, start_py, start_gcsstd));
                        gpstart_index = gptrace.size()-1;
                    }
                }

                ////////////////////////////////////////////////
                // trace 2
                ////////////////////////////////////////////////
                bt.trace(start_px, start_py, -start_vx, -start_vy, start_gcsstd,
                        inimg_xy,
                        std_angle_deg,
                        std_gcsstd_pix,
                        node2lab,
                        fe.mask_xy.get(0),
                        gp_label_map,
                        start_GpLab,
                        t2_NODE_GDPNT
                );

//                System.out.print("->");

                t2_reached_NODE     = t2_NODE_GDPNT[0] != Integer.MIN_VALUE;
                t2_reached_GDPNT    = t2_NODE_GDPNT[1] != Integer.MIN_VALUE;

                if (!t1_completed) {
                    gptrace.add(new Node(start_px, start_py, start_gcsstd));
                    gpstart_index = gptrace.size()-1;
                }

                int LL = bt.iter_counter-1;
                for (int j = 0; j <= LL; j++) gptrace.add(new Node(bt.xc[j][0], bt.xc[j][1], bt.rc[j]));

                // gptrace contains the set of Nodes that are added to the stalker list

                // check if there is a closed loop and if so skip adding this trace
                Arrays.fill(trace_tag_map, 0);

                int mtj = maxTagJump(gptrace, trace_tag_map, W, H);
                // add Node, add the label of the guidepoint to the table, add linking, add the trace(list of Nodes) to the map
                if (mtj>max_tag_disc) {
                    System.out.print(" ("+mtj+") ");
                    // only add the existing guidepoint without the trace and any further connections, add it to stalker list
                    stalker.add(new Node(start_px, start_py, start_gcsstd));
                    gpLab2NodeLab[qidx+1] = stalker.size()-1;
                    // no Linker for this point
                    addToMap2D(stalker.size() - 1, node2lab, W, H, start_GpLab, gp_label_map);
                }
                else { // trace was ok

                    for (int tidx = 0; tidx < gptrace.size(); tidx++) {
                        stalker.add(gptrace.get(tidx));
                        if (gpstart_index==tidx) gpLab2NodeLab[qidx+1] = stalker.size()-1;
                        if (tidx==0) {
                            // just link if some of the regions reached
                            if (t1_reached_NODE)        stalkLinker(t1_NODE_GDPNT[0]);
                            else if (t1_reached_GDPNT)  gpntLinker(t1_NODE_GDPNT[1]);

                        }
                        if (tidx==gptrace.size()-1) { // last element

                            if (tidx>0) stalkLinker(stalker.size()-2);

                            if (t2_reached_NODE)        stalkLinker(t2_NODE_GDPNT[0]);
                            else if (t2_reached_GDPNT)  gpntLinker(t2_NODE_GDPNT[1]);
                        }
                        else {
                            if (tidx>0) stalkLinker(stalker.size()-2);
                        }

                        if (tidx==0)    addToMap2D(stalker.size() - 1, node2lab, W, H, start_GpLab, gp_label_map);
                        else            addToMapWithLinking2D(stalker.size() - 1, node2lab, W, H, start_GpLab, gp_label_map);

                    }
                }

            }


//        System.out.println("CHECKING>>>");
//        for (int i = 0; i < gpLab2NodeLab.length; i++) {
//            if (gpLab2NodeLab[i]==0) System.out.println(i + " ::-:: " + gpLab2NodeLab[i]);
//        }
//        System.out.println("done checking");

        complete_neighbourhood(stalker, gpLab2NodeLab);

        if (save_midresults) IJ.saveAs(getTagMap(node2lab, W, H), "Tiff", output_dir_midresults + "tagMap" + String.format("%05d", stalker.size()) + ".tif");

    }

    private void complete_neighbourhood(ArrayList<Node> _stalker, int[] _gpLab2NodeLab) {

        // linking guidepoint label and the stalker index, so that all neighbours are in stalker indexes nbrGpnt->nbrNode
        for (int i = 1; i < _stalker.size(); i++) {
            for (int j = 0; j < _stalker.get(i).nbrGpnt.size(); j++) {
                int curr_node_idx = i;

                int nbgr_node_idx = _gpLab2NodeLab[_stalker.get(i).nbrGpnt.get(j)]; // gp label will be used as the index in the table, value read at the index will be node index

//                System.out.println("guidepoint label" + _stalker.get(i).nbrGpnt.get(j));
//                System.out.println("form table:");
//                System.out.println("nbgr_node_idx = " + nbgr_node_idx);
                // means that curr_node_idx has nbgr_node_idx as a neighbour, and the other way around
                _stalker.get(curr_node_idx).nbrNode.add(nbgr_node_idx);
                _stalker.get(nbgr_node_idx).nbrNode.add(curr_node_idx);
            }
        }

        // remove duplicate neighbourhood links
        for (int i = 1; i < _stalker.size(); i++) {
            Set<Integer> set = new HashSet<Integer>();
            set.addAll(_stalker.get(i).nbrNode);
            _stalker.get(i).nbrNode.clear();
            _stalker.get(i).nbrNode.addAll(set);
        }

        // check if the neighborhoods are conistent
        System.out.print("\nchecking neighbourhood consistency...");
        for (int i = 1; i < _stalker.size(); i++) {
            for (int j = 0; j < _stalker.get(i).nbrNode.size(); j++) {
                int nbr_idx = _stalker.get(i).nbrNode.get(j);
                if (Collections.frequency(_stalker.get(nbr_idx).nbrNode, i)!=1) System.out.println("ERROR: " + i + " -- " + nbr_idx);
            }
        }
        System.out.println("done.");

    }



    private float point2line(float n1x, float n1y,  // float n1z,
                             float n2x, float n2y,  // float n2z,
                             float px,  float py    //, float pz
    )
    {

        float d = 0;

        float[] p_b = new float[2];

        float n21Len = (float) Math.sqrt(Math.pow(n2x-n1x,2)+Math.pow(n2y-n1y,2)); // +Math.pow(n2z-n1z,2)
        float n21x = (n2x-n1x)/n21Len;
        float n21y = (n2y-n1y)/n21Len;
//        float n21z = (n2z-n1z)/n21Len;

        float proj = (px - n1x) * n21x + (py - n1y) * n21y; // + (pz - n1z) * n21z; // dot prod

        p_b[0] = -(px - n1x) + proj * n21x;
        p_b[1] = -(py - n1y) + proj * n21y;
//        p_b[2] = -(pz - n1z) + proj * n21z;

        return (float) Math.sqrt(p_b[0]*p_b[0] + p_b[1]*p_b[1]); // + p_b[2]*p_b[2]

    }

    private float point2line(float n1x, float n1y, float n1z,
                             float n2x, float n2y, float n2z,
                             float px,  float py , float pz
    )
    {

        float d = 0;

        float[] p_b = new float[3];

        float n21Len = (float) Math.sqrt(Math.pow(n2x-n1x,2)+Math.pow(n2y-n1y,2)+Math.pow(n2z-n1z,2));
        float n21x = (n2x-n1x)/n21Len;
        float n21y = (n2y-n1y)/n21Len;
        float n21z = (n2z-n1z)/n21Len;

        float proj = (px - n1x) * n21x + (py - n1y) * n21y + (pz - n1z) * n21z; // dot prod

        p_b[0] = -(px - n1x) + proj * n21x;
        p_b[1] = -(py - n1y) + proj * n21y;
        p_b[2] = -(pz - n1z) + proj * n21z;

        return (float) Math.sqrt(p_b[0]*p_b[0] + p_b[1]*p_b[1] + p_b[2]*p_b[2]); 

    }

    private Overlay getStalkerNodes(){

        // will output the stalker nodes and stalker_gps with groundpoints marked

        Overlay ov = new Overlay();
        ov.drawNames(true);

        for (int i = 1; i < stalker.size(); i++) {

            PointRoi ovr = new PointRoi(stalker.get(i).loc[0]+.5, stalker.get(i).loc[1]+.5);
            ovr.setName(Integer.toString(i));
            ovr.setStrokeColor(Color.RED);
            ov.add(ovr);

            ArrayList<Integer> nbrs = stalker.get(i).nbrNode;
            for (int j = 0; j < nbrs.size(); j++) {

//                boolean chk = false;
//                for (int k = 0; k < stalker.get(nbrs.get(j)).nbrNode.size(); k++) {
//                    if (stalker.get(nbrs.get(j)).nbrNode.get(k)==i) chk = true;
//                }
//                if(!chk) {}//System.out.println("BUG IN NEIGHBOURHOOD");

                float x1 = stalker.get(i).loc[0]+.5f;
                float y1 = stalker.get(i).loc[1]+.5f;
                float x2 = stalker.get(nbrs.get(j)).loc[0]+.5f;
                float y2 = stalker.get(nbrs.get(j)).loc[1]+.5f;

                Line l = new Line(x1, y1, x2, y2);
                l.setStrokeColor(Color.YELLOW);
                ov.add(l);

                float vnorm = (float) Math.sqrt(Math.pow(x2-x1,2)+Math.pow(y2-y1,2));

                float x3 = x2 + 0.25f * vnorm * (x1-x2)/vnorm + .05f * vnorm * (y1-y2)/vnorm;
                float y3 = y2 + 0.25f * vnorm * (y1-y2)/vnorm - .05f * vnorm * (x1-x2)/vnorm;

                float x4 = x2 + 0.25f * vnorm * (x1-x2)/vnorm - .05f * vnorm * (y1-y2)/vnorm;
                float y4 = y2 + 0.25f * vnorm * (y1-y2)/vnorm + .05f * vnorm * (x1-x2)/vnorm;

                l = new Line(x1, y1, x2, y2);
                l.setStrokeColor(Color.YELLOW);
                ov.add(l);

                l = new Line(x3, y3, x2, y2);
                l.setStrokeColor(Color.YELLOW);
                ov.add(l);

                l = new Line(x4, y4, x2, y2);
                l.setStrokeColor(Color.YELLOW);
                ov.add(l);


            }

        }

        for (int j = 0; j < stalker_gps.size(); j++) ov.add(stalker_gps.get(j)); //

        return ov;

    }

}

class BfsQueue<E> {
    private LinkedList<E> list = new LinkedList<E>();
    public void enqueue(E item) {
        list.addLast(item);
    }
    public E dequeue() {
        return list.poll();
    }
    public boolean hasItems() {
        return !list.isEmpty();
    }
    public int size() {
        return list.size();
    }
//    public void addItems(BfsQueue<? extends E> q) {
//        while (q.hasItems())
//            list.addLast(q.dequeue());
//    }
}