import java.io.*;
import java.util.*;

/**
 * reader for swc reconstruction file
 *
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 10/23/13
 * Time: 5:02 PM
 */
public class ReadSWC {

    // IDs have to be unique in swc format, therefore, IDs cannot repeat
    // depending on the export tool - it can happen that the IDs repeat or missing
    // input is read as the list of nodes (Node class) taking into account the linking
    // Node has the geometry + the link towards the neighbour

//    Standardized swc files (www.neuromorpho.org) -
//    0 - undefined
//    1 - soma
//    2 - axon
//    3 - (basal) dendrite
//    4 - apical dendrite
//    5+ - custom


    // new list with nodes
    public ArrayList<Node> nnodes = new ArrayList<Node>();
    public ArrayList<float[]> nodes  = new ArrayList<float[]>();

    public int maxID = Integer.MIN_VALUE;

	public float minR = Float.POSITIVE_INFINITY, maxR = Float.NEGATIVE_INFINITY;
	public float minX = Float.POSITIVE_INFINITY, maxX = Float.NEGATIVE_INFINITY;
	public float minY = Float.POSITIVE_INFINITY, maxY = Float.NEGATIVE_INFINITY;
	public float minZ = Float.POSITIVE_INFINITY, maxZ = Float.NEGATIVE_INFINITY;

	// indexes
	public static int ID 		= 0;
	public static int TYPE 		= 1;
	public static int XCOORD 	= 2;
	public static int YCOORD 	= 3;
	public static int ZCOORD 	= 4;
	public static int RADIUS 	= 5;
	public static int MOTHER 	= 6;

    public static int SWC_LINE_LENGTH = 7;

    // variables that count connections
//    int[] cnt_conn; // count connections
//    int[] cnt_stps; // count steps (only for terminal points)

    public ReadSWC(String _swcFilePath) {

        String swcFilePath = new File(_swcFilePath).getAbsolutePath();// path to swc file

        if (!(new File(swcFilePath).exists())) {
            System.err.println(swcFilePath+" does not exist!");
            return;
        }

        System.out.println("reading...\t"+_swcFilePath);

        ArrayList<float[]> nodes_load = new ArrayList<float[]>(); // 1x7 rows (swc format)

        // read the node list of line elements - it's not guaranteed that the node IDs will be compatible with indexing
        // in sense that they are arranged ascending and even that all will exist alltogether
        // also can happen (although illegal to have doubling of the IDs)
        // read it first all to see what th efull range if IDs would be to know how much to allocate
        try { // scan the file

			FileInputStream fstream 	= new FileInputStream(swcFilePath);
            BufferedReader br 			= new BufferedReader(new InputStreamReader(new DataInputStream(fstream)));
            String read_line;

			while ( (read_line = br.readLine()) != null ) { // it will break on the empty line !!!

//              System.out.println("happened that it was empty ["+read_line+"]");//+br.readLine());//+"----->"+br.readLine()+"----"+(read_line = br.readLine()) != null);

                if (read_line.isEmpty()) continue;

                if(!read_line.trim().startsWith("#")) { // # are comments

//                    fileLength++;
                    // split values

                    String[] 	readLn = 	read_line.trim().replaceAll("," , ".").split("\\s+");

                    if (readLn.length!=SWC_LINE_LENGTH) continue; // skip the line that did not have enough values

					float[] 	valsLn = 	new float[SWC_LINE_LENGTH]; // x, y, z, mother_index

                    valsLn[0] = Integer.valueOf(readLn[ID].trim()).intValue();      // id
                    valsLn[1] = Integer.valueOf(readLn[TYPE].trim()).intValue();    // type

					valsLn[2] = Float.valueOf(readLn[XCOORD].trim()).floatValue();  // x, y, z
					valsLn[3] = Float.valueOf(readLn[YCOORD].trim()).floatValue();
					valsLn[4] = Float.valueOf(readLn[ZCOORD].trim()).floatValue();

					valsLn[5] = Float.valueOf(readLn[RADIUS].trim()).floatValue();  // radius

					valsLn[6] = Integer.valueOf(readLn[MOTHER].trim()).intValue();  // mother idx

					nodes_load.add(valsLn);

                    maxID = (valsLn[0]>maxID)? (int) valsLn[0] : maxID;

					minR = (valsLn[RADIUS]<minR)? valsLn[RADIUS] : minR;
					maxR = (valsLn[RADIUS]>maxR)? valsLn[RADIUS] : maxR;

					minX = (valsLn[XCOORD]<minX)? valsLn[XCOORD] : minX;
					maxX = (valsLn[XCOORD]>maxX)? valsLn[XCOORD] : maxX;

					minY = (valsLn[YCOORD]<minY)? valsLn[YCOORD] : minY;
					maxY = (valsLn[YCOORD]>maxY)? valsLn[YCOORD] : maxY;

					minZ = (valsLn[ZCOORD]<minZ)? valsLn[ZCOORD] : minZ;
					maxZ = (valsLn[ZCOORD]>maxZ)? valsLn[ZCOORD] : maxZ;

                }

            }

			br.close();
            fstream.close();

		}
        catch (Exception e) {
            System.err.println("Error: " + e.getMessage());
        }

        // analyze
        System.out.println(nodes_load.size() + " lines");
        System.out.println("maxID= " + maxID);

        // initialize all with null, ID will correspond to the index in nnodes list
        nnodes = new ArrayList<Node>(maxID + 1); // zeroth is dummy one (so that indexing can work with IDs) hence there is one more allocated
        for (int i = 0; i <= maxID; i++) nnodes.add(i, null);

        for (int i = 0; i < nodes_load.size(); i++) { // fill the nnodes list elements only

            int     currId      = Math.round(nodes_load.get(i)[ID]);
            int     currType    = Math.round(nodes_load.get(i)[TYPE]);

            float   currX       = nodes_load.get(i)[XCOORD];
            float   currY       = nodes_load.get(i)[YCOORD];
            float   currZ       = nodes_load.get(i)[ZCOORD];

            float   currR       = nodes_load.get(i)[RADIUS];

            int     prevId      = Math.round(nodes_load.get(i)[MOTHER]);

            if (nnodes.get(currId)==null) { // add the node
                nnodes.set(currId, new Node(currX, currY, currZ, currR, currType));
//                nnodes.get(i).nbrNode.add(prevId);
            }
            else {
                System.out.println("DOUBLING NODE ID HAPPENED!");
                // something is there, this means that we're doubling the same node (illegal but can happen)
                // by convention we'll assume that the doubled node is the same x,y,z,r - so we won't read the new coordinates again but the neighbours only
//                nnodes.get(i).nbrNode.add(prevId);
            }

        }

        // once the nodes are added add 2-directional connections read from the swc file (both neighbour references can be added as the nodes exist)
        for (int i = 0; i < nodes_load.size(); i++) {

            int     currId      = Math.round(nodes_load.get(i)[ID]);
            int     prevId      = Math.round(nodes_load.get(i)[MOTHER]);

//            System.out.println(currId + " -- " + prevId + " --- " + nnodes.size());
//            if (nnodes.get(currId)==null) System.out.println("it was null");

            if (prevId!=-1) {
                nnodes.get(currId).nbrNode.add(prevId);
                nnodes.get(prevId).nbrNode.add(currId);
            }

        }

        // remove duplicate neighbouring links so that it does not confuse the bfs when checking connectivity
        remove_duplicate_neighbourhoods(nnodes);

        // calculate the connectivity using bfs expanding from each individual node giving the amount of nodes
        // necessary to reach the other node and set the nr. steps that was necessary to pass
        // 
        int[][]     nnodes_connt = new int[nnodes.size()][nnodes.size()];
        get_connectivity(nnodes, nnodes_connt);

        float[][]   nnodes_ovlps = new float[nnodes.size()][nnodes.size()];
        get_overlapping(nnodes, nnodes_ovlps);

        // link those node pairs that overlap but are not connected
        ArrayList<int[]> linkings = new ArrayList<int[]>();
        for (int i = 0; i < nnodes.size(); i++) {
//            System.out.println(Arrays.toString(nnodes_ovlps[i]));
            for (int j = i+1; j < nnodes.size(); j++) {

                // overlapping but not connected
                if (nnodes_ovlps[i][j]>Float.MIN_VALUE && nnodes_connt[i][j]==-1 && nnodes_connt[j][i]==-1) { //
                    int[] link = new int[]{i, j};
                    linkings.add(link);
                    System.out.println(Arrays.toString(link));
                }
                // overlapping and connected in N steps

            }
        }

        // link those that overlap but are connected with large amount of steps



        // remove duplicate connections


        // extract the tree of nodes (list of nodes) by travesrsing the whole tree once the connections have been made



//        float min_overlap = Float.MAX_VALUE;
//        boolean[][] checked = new boolean[nnodes.size()][nnodes.size()];
//        boolean found_overlap = false;
//        int i_overlap, j_overlap;
//        int cnt = 0;
//        do {
//
//            min_overlap = Float.MAX_VALUE;
//            i_overlap = -1;
//            j_overlap = -1;
//            found_overlap = false;
//
//            for (int i = 0; i < nnodes.size(); i++) {
//                for (int j = 0; j < nnodes.size(); j++) {
//                    if (nnodes.get(i)!=null && nnodes.get(j)!=null) {
//
//                        float oo = nnodes.get(i).overlap(nnodes.get(j));
//
//                        if ( j!=i && !nnodes.get(i).nbrNode.contains(j) && oo>0 && !checked[i][j] ) {
//                            min_overlap = oo;
//                            i_overlap = i;
//                            j_overlap = j;
//                            found_overlap = true;
//                        }
//
//                    }
//                }
//            }
//
//            if (found_overlap) {
//
//                checked[i_overlap][j_overlap] = true;
//
//                ArrayList<Integer> itersec = Node.intersection(nnodes.get(i_overlap).nbrNode, nnodes.get(j_overlap).nbrNode);
//
//                if (itersec.size()==0) { // there is no connection with the neighbours of i_overlap
//                    nnodes.get(i_overlap).nbrNode.add(j_overlap); // add bidirectional neighbourhood
//                    nnodes.get(j_overlap).nbrNode.add(i_overlap);
//                }
//                else {
//
//                    // find that node, neighbour of i_overlap that is also the neighbour of j_overlap
//                    Node nbridge = nnodes.get(itersec.get(0));
//
//                    if (min_overlap > nnodes.get(i_overlap).overlap(nbridge)) {
//                        nnodes.get(i_overlap).nbrNode.add(j_overlap); // add bidirectional neighbourhood
//                        nnodes.get(j_overlap).nbrNode.add(i_overlap);
//                    }
//
//                }
//
//            }
//
//            cnt++;
//
//            System.out.println("cnt=" + cnt);
//
//        }
//        while (found_overlap && cnt<=10);


        // loop through the whole sequence to find the terminal points and junctions
//        cnt_conn = new int[nodes.size()]; // store the number of connections for that node
//        cnt_stps = new int[nodes.size()]; // count the number of steps towards the root critical point (only for critical points)

//        Arrays.fill(cnt_conn, 0);
//        Arrays.fill(cnt_stps, -1);

//        for (int i = 0; i < nodes.size(); i++) { // loop all the nodes to fill the table
//
//            int id 		= Math.round(nodes.get(i)[ReadSWC.ID]);
//
//            if (id<=0) {
//                System.out.println("found node id <= 0   " + id);
//                return;
//            }
//
//            int id_mother 	= Math.round(nodes.get(i)[ReadSWC.MOTHER]);
//
//            // store both values in the table - table index corresponds to the ID
//            cnt_conn[id-1]++;
//            if (id_mother!=-1) {
////                System.out.println("it was---> " + id_mother);
//                cnt_conn[id_mother-1]++;
//            }
//
//        }

//        for (int i = 0; i < nodes.size(); i++) { // depending on the table value, recursively go back to the nearest critical point node
//
//            if (isEndPoint(i)) { // only count the length for endpoints
//
//                // trace back
//                int next_id = Math.round(nodes.get(i)[ReadSWC.MOTHER]);
//
//                boolean isCP = false;
//                cnt_stps[i] = 0;
//
//                while (!isCP) {
//
//                    if (next_id!=-1) {
//                        isCP = isCriticalPoint(next_id-1);
//                        next_id = Math.round(nodes.get(next_id-1)[ReadSWC.MOTHER]);
//                    }
//                    else {
//                        isCP = true; // will break the loop
//                    }
//
//                    cnt_stps[i]++;
//
//                }
//            }
//
//        }

    }

    private void get_overlapping(ArrayList<Node> _nnodes, float[][] _nnodes_ovlps) {
        // will calculate volumetric overlap between nodes

        System.out.print("\nextracting node overlap... ");

        for (int i = 0; i < _nnodes.size(); i++) {
            for (int j = i; j < _nnodes.size(); j++) {
                if (_nnodes.get(i)!=null && _nnodes.get(j)!=null) {
                    if (i==j) {
                        _nnodes_ovlps[i][j] = -2;
                    }
                    else {
                        float vol_ovlp = _nnodes.get(i).overlap(_nnodes.get(j));
                        _nnodes_ovlps[i][j] = vol_ovlp;
                        _nnodes_ovlps[j][i] = vol_ovlp;
                    }
                }
                else {
                    _nnodes_ovlps[i][j] = -2;
                    _nnodes_ovlps[j][i] = -2;
                }

            }
        }

        System.out.println("done.");

    }

    private void get_connectivity(ArrayList<Node> _nnodes, int[][] _nnodes_connt) {
        // will expand BFS from each node and assign every other node with the connectivity info
        // nr nodes to get from each other node to this one

        System.out.print("\nextracting connectivity... ");

        for (int i = 0; i < _nnodes.size(); i++) {
            if (_nnodes.get(i)!=null) {
                // initialize connectivity info: -2: the other one is null, -1: they are not overlapping
                for (int j = 0; j < _nnodes.size(); j++) _nnodes_connt[i][j] = (_nnodes.get(j) == null) ? -2 : -1;
                // extract connectivity by traversing from the current node
                bfs_connectivity(i, _nnodes, _nnodes_connt[i]); // seed, tree list, nr nodes to the seed node
            }
            else Arrays.fill(_nnodes_connt[i], -2);
        }

        System.out.println("done.");

    }

    private ArrayList[] init_discovered_list(ArrayList<Node> _nnodes) {
        // will be used to book-keep discovered traces (graph vertices)
        ArrayList[] discovered = new ArrayList[_nnodes.size()];
        for (int i = 0; i < _nnodes.size(); i++) {

            discovered[i] = new ArrayList<Boolean>(); // _nnodes.get(i).nbrNode.size()

            if (_nnodes.get(i)!=null) {
                for (int j = 0; j < _nnodes.get(i).nbrNode.size(); j++) {
                    discovered[i].add(false);
                }
            }


        }
        return discovered;
    }

    private void bfs_connectivity(int seed_node_idx, ArrayList<Node> linked_nodes, int[] connectivity){

        /**
         *  breadth-first search (BFS) to traverse the tree from seed node jotting the number of nodes to reach seed node
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

        ArrayList[] discovered = init_discovered_list(linked_nodes); // initialize discovered list, will be reset before each bfs

        BfsQueue bfsQueue = new BfsQueue();

        connectivity[seed_node_idx] = 0;

        // add the neighbors to the queue and label them as discovered
        for (int j = 0; j <linked_nodes.get(seed_node_idx).nbrNode.size(); j++) {
            int next = linked_nodes.get(seed_node_idx).nbrNode.get(j);
            // enqueue(), add to FIFO structure, http://en.wikipedia.org/wiki/Queue_%28abstract_data_type%29
            bfsQueue.enqueue(new int[]{seed_node_idx, next});
            discovered[seed_node_idx].set(j, true);                                          // set label to discovered in both neighbouting index lists
            discovered[next].set(linked_nodes.get(next).nbrNode.indexOf(seed_node_idx), true);   // index where the background link was found
        }

//        System.out.println("BFS expands from node " + seed_node_idx + " --- |Q| = " + bfsQueue.size() + " elements.");

        while (bfsQueue.hasItems()) {

            // dequeue(), take from FIFO structure, http://en.wikipedia.org/wiki/Queue_%28abstract_data_type%29
            int [] getLnk = (int[]) bfsQueue.dequeue();

            // next neighbour at the time it was added to the queue becomes current
            int prev = getLnk[0]; //Q.get(Q.size()-1)[0];
            int curr = getLnk[1]; //Q.get(Q.size()-1)[1];

            connectivity[curr] = connectivity[prev] + 1;

            while(Collections.frequency(discovered[curr], false)==1) { // step further if only one undiscovered

                prev = curr;
                curr = linked_nodes.get(curr).nbrNode.get(discovered[curr].indexOf(false));

//                System.out.println(prev + " -- " + curr);
                connectivity[curr] = connectivity[prev] + 1;

                // mark as discovered the connections curr--prev and prev--curr
                discovered[curr].set(linked_nodes.get(curr).nbrNode.indexOf(prev), true);
                discovered[prev].set(linked_nodes.get(prev).nbrNode.indexOf(curr), true);

            }

            // there is !=1 neighbour
            for (int i = 0; i < discovered[curr].size(); i++) {
                boolean isDiscovered =  (Boolean) discovered[curr].get(i);
                if (!isDiscovered) { // if it was not discovered (0,2,3...)
                    int next = linked_nodes.get(curr).nbrNode.get(i);

                    bfsQueue.enqueue(new int[]{curr, next});   // enqueue()

                    discovered[curr].set(i, true); // label as discovered
                    discovered[next].set(linked_nodes.get(next).nbrNode.indexOf(curr), true);

                }
            }

        }

    }

    private void remove_duplicate_neighbourhoods(ArrayList<Node> _nnodes) {

        // remove duplicate neighbourhood links
        System.out.print("\nremove duplicated links... ");
        for (int i = 0; i < _nnodes.size(); i++) {
            if (_nnodes.get(i)!=null) {
                Set<Integer> set = new HashSet<Integer>();
                set.addAll(_nnodes.get(i).nbrNode);
                _nnodes.get(i).nbrNode.clear();
                _nnodes.get(i).nbrNode.addAll(set);
            }
        }
        System.out.println("done.");

        // check if the neighborhoods are conistent
        System.out.print("\nchecking neighbourhood consistency... ");
        for (int i = 0; i < _nnodes.size(); i++) {
            if (_nnodes.get(i)!=null) {
                for (int j = 0; j < _nnodes.get(i).nbrNode.size(); j++) {
                    int nbr_idx = _nnodes.get(i).nbrNode.get(j);
                    if (nbr_idx>0) {
                        if (Collections.frequency(_nnodes.get(nbr_idx).nbrNode, i)!=1)
                            System.out.println("ERROR: " + i + " -- " + nbr_idx);
                    }
                    else {
//                        System.out.println("ERROR: parent id was " + nbr_idx);
                    }
                }
            }
        }
        System.out.println("done.");
    }

    public boolean isEndPoint(int node_idx)
    {
        return true; // cnt_conn[node_idx]==1;
    }

    public boolean isBifPoint(int node_idx)
    {
        return true; // cnt_conn[node_idx]==3;
    }

    public boolean isCrsPoint(int node_idx)
    {
        return true; // cnt_conn[node_idx]==4;
    }

    public boolean isCriticalPoint(int node_idx)
    {
        return true; // cnt_conn[node_idx]!=2;
    }

//    public boolean isSomaNode(int node_idx)
//    {
//        return Math.round(nodes.get(node_idx)[TYPE])==SOMA;
//    }

    public void printCritpointCounts()
    {
//        for (int i = 0; i < cnt_conn.length; i++) {
//            if (isCriticalPoint(i)) {
//                System.out.println( "id" + (i+1) + " --> " + cnt_conn[i] + " connections\t" + cnt_stps[i] + " steps");
//            }
//        }
    }

//    public int getSteps(int node_idx)
//    {
//        return 1; // cnt_stps[node_idx]; // -1, unless it is an endpoint, then it counts the path towards nearest critpoint
//    }

//    public float maxDiameter()
//    {
//
//        float mx = Float.NEGATIVE_INFINITY;
//        for (int i = 0; i < nodes.size(); i++) {
//            if (Math.round(nodes.get(i)[TYPE]) != SOMA) {
//                if (nodes.get(i)[RADIUS]>mx)
//                    mx = nodes.get(i)[RADIUS];
//            }
//        }
//
//        return mx;
//
//    }

//    public void lowRadiusBoundary(float low_bdry)
//    {
//
//        for (int i = 0; i < nodes.size(); i++) {
//            if (Math.round(nodes.get(i)[TYPE]) != SOMA) {
//                if (nodes.get(i)[RADIUS]<low_bdry)
//                    nodes.get(i)[RADIUS] = low_bdry;
//            }
//        }
//
//    }

//    public float averageDiameter()
//    {
//
//        float avg = 0;
//
//        if (nodes.size()>0) {
//
//            int cnt = 0;
//
//            for (int i = 0; i < nodes.size(); i++) {
//                if (Math.round(nodes.get(i)[TYPE]) != SOMA) {
//                    cnt++;
//                    avg += 2 * nodes.get(i)[RADIUS];
//                }
//            }
//
//            avg = avg / cnt;
//
//        }
//        else {
//            avg = Float.NaN;
//        }
//
//        return avg;
//
//    }

//    public float medianDiameter()
//    {
//
//        int count = 0;
//
//        for (int i = 0; i < nodes.size(); i++)
//            if (Math.round(nodes.get(i)[TYPE]) != SOMA) count++;
//
//        float[] vals = new float[count];
//
//        count = 0;
//        for (int i = 0; i < nodes.size(); i++) {
//            if (Math.round(nodes.get(i)[TYPE]) != SOMA) vals[count++] = 2 * nodes.get(i)[RADIUS];
//        }
//
//        return Toolbox.median(vals);
//
//    }

//    public void print(){
//
//				for (int ii=0; ii<nodes.size(); ii++) {
//		            for (int jj=0; jj<nodes.get(ii).length; jj++) {
//		                System.out.print(nodes.get(ii)[jj]+"   ");
//		            }
//		            System.out.println();
//				}
//
//	}

//    /*
//        export SWC critical points (junctions and end-points) in DET (2D) format (same as Critpoint2D export)
//        or SWC with the critical points only,
//        scores will be set to 1
//     */
//    public void exportCritpoint(
//            String                  _exportPath,
//            ExportCritpointFormat   _export_format,
//            ExportCritpointType     _export_type
//    )
//    {
//
//        float SCALE = 2; // defines how much it will scale current radius read from swc at that node
//
//        PrintWriter logWriter = null;
//
//        try {
//            logWriter = new PrintWriter(_exportPath); logWriter.print("" +
//                    "# ADVANTRA: exportCritpoint()    \n" +
//                    "# format: " +          _export_format  + "\n" +
//                    "# critpoint type: " +  _export_type  +   "\n" +
//                    "# author: miroslavr\n");
//            logWriter.close();
//        } catch (FileNotFoundException ex) {}
//
//        try {
//            logWriter = new PrintWriter(new BufferedWriter(new FileWriter(_exportPath, true)));
//        } catch (IOException e) {}
//
//        // will take loaded swc file and create a new one that extracts the critical points
//        int         currId, count;
//        boolean     isBif,  isEnd;
//
//        // extraction, loop the nodes of swc reconstruction
//        for (int ii=0; ii<nodes.size(); ii++) {
//
//            currId = Math.round(nodes.get(ii)[ID]);
//
//            if (Math.round(nodes.get(ii)[MOTHER])==-1) continue; // root point is not either of those (can be endpoint potentially though)
//
//            /*
//                check if it is END (algorithm: no one referred to it)
//             */
//
//            isEnd = true;
//            for (int jj=ii+1; jj<nodes.size(); jj++) { // loop onwards the rest of the nodes
//                if (currId==Math.round(nodes.get(jj)[MOTHER])) {
//                    isEnd = false;
//                    break; // stop looping further
//                }
//            }
//
//            if (isEnd) {
//                // double check the first part of the list to see if there was some that was earlier and that referred to this one
//                // the one that was endpoint will need to check all of the remaining nodes
//                for (int jjj = 0; jjj < ii; jjj++) {
//                    if (currId==Math.round(nodes.get(jjj)[MOTHER])) {
//                        isEnd = false;
//                        break;
//                    }
//                }
//            }
//
//            if (isEnd && (_export_type==ExportCritpointType.ALL || _export_type==ExportCritpointType.END)) { // only in case it was end
//
//                // add ii node to the output det (check output format description in Detector2D.saveDetection())
//                switch (_export_format) {
//                    case DET_2D:
//                        logWriter.println(
//                                IJ.d2s(nodes.get(ii)[XCOORD], 2)+", "+IJ.d2s(nodes.get(ii)[YCOORD],2)+", "+
//                                        IJ.d2s(SCALE*nodes.get(ii)[RADIUS],2)+", "+
//                                        IJ.d2s(1.0,2)+", "+ CritpointRegion.RegionType.END+", "+
//                                        IJ.d2s(Float.NaN)+", "+IJ.d2s(Float.NaN)
//                        );
//                        break;
//                    case SWC:
//                        logWriter.println(String.format("%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f -1",  // -1 will be the mother index !!!!
//                                currId,
//                                6,  // type = 6 (yellow in vaa3d)
//                                nodes.get(ii)[XCOORD],
//                                nodes.get(ii)[YCOORD],
//                                nodes.get(ii)[ZCOORD],
//                                SCALE*nodes.get(ii)[RADIUS]));
//                        break;
//                    default:
//                        break;
//                }
//
//            }
//
//            if (isEnd) continue; // no need to check if it is junction then otherwise continue
//
//            /*
//                check if it is BIF (JUN) (there were 2+ that referred to it)
//             */
//
//            count = 0;
//            isBif = false;
//            for (int jj=ii+1; jj<nodes.size(); jj++) { // check the second half of the node list
//                if (currId==Math.round(nodes.get(jj)[MOTHER])) {
//                    count++;
//                    if (count==2) {
//                        isBif = true;
//                        break;
//                    }
//                }
//            }
//
//            // consider the first half of the node list as well - can happen that there is a continuation with a lower index referring to this one
//            if (!isBif) {
//                // here is a computational issue: the one that is not the junction will surely be checked all along and tak a lot of time to confirm that it was not a junction
//                for (int jjj = 0; jjj < ii; jjj++) {
//                    if (currId==Math.round(nodes.get(jjj)[MOTHER])) {
//                        count++;
//                        if (count==2) {
//                            isBif = true;
//                            break;
//                        }
//                    }
//                }
//            }
//
//            if (isBif && (_export_type==ExportCritpointType.ALL || _export_type==ExportCritpointType.JUN)) {
//
//                // add ii node to the output det (check output format description in Detector2D.saveDetection())
//                switch (_export_format) {
//                    case DET_2D:
//                        logWriter.println( // add ii node to the output det - since swc does not contain directional info - it will be added as NaN
//                                IJ.d2s(nodes.get(ii)[XCOORD], 2)+", "+IJ.d2s(nodes.get(ii)[YCOORD],2)+", "+
//                                        IJ.d2s(SCALE*nodes.get(ii)[RADIUS],2)+", "+
//                                        IJ.d2s(1.0,2)+", "+ CritpointRegion.RegionType.BIF+", "+
//                                        IJ.d2s(Float.NaN)+", "+IJ.d2s(Float.NaN)+", "+
//                                        IJ.d2s(Float.NaN)+", "+IJ.d2s(Float.NaN)+", "+
//                                        IJ.d2s(Float.NaN)+", "+IJ.d2s(Float.NaN)
//                        );
//                        break;
//                    case SWC:
//                        logWriter.println(String.format("%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f -1", // add ii node to the output swc
//                                currId,
//                                2,  // type = 2 (red in vaa3d)
//                                nodes.get(ii)[XCOORD],
//                                nodes.get(ii)[YCOORD],
//                                nodes.get(ii)[ZCOORD],
//                                SCALE*nodes.get(ii)[RADIUS]));
//                        break;
//                    default:
//                        break;
//                }
//
//            }
//        }
//
//        logWriter.close();
//
//        System.out.println(_exportPath + " exported in "+_export_format+" format.");
//
//    }

    /*
        include shift into any of the components
     */
//    public void modifySwc(
//            String _exportPath,
//            float dx,
//            float dy,
//            float dz,
//            float dr
//    )
//    {
//
//        PrintWriter logWriter = null;
//
//        try {
//            logWriter = new PrintWriter(_exportPath); logWriter.print("# ADVANTRA: modifySwc()\n# author: miroslavr\n"); logWriter.close();
//        } catch (FileNotFoundException ex) {}
//
//        try {
//            logWriter = new PrintWriter(new BufferedWriter(new FileWriter(_exportPath, true)));
//        } catch (IOException e) {}
//
//        // extraction, loop the nodes of swc reconstruction
//        for (int ii=0; ii<nodes.size(); ii++) {
//            logWriter.println(String.format("%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d",  // -1 will be the mother index !!!!
//                            Math.round(nodes.get(ii)[ID]),
//                            Math.round(nodes.get(ii)[TYPE]),
//                            nodes.get(ii)[XCOORD]   + dx,
//                            nodes.get(ii)[YCOORD]   + dy,
//                            nodes.get(ii)[ZCOORD]   + dz,
//                            nodes.get(ii)[RADIUS]   + dr,
//                            Math.round(nodes.get(ii)[MOTHER])
//                    )
//                    );
//
//        }
//
//        logWriter.close();
//
//        System.out.println(_exportPath + " exported.");
//
//    }

//    public void umToPix(float _um_pix_xy, float _um_pix_z)
//    {
//
//        // will convert all the metric about neuron from um into pix using the constant saying how many um fit into one pixel
//        minX = Float.POSITIVE_INFINITY;
//        maxX = Float.NEGATIVE_INFINITY;
//
//        minY = Float.POSITIVE_INFINITY;
//        maxY = Float.NEGATIVE_INFINITY;
//
//        minZ = Float.POSITIVE_INFINITY;
//        maxZ = Float.NEGATIVE_INFINITY;
//
//        minR = Float.POSITIVE_INFINITY;
//        maxR = Float.NEGATIVE_INFINITY;
//
//        for (int i = 0; i < nodes.size(); i++) {
//
//            nodes.get(i)[XCOORD] /= _um_pix_xy;
//            nodes.get(i)[YCOORD] /= _um_pix_xy;
//            nodes.get(i)[ZCOORD] /= _um_pix_z;
//            nodes.get(i)[RADIUS] /= _um_pix_xy; // use the same one as for xy (not sure if it should be separate parameter)
//
//            minR = (nodes.get(i)[RADIUS]<minR)? nodes.get(i)[RADIUS] : minR;
//            maxR = (nodes.get(i)[RADIUS]>maxR)? nodes.get(i)[RADIUS] : maxR;
//
//            minX = (nodes.get(i)[XCOORD]<minX)? nodes.get(i)[XCOORD] : minX;
//            maxX = (nodes.get(i)[XCOORD]>maxX)? nodes.get(i)[XCOORD] : maxX;
//
//            minY = (nodes.get(i)[YCOORD]<minY)? nodes.get(i)[YCOORD] : minY;
//            maxY = (nodes.get(i)[YCOORD]>maxY)? nodes.get(i)[YCOORD] : maxY;
//
//            minZ = (nodes.get(i)[ZCOORD]<minZ)? nodes.get(i)[ZCOORD] : minZ;
//            maxZ = (nodes.get(i)[ZCOORD]>maxZ)? nodes.get(i)[ZCOORD] : maxZ;
//
//        }
//
//    }

    public enum ExportCritpointFormat {SWC, DET_2D}

    public enum ExportCritpointType {ALL, JUN, END}

}
