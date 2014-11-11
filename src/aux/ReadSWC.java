package aux;

import detection2d.CritpointRegion;
import ij.IJ;

import java.io.*;
import java.util.ArrayList;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 10/23/13
 * Time: 5:02 PM
 */
public class ReadSWC {

	// list with nodes
	public ArrayList<float[]> nodes = new ArrayList<float[]>(); // 1x7 rows (swc format)

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

    public ReadSWC(String _swcFilePath) {

        String swcFilePath = new File(_swcFilePath).getAbsolutePath();// path to swc file

        if (!(new File(swcFilePath).exists())) {
            System.err.println(swcFilePath+" does not exist! class not initialized...");
            return;
        }

//        int fileLength = 0; // counter

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

					float[] 	valsLn = 	new float[7]; // x, y, z, mother_index

                    valsLn[0] = Float.valueOf(readLn[ID].trim()).floatValue();  // id
                    valsLn[1] = Float.valueOf(readLn[TYPE].trim()).floatValue();  // type

					valsLn[2] = Float.valueOf(readLn[XCOORD].trim()).floatValue();  // x, y, z
					valsLn[3] = Float.valueOf(readLn[YCOORD].trim()).floatValue();
					valsLn[4] = Float.valueOf(readLn[ZCOORD].trim()).floatValue();

					valsLn[5] = Float.valueOf(readLn[RADIUS].trim()).floatValue();  // radius

					valsLn[6] = Float.valueOf(readLn[MOTHER].trim()).floatValue();  // mother idx

					nodes.add(valsLn);

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

    }

	public void print(){

		//		for (int ii=0; ii<nodes.size(); ii++) {
		//            for (int jj=0; jj<7; jj++) {
		//                System.out.print(nodes.get(ii)[jj]+"\t");
		//            }
		//            System.out.println();
		//		}

	}

    /*
        export SWC critical points (junctions and end-points) in DET (2D) format (same as Critpoint2D export)
        or SWC with the critical points only,
        scores will be set to 1
     */
    public void exportCritpoint(
            String                  _exportPath,
            ExportCritpointFormat   _export_format,
            ExportCritpointType     _export_type
    )
    {

        float SCALE = 2; // defines how much it will scale current radius read from swc at that node

        PrintWriter logWriter = null;

        try {
            logWriter = new PrintWriter(_exportPath); logWriter.print("" +
                    "# ADVANTRA: exportCritpoint()    \n" +
                    "# format: " +          _export_format  + "\n" +
                    "# critpoint type: " +  _export_type  +   "\n" +
                    "# author: miroslavr\n");
            logWriter.close();
        } catch (FileNotFoundException ex) {}

        try {
            logWriter = new PrintWriter(new BufferedWriter(new FileWriter(_exportPath, true)));
        } catch (IOException e) {}

        // will take loaded swc file and create a new one that extracts the critical points
        int         currId, count;
        boolean     isBif,  isEnd;

        // extraction, loop the nodes of swc reconstruction
        for (int ii=0; ii<nodes.size(); ii++) {

            currId = Math.round(nodes.get(ii)[ID]);

            if (Math.round(nodes.get(ii)[MOTHER])==-1) continue; // root point is not either of those (can be endpoint potentially though)

            /*
                check if it is END (algorithm: no one referred to it)
             */

            isEnd = true;
            for (int jj=ii+1; jj<nodes.size(); jj++) { // loop onwards the rest of the nodes
                if (currId==Math.round(nodes.get(jj)[MOTHER])) {
                    isEnd = false;
                    break; // stop looping further
                }
            }

            if (isEnd) {
                // double check the first part of the list to see if there was some that was earlier and that referred to this one
                // the one that was endpoint will need to check all of the remaining nodes
                for (int jjj = 0; jjj < ii; jjj++) {
                    if (currId==Math.round(nodes.get(jjj)[MOTHER])) {
                        isEnd = false;
                        break;
                    }
                }
            }

            if (isEnd && (_export_type==ExportCritpointType.ALL || _export_type==ExportCritpointType.END)) { // only in case it was end

                // add ii node to the output det (check output format description in Detector2D.saveDetection())
                switch (_export_format) {
                    case DET_2D:
                        logWriter.println(
                                IJ.d2s(nodes.get(ii)[XCOORD], 2)+", "+IJ.d2s(nodes.get(ii)[YCOORD],2)+", "+
                                        IJ.d2s(SCALE*nodes.get(ii)[RADIUS],2)+", "+
                                        IJ.d2s(1.0,2)+", "+ CritpointRegion.RegionType.END+", "+
                                        IJ.d2s(Float.NaN)+", "+IJ.d2s(Float.NaN)
                        );
                        break;
                    case SWC:
                        logWriter.println(String.format("%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f -1",  // -1 will be the mother index !!!!
                                currId,
                                6,  // type = 6 (yellow in vaa3d)
                                nodes.get(ii)[XCOORD],
                                nodes.get(ii)[YCOORD],
                                nodes.get(ii)[ZCOORD],
                                SCALE*nodes.get(ii)[RADIUS]));
                        break;
                    default:
                        break;
                }

            }

            if (isEnd) continue; // no need to check if it is junction then otherwise continue

            /*
                check if it is BIF (JUN) (there were 2+ that referred to it)
             */

            count = 0;
            isBif = false;
            for (int jj=ii+1; jj<nodes.size(); jj++) { // check the second half of the node list
                if (currId==Math.round(nodes.get(jj)[MOTHER])) {
                    count++;
                    if (count==2) {
                        isBif = true;
                        break;
                    }
                }
            }

            // consider the first half of the node list as well - can happen that there is a continuation with a lower index referring to this one
            if (!isBif) {
                // here is a computational issue: the one that is not the junction will surely be checked all along and tak a lot of time to confirm that it was not a junction
                for (int jjj = 0; jjj < ii; jjj++) {
                    if (currId==Math.round(nodes.get(jjj)[MOTHER])) {
                        count++;
                        if (count==2) {
                            isBif = true;
                            break;
                        }
                    }
                }
            }

            if (isBif && (_export_type==ExportCritpointType.ALL || _export_type==ExportCritpointType.JUN)) {

                // add ii node to the output det (check output format description in Detector2D.saveDetection())
                switch (_export_format) {
                    case DET_2D:
                        logWriter.println( // add ii node to the output det - since swc does not contain directional info - it will be added as NaN
                                IJ.d2s(nodes.get(ii)[XCOORD], 2)+", "+IJ.d2s(nodes.get(ii)[YCOORD],2)+", "+
                                        IJ.d2s(SCALE*nodes.get(ii)[RADIUS],2)+", "+
                                        IJ.d2s(1.0,2)+", "+ CritpointRegion.RegionType.BIF+", "+
                                        IJ.d2s(Float.NaN)+", "+IJ.d2s(Float.NaN)+", "+
                                        IJ.d2s(Float.NaN)+", "+IJ.d2s(Float.NaN)+", "+
                                        IJ.d2s(Float.NaN)+", "+IJ.d2s(Float.NaN)
                        );
                        break;
                    case SWC:
                        logWriter.println(String.format("%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f -1", // add ii node to the output swc
                                currId,
                                2,  // type = 2 (red in vaa3d)
                                nodes.get(ii)[XCOORD],
                                nodes.get(ii)[YCOORD],
                                nodes.get(ii)[ZCOORD],
                                SCALE*nodes.get(ii)[RADIUS]));
                        break;
                    default:
                        break;
                }

            }
        }

        logWriter.close();

        System.out.println(_exportPath + " exported in "+_export_format+" format.");

    }

    /*
        include shift into any of the components
     */
    public void modifySwc(
            String _exportPath,
            float dx,
            float dy,
            float dz,
            float dr
    )
    {

        PrintWriter logWriter = null;

        try {
            logWriter = new PrintWriter(_exportPath); logWriter.print("# ADVANTRA: modifySwc()\n# author: miroslavr\n"); logWriter.close();
        } catch (FileNotFoundException ex) {}

        try {
            logWriter = new PrintWriter(new BufferedWriter(new FileWriter(_exportPath, true)));
        } catch (IOException e) {}

        // extraction, loop the nodes of swc reconstruction
        for (int ii=0; ii<nodes.size(); ii++) {
            logWriter.println(String.format("%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d",  // -1 will be the mother index !!!!
                            Math.round(nodes.get(ii)[ID]),
                            Math.round(nodes.get(ii)[TYPE]),
                            nodes.get(ii)[XCOORD]   + dx,
                            nodes.get(ii)[YCOORD]   + dy,
                            nodes.get(ii)[ZCOORD]   + dz,
                            nodes.get(ii)[RADIUS]   + dr,
                            Math.round(nodes.get(ii)[MOTHER])
                    )
                    );

        }

        logWriter.close();

        System.out.println(_exportPath + " exported.");

    }

    public enum ExportCritpointFormat {SWC, DET_2D}

    public enum ExportCritpointType {ALL, JUN, END}

}
