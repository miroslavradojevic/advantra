package aux;

import java.io.*;
import java.util.ArrayList;

/**
 * Created by miroslav on 26-7-14.
 *
 * read .det  detection file, assumes it is csv with # marking comments
 * line format (2d detection):
 * x, y, r, score, type, outwards_directions
 *
 */
public class ReadDET {

	public ArrayList<Float> 		x = new ArrayList<Float>(); // x
	public ArrayList<Float> 		y = new ArrayList<Float>(); // y
	public ArrayList<Float> 		r = new ArrayList<Float>(); // radius
	public ArrayList<Float> 		s = new ArrayList<Float>(); // score
	public ArrayList<String> 		t = new ArrayList<String>(); // type
	public ArrayList<float[][]> 	v = new ArrayList<float[][]>(); // outwards direcitons

	// indexes
	public static int XCOORD 	= 0;
	public static int YCOORD 	= 1;
	public static int RADIUS 	= 2;
	public static int SCORE		= 3;
	public static int TYPE 		= 4;

	public ReadDET(String _detFilePath) {

		String detFilePath = new File(_detFilePath).getAbsolutePath();// path to det file

		if (!(new File(detFilePath).exists())) {
			System.err.println(detFilePath+" does not exist! class not initialized...");
			return;
		}

		try {

			FileInputStream fstream 	= new FileInputStream(detFilePath);
			BufferedReader br 			= new BufferedReader(new InputStreamReader(new DataInputStream(fstream)));
			String read_line;

			while ( (read_line = br.readLine()) != null ) {
				if (!read_line.trim().startsWith("#")) { // # are comments

					String[] 	readLn = 	read_line.trim().split(",\\s*"); // comma separated

					if (readLn.length==7 || readLn.length==9 || readLn.length==11 || readLn.length==13) { // legal

						x.add(Float.valueOf(readLn[XCOORD].trim()).floatValue());
						y.add(  Float.valueOf(readLn[YCOORD].trim()).floatValue()   );
						r.add(  Float.valueOf(readLn[RADIUS].trim()).floatValue()   );
						s.add(	Float.valueOf(readLn[SCORE].trim()).floatValue()    );
						t.add(	readLn[TYPE].trim()                                 );

						float[][] vxy=null;

						if (readLn.length==7) {
							vxy = new float[1][2];
							vxy[0][0] = Float.valueOf(readLn[5].trim()).floatValue();
							vxy[0][1] = Float.valueOf(readLn[6].trim()).floatValue();
						}
						else if (readLn.length==9) {
							vxy = new float[2][2];
							vxy[0][0] = Float.valueOf(readLn[5].trim()).floatValue();
							vxy[0][1] = Float.valueOf(readLn[6].trim()).floatValue();
							vxy[1][0] = Float.valueOf(readLn[7].trim()).floatValue();
							vxy[1][1] = Float.valueOf(readLn[8].trim()).floatValue();
						}
						else if (readLn.length==11) {
							vxy = new float[3][2];
							vxy[0][0] = Float.valueOf(readLn[5].trim()).floatValue();
							vxy[0][1] = Float.valueOf(readLn[6].trim()).floatValue();
							vxy[1][0] = Float.valueOf(readLn[7].trim()).floatValue();
							vxy[1][1] = Float.valueOf(readLn[8].trim()).floatValue();
							vxy[2][0] = Float.valueOf(readLn[9].trim()).floatValue();
							vxy[2][1] = Float.valueOf(readLn[10].trim()).floatValue();
						}
						else if (readLn.length==13) {
							vxy = new float[4][2];
							vxy[0][0] = Float.valueOf(readLn[5].trim()).floatValue();
							vxy[0][1] = Float.valueOf(readLn[6].trim()).floatValue();
							vxy[1][0] = Float.valueOf(readLn[7].trim()).floatValue();
							vxy[1][1] = Float.valueOf(readLn[8].trim()).floatValue();
							vxy[2][0] = Float.valueOf(readLn[9].trim()).floatValue();
							vxy[2][1] = Float.valueOf(readLn[10].trim()).floatValue();
							vxy[3][0] = Float.valueOf(readLn[11].trim()).floatValue();
							vxy[3][1] = Float.valueOf(readLn[12].trim()).floatValue();
						}

						v.add(vxy);

					}

				}
			}

			br.close();
			fstream.close();

		}
		catch (Exception e) {
			System.err.println(e.getMessage());
		}

//		System.out.println(x.size() + " detections found!");

	}

	public void print() {

		System.out.println(x.size() + " detections found : ");

		for (int i = 0; i < x.size(); i++) {
			System.out.println(">> " + (i+1) + "  " + t.get(i) + " :\t " + x.get(i) + " , " + y.get(i) + " (" + r.get(i) + ")");
		}

	}

}
