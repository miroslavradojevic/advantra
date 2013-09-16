import ij.IJ;
import ij.ImagePlus;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import aux.Tools;

import java.io.*;
import java.util.ArrayList;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 7/2/13
 * Time: 4:25 AM
 */
public class ExtractFeatures implements PlugInFilter {

	ImagePlus 		inimg;
	String      	inimgPath;
	ImagePlus   	inmask;
	ArrayList<Feat> f;
	String 			featFile = "feats.csv";

	public int setup(String s, ImagePlus imagePlus)
	{
		if(imagePlus==null) {
			IJ.showMessage("needs opened image"); return DONE;
		}
		inimg = Tools.convertToFloatImage(imagePlus);
		inimg.setTitle("inimg");
		inimgPath = imagePlus.getOriginalFileInfo().directory+imagePlus.getOriginalFileInfo().fileName;
		return DOES_8G+DOES_32+NO_CHANGES;
	}

	public void run(ImageProcessor imageProcessor) {

		int[] t = new int[]{4};
		double[] s = new double[]{1.5, 2.0, 3.0};

		f = new ArrayList<Feat>();
		for (int i=0; i<t.length; i++) {
			for (int j=0; j<s.length; j++) {
				f.add(new Feat(t[i], s[j]));
			}
		}

        // to empty the file
        PrintWriter writer = null;
        try {
            writer = new PrintWriter(featFile);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        writer.print("");
        writer.close();

		inimg.show();

		String inmaskPath = Tools.removeExtension(inimgPath)+".mask";

		if (new File(inmaskPath).exists()) {
			inmask = new ImagePlus(inmaskPath);
			inmask.setTitle("inmask");
		}
		else {
			byte[] array = new byte[inimg.getWidth() * inimg.getHeight()];
			for (int a = 0; a < array.length; a++) array[a] = (byte) 255;
			inmask = new ImagePlus("inmask", new ByteProcessor(inimg.getWidth(), inimg.getHeight(), array));
		}

		inmask.show();

        ArrayList<ArrayList<Double>> allLocs = new ArrayList<ArrayList<Double>>();

        long t1 = System.currentTimeMillis();

		for (int x=0; x<inimg.getWidth(); x++) {
			for (int y=0; y<inimg.getHeight(); y++) {
				if (inmask.getProcessor().getf(x, y)==255) {

                    ArrayList<Double> oneLoc = new ArrayList<Double>();

                    oneLoc.add(Double.valueOf(x));
                    oneLoc.add(Double.valueOf(y));

					int cnt = 0;
					for (int fIdx = 0 ; fIdx < f.size(); fIdx++) {
						double[] featureContainer = new double[15];
                        //f.get(fIdx).extractFeatures(x, y, (FloatProcessor) inimg.getProcessor());
						if (f.get(fIdx).ap!=null) {
							cnt++;
                            // add the features
                            for (int w=0; w<15; w++)
                                oneLoc.add(featureContainer[w]);
						}
                        else {
                            for (int w=0; w<15; w++)
                                oneLoc.add(Double.NaN);
                        }
					}

					if (cnt>=2) {    // at least 2 features scored
                        allLocs.add(oneLoc);
					}

				}

			}
		}

        long t2 = System.currentTimeMillis();

		IJ.showMessage("extracted! "+allLocs.size()+" locations "+((t2-t1)/1000f)+" sec.");

        // save it
        // append line to output file (features to be used)
        try {
            PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(featFile, true)));

            for (int i=0; i<allLocs.size(); i++) {
                for (int j=0; j<allLocs.get(i).size(); j++) {
                    out.print(""+allLocs.get(i).get(j));
                    if(j<allLocs.get(i).size()-1) {
                        out.print(", ");
                    }
                    else {
                        out.print("\n");
                    }
                }
            }
            out.close();

        } catch (IOException e1) {}

	}
}
