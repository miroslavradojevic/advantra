import ij.IJ;
import ij.ImagePlus;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

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

		int[] t = new int[]{3, 4};
		double[] s = new double[]{1.5, 2.0, 2.5};

		f = new ArrayList<Feat>();
		for (int i=0; i<t.length; i++) {
			for (int j=0; j<s.length; j++) {
				f.add(new Feat(t[i], s[j]));
			}
		}

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

		String featScore = "";

		for (int x=0; x<inimg.getWidth(); x++) {
			for (int y=0; y<inimg.getHeight(); y++) {
				if (inmask.getProcessor().getf(x, y)==255) {

					int cnt = 0;
					String pScore = "";
					for (int fIdx = 0 ; fIdx < f.size(); fIdx++) {
						String fScore = f.get(fIdx).extractFeatures(x, y, (FloatProcessor) inimg.getProcessor());
						if (fScore!="") {
							cnt++;
							pScore+= fScore;
						}
					}

					if (cnt>=2) {
						featScore+=pScore+"\n";
					}

				}

			}
		}

		IJ.showMessage("extracted!");

		IJ.showMessage(featScore);

	}
}
