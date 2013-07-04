import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.*;
import java.io.File;
import java.util.ArrayList;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 6/22/13
 * Time: 12:21 PM
 */
public class BifDetect implements PlugInFilter {

	ImagePlus 	inimg;
	String      inimgPath;
	ImagePlus   inmask;
    ArrayList<Feat> f;

    // variables used for calculating the score
    double[] sumCluster;
    int[][] map;// = new int[3][];

	public void run(ImageProcessor imageProcessor) {

//		inmaskPath = MyOpener.open("Open mask file", false);
//		if (inmaskPath==null) {
//			IJ.log(inmaskPath);
//			return;
//		}

        f= new ArrayList<Feat>();
        f.add(new Feat(3, 1.0));
        f.add(new Feat(3, 2.0));
        f.add(new Feat(3, 3.0));
        f.add(new Feat(3, 4.0));
//        f.add(new Feat(4, 2.5));
//        f.add(new Feat(4, 3.5));

        // for score calculation
        map = new int[3][f.size()];
        for (int i=0; i<map.length; i++) {
            for (int j=0; j<map[0].length; j++) {
                map[i][j] = -1;
            }
        }
        sumCluster = new double[3];

		int     t       		= 4;
		double  scale   		= 2.0;
		String 	inmaskPath		= Tools.removeExtension(inimgPath)+".mask";
		boolean useMask 		= true;

		t     	= (int) Prefs.get("advantra.critpoint.neuron_diam", t);
		scale   = Prefs.get("advantra.critpoint.scale", scale);

		GenericDialog gd = new GenericDialog("Fit Features");
		gd.addMessage("feature parameters");
		gd.addNumericField("neuron diameter min", t, 0, 5, "pix");
		gd.addNumericField("n'hood", scale, 1, 5, "x diameter");

        gd.addMessage("detection parameters - configuration file");

        gd.addMessage("mask (avoid processing all)");
		gd.addStringField("mask path", inmaskPath, 50);
		gd.addCheckbox("", useMask);

		gd.showDialog();
		if (gd.wasCanceled()) return;

		t 		=  	(int)gd.getNextNumber();
		Prefs.set("advantra.critpoint.neuron_diam", 	t);
		scale   =   gd.getNextNumber();
		Prefs.set("advantra.critpoint.scale", 	scale);
//		D   	=   gd.getNextNumber();
//		E   	=   gd.getNextNumber();
		inmaskPath = gd.getNextString();
		useMask = gd.getNextBoolean();

//		f= new Feat(t, scale);

		inimg.show();
        inimg.getCanvas().zoomIn(0,0);
        inimg.getCanvas().zoomIn(0,0);
        inimg.getCanvas().zoomIn(0,0);
        inimg.getCanvas().zoomIn(0,0);


		if (new File(inmaskPath).exists() && useMask) {
			inmask = new ImagePlus(inmaskPath);
			inmask.setTitle("inmask");
		}
		else {
			byte[] array = new byte[inimg.getWidth() * inimg.getHeight()];
			for (int a = 0; a < array.length; a++) array[a] = (byte) 255;
			inmask = new ImagePlus("inmask", new ByteProcessor(inimg.getWidth(), inimg.getHeight(), array));
		}

//		inmask.show();
//		IJ.selectWindow("inimg");
//		IJ.run("Add Image...", "image=inmask x="+0+" y="+0+" opacity=30");
//		inmask.close();

		IJ.log("filtering...");
		long t1 = System.currentTimeMillis();
		ImagePlus score = new ImagePlus("detection", filterWithMask(f, inmask.getProcessor(), (FloatProcessor) inimg.getProcessor()));
        long t2 = System.currentTimeMillis();
        score.show();
		IJ.log("done. "+((t2-t1)/1000f)+" sec.");

        IJ.log("detection");
        t1 = System.currentTimeMillis();
        MeanShift2D ms2d = new MeanShift2D((FloatProcessor) score.getProcessor(), 5);
        ms2d.run(150, 0.0001);
        t2 = System.currentTimeMillis();
        IJ.log("done. "+((t2-t1)/1000f)+" sec.");

        // show
        Overlay ov = new Overlay();
        for (int i = 0; i < ms2d.S.length; i++) {
            PointRoi p = new PointRoi(ms2d.T[i][1]+0.5, ms2d.T[i][0]+0.5);
            p.setStrokeColor(Color.RED);
            ov.add(p);
        }
        inimg.setOverlay(ov);
        score.setOverlay(ov);

        // clustering
        //int[][] clust = ms2d.extractClust(5);
        double[][] clust1 = ms2d.extractConvPoints(1.0, 20);
        ov = new Overlay();
        for (int i = 0; i < clust1.length; i++) {
            ov.add(new PointRoi(clust1[i][1]+0.5, clust1[i][0]+0.5));
        }
        inimg.setOverlay(ov);

	}

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

    public int mapTracker(int[][] map, int i, int q)
    {

        int idx = -1;

        if (q>=0 && q<map[0].length) {

            idx = map[i][q];

            while (q>=0) {

                q--;

                if (idx!=-1) {

                    idx = map[idx][q];

                }

            }

            return idx;

        }
        else
            return idx;

    }

	public FloatProcessor filterWithMask(ArrayList<Feat> feats, ImageProcessor msk, FloatProcessor input)
	{

		FloatProcessor ipOut = new FloatProcessor(input.getWidth(), input.getHeight());

		for (int x=0; x<input.getWidth(); x++) {
			for (int y=0; y<input.getHeight(); y++) {
				if (msk.getf(x, y)==255) {

                    // geom. mean of the sum values from each cluster!!!
                    double      sum = 0;
                    int         cnt = 0;
                    boolean     first = true;

                    for (int q=0; q<feats.size(); q++) {

                        feats.get(q).getAngles(x, y, input, false);

                        // feats.get(q).sum[0],  feats.get(q).lp[0][0,1]

                        if (feats.get(q).sum!=null) { // there was something extracted

                            if (first) {   // no need to match them

								map[0][q] = 0;
								map[1][q] = 1;
								map[2][q] = 2;

                                sumCluster[0] = feats.get(q).sum[0];
                                sumCluster[1] = feats.get(q).sum[1];
                                sumCluster[2] = feats.get(q).sum[2];

								first = false;

                            }
                            else { // match using Hungarian algorithm

								//int[] conn = Tools.match(feats.get(q-1).lp, feats.get(q).lp);

								boolean[][] chkd = new boolean[3][3];
								double[][] 	dst2 = new double[3][3];

								for (int i=0; i<3; i++) {
									for (int j=0; j<3; j++) {
										dst2[i][j] = Math.pow(feats.get(q).lp[i][0]-feats.get(q).lp[j][0], 2) + Math.pow(feats.get(q).lp[i][1]-feats.get(q).lp[j][1], 2);
									}
								}

								int[] matches = new int[3];

								for (int check=0; check<3; check++) {

									double dst2Min = Double.MAX_VALUE;
									int imin = -1;
									int jmin = -1;

									for (int i=0; i<3; i++) {

										for (int j=0; j<3; j++) {
											if (!chkd[i][j] && dst2[i][j]<dst2Min) {
												dst2Min = dst2[i][j];
												imin = i;
												jmin = j;
											}
										}

										// row imin in chkd to true
										for (int w=0; w<3; w++) chkd[imin][w] = true;
										// col jmin in chkd to true
										for (int w=0; w<3; w++) chkd[w][jmin] = true;

										// imin-jmin pair
										map[imin][q] = jmin;

									}

								}

								// formed matches

//								map[0][q] = 0;
//								map[1][q] = 1;
//								map[2][q] = 2;

								// finish here ...

								sumCluster[0] = feats.get(q).sum[0];
								sumCluster[1] = feats.get(q).sum[1];
								sumCluster[2] = feats.get(q).sum[2];



                            }

                            sum += Math.log(feats.get(q).centralAvg(x, y, input));
                            sum += Math.log(feats.get(q).sum[0]);
                            sum += Math.log(feats.get(q).sum[1]);
                            sum += Math.log(feats.get(q).sum[2]);
                            cnt ++;

                        }

                    }

                    sum = (cnt>0)? Math.exp(sum/cnt) : 0 ;

                    ipOut.setf(x, y, (float) sum);

                    // just profile entropy
                    //ipOut.setf(x, y, feats.get(0).entropy);

				}
			}
		}
		return ipOut;

	}

}
