import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.*;
import ij.plugin.PlugIn;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 6/11/13
 * Time: 2:53 PM
 */
public class TryFeatures implements PlugInFilter, MouseListener {

    ImagePlus   inimg;
	ImagePlus 	profileImg;
    String      inimgPath;
    String      confFile;

	ArrayList<Feat> f;

    public void run(ImageProcessor imageProcessor)
	{

        int     t       = 4;
        double  scale   = 2.0;
        double  D       = 10;

/*        GenericDialog gd = new GenericDialog("Fit Features");
        gd.addMessage("feature parameters");
        gd.addNumericField("neuron diameter", t, 0, 5, "pix");
        gd.addNumericField("check range radius", scale, 1, 5, "x diameter");
        gd.showDialog();
        if (gd.wasCanceled()) return;
		t 		=  	(int)gd.getNextNumber();
        scale   =   gd.getNextNumber();*/

		f= new ArrayList<Feat>();
        f.add(new Feat(3, 1.5));
        f.add(new Feat(3, 2.0));
        f.add(new Feat(3, 3.0));
        f.add(new Feat(4, 1.5));
        f.add(new Feat(4, 2.0));
        f.add(new Feat(4, 3.0));

		//new ImagePlus("offsets",f.plotOffsets()).show();
		//new ImagePlus("pattern", f.plotTemplate(new double[]{0, (1f/2)*Math.PI, (2f/2)*Math.PI})).show();

        inimg.show();

        inimg.getCanvas().addMouseListener(this);

        confFile = "current.conf";
        PrintWriter writer = null;
        try {
            writer = new PrintWriter(confFile);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        writer.print("");
        writer.close();

		profileImg = new ImagePlus();

    }

    public void mouseClicked(MouseEvent e)
    {

        ImageCanvas srcCanv = (ImageCanvas) e.getSource();
        int atX = 	srcCanv.offScreenX(e.getX());
		int atY = 	srcCanv.offScreenY(e.getY());

		FloatProcessor inip = (FloatProcessor) inimg.getProcessor().duplicate(); // make it possible to cast here  (FloatProcessor)

        Overlay msOverlay = new Overlay();
        ImageStack isPlot;
        ImageProcessor ipPlot;

        ipPlot = f.get(0).getAngles(atX, atY, inip, true);
        if (f.get(0).ap!=null) {
            PointRoi p;
            for (int q=0; q<3; q++) {
                p = new PointRoi(
                        atX+0.5*(f.get(0).r+f.get(0).rInner)*Math.cos(f.get(0).ap[q])+0.5,
                        atY+0.5*(f.get(0).r+f.get(0).rInner)*Math.sin(f.get(0).ap[q])+0.5);
                msOverlay.add(p);
            }
        }

        isPlot = new ImageStack(ipPlot.getWidth(), ipPlot.getHeight());
        isPlot.addSlice(ipPlot);

        for (int fIdx = 1; fIdx<f.size(); fIdx++) {

            ipPlot = f.get(fIdx).getAngles(atX, atY, inip, true);
            if (f.get(fIdx).ap!=null) {
                PointRoi p;
                for (int q=0; q<3; q++) {
                    p = new PointRoi(
                            atX+0.5*(f.get(fIdx).r+f.get(fIdx).rInner)*Math.cos(f.get(fIdx).ap[q])+0.5,
                            atY+0.5*(f.get(fIdx).r+f.get(fIdx).rInner)*Math.sin(f.get(fIdx).ap[q])+0.5);
                    p.setStrokeColor(Color.BLUE);
                    msOverlay.add(p);
                }
            }
            isPlot.addSlice(ipPlot);

        }

        inimg.setOverlay(msOverlay);


		profileImg.setStack(isPlot);
		profileImg.updateAndDraw();
		profileImg.show();

//		if (f.ap!=null) {
//			IJ.log("final angles:");
//			for (int g = 0; g<f.ap.length; g++) {
//				IJ.log("angle"+g+" > "+(f.ap[g]*(360f/((float)Math.PI*2)))+" degrees");
//			}
//		}

		//f.regionScores(atX, atY, inip, f.ap, true);  // true means that it will add locations to Overlay() if there are enough directions

        //inimg.setOverlay(f.ov);

//		if (f.ap!=null) {
//			IJ.log("ap (sorted) angles:");
//			for (int g = 0; g<f.ap.length; g++) {
//				IJ.log("angle"+g+" > "+(f.ap[g]*(360f/((float)Math.PI*2)))+" degrees");
//			}
//		}



        // append line to output file (features to be used)
        try {
            PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(confFile, true)));
            out.println("nothing...");
            out.close();
        } catch (IOException e1) {
            //oh noes!
        }

/*        if (f.ap!=null && f.ap.length>=3) {

//            ImagePlus templateFit = new ImagePlus("template", f.plotTemplate(f.ap));//   //  // new double[]{0, 0.5*Math.PI, Math.PI}
//			templateFit.show();
//            IJ.selectWindow("inimg");
//            IJ.run("Add Image...", "image=template x="+(atX-templateFit.getWidth()/2)+" y="+(atY-templateFit.getHeight()/2)+" opacity=50");
//			templateFit.close();

            // calculate the moments for each patch and show
            double[] c = new double[2];
            double[] th = new double[1];
            Overlay ovPtch = new Overlay();

            for (int dir=0; dir<3; dir++) {

                Tools.extractMoments2D(f.patches3[dir], c, th);
                IJ.log("cen: "+Arrays.toString(c)+" theta: "+Arrays.toString(th));
                PointRoi p = new PointRoi(c[0]+0.5, c[1]+0.5);
                p.setStrokeColor(Color.RED);
                p.setPosition(dir+1);
                ovPtch.add(p);
                Line l = new Line(c[0]+0.5, c[1]+0.5, c[0]+0.5+5*Math.sin(th[0]),  c[1]+0.5+5*(-1)*Math.cos(th[0]));
                l.setStrokeColor(Color.YELLOW);
                l.setPosition(dir+1);
                ovPtch.add(l);

            }

            ImagePlus showC;
            showC = new ImagePlus("ptchs", f.showRegionPatches());
            showC.setOverlay(ovPtch);
            showC.show();
            for (int q=0; q<5; q++) showC.getCanvas().zoomIn(0, 0);

		}*/

		IJ.selectWindow("inimg");

    }

    public void mousePressed(MouseEvent e) {
    }

    public void mouseReleased(MouseEvent e) {
    }

    public void mouseEntered(MouseEvent e) {
    }

    public void mouseExited(MouseEvent e) {
    }

	public int setup(String s, ImagePlus imagePlus) {
        if(imagePlus==null) {IJ.showMessage("needs image opened"); return DONE; }
		inimg = Tools.convertToFloatImage(imagePlus);
		inimg.setTitle("inimg");
		inimgPath = imagePlus.getOriginalFileInfo().directory+imagePlus.getOriginalFileInfo().fileName;
		return DOES_8G+DOES_32+NO_CHANGES;
	}

}