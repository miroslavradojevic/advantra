package detection;

import aux.Tools;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.*;
import ij.plugin.filter.PlugInFilter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 10/8/13
 * Time: 2:10 PM
 */
public class JunctionTrace implements PlugInFilter, MouseListener {

    private ImagePlus   imp;
    private ImageCanvas canvas;
    //private ImageWindow win;
    private double	D = 3;
    private double	s = 1.5;
    private int	    N = 1;
    private int	    Loops = 1;
    private float[] exProf = null;
    private int		wStdRatioToD = 4;
    private int pointsLS = 150;
    private double[] startLS 	= new double[pointsLS];
    private double[] finishLS 	= new double[pointsLS];
    float[] peakIdxs = null;
    Overlay juncTrace = new Overlay();

    private static float 	Deg2Rad = (float) (Math.PI/180f);
    private static float 	Rad2Deg = (float) (180f/Math.PI);

    public int setup(String s, ImagePlus imagePlus) {

        // parameters
        D       =  	Prefs.get("advantra.critpoint.D", 3);
        //s       = 	Prefs.get("advantra.critpoint.scale", 			1.5);

        GenericDialog gd = new GenericDialog("JUNCTION TRACE");
        gd.addNumericField("neuron diameter ",  D, 0, 10, " pix");
        //gd.addNumericField("scale ",  s, 1, 10, " x(neuron diameter)");

        gd.showDialog();
        if (gd.wasCanceled()) return DONE;

        D       =  gd.getNextNumber();
        Prefs.set("advantra.critpoint.D",   D);

        //s               = (float) gd.getNextNumber();
        //Prefs.set("advantra.critpoint.scale",   s);

        //turnOn();

        IJ.setTool("hand");

        if(imagePlus!=null) {
            int dim[] = imagePlus.getDimensions();
            //this.dim = new int[dim.length];
            //for (int i=0; i<dim.length; i++) this.dim[i] = dim[i];
            boolean is2d = dim[0] > 0 && dim[1] > 0 && dim[2] == 1 && dim[3] == 1 && dim[4] == 1;
            if(!is2d) {
                IJ.error("This plugin works with 2d images without channels or frames.");
                return DONE;
            }
            else {
                this.imp = Tools.convertToFloatImage(imagePlus);
                canvas = imagePlus.getCanvas();
                //win = new ImageWindow(imagePlus);
                return DOES_ALL+NO_CHANGES;
            }
        }
        else {
            return DONE;
        }
    }

    public void run(ImageProcessor imageProcessor) {
        canvas.addMouseListener(this);
        IJ.log("junction tracer...\nstart with mouse click!");
    }

    public void mouseClicked(MouseEvent e) {

        int offscreenX = canvas.offScreenX(e.getX());
        int offscreenY = canvas.offScreenY(e.getY());

        //canvas.zoomIn(offscreenX, offscreenY);

//        IJ.log(offscreenX+" : "+offscreenY);
//        long t1, t2;
//        t1 = System.currentTimeMillis();
//        t2 = System.currentTimeMillis();
//        IJ.log("found "+peakIdxs.length+" peaks, "+((t2-t1)/1000f)+" sec. elapsed");

        juncTrace.clear();

        // overlay central location
        double R = .5;
        OvalRoi ovroi = new OvalRoi(offscreenX-R+.5, offscreenY-R+.5, 2*R, 2*R);
        ovroi.setStrokeWidth(1);
        ovroi.setStrokeColor(Color.YELLOW);
        juncTrace.add(ovroi);
        int cc = 0;

        while(cc<N) {

            exProf = Profiler.extractProfile(D, s, wStdRatioToD, offscreenX,  offscreenY, (FloatProcessor) imp.getProcessor());
            peakIdxs = Analyzer.extractPeakIdxs(exProf, startLS, finishLS);
            cc++;

            if (peakIdxs!=null) {
                for (int ii=0; ii<peakIdxs.length; ii++) {

                    float resolDeg = Profiler.getResolDeg(s);
                    float angleRad = peakIdxs[ii] * resolDeg * Deg2Rad;

                    double atX = offscreenX + s*D * Math.cos( angleRad );
                    double atY = offscreenY - s*D * Math.sin( angleRad );

                    ovroi = new OvalRoi(atX-R+.5, atY-R+.5, 2*R, 2*R);
                    ovroi.setStrokeWidth(2);
                    ovroi.setStrokeColor(Color.RED);
                    juncTrace.add(ovroi);


                }
            }

        }



        IJ.log("window width: "+canvas.getImage().getWindow().getMaximumBounds().width);
        IJ.log("window height: "+canvas.getImage().getWindow().getMaximumBounds().height);




        // overlay those peaks

        canvas.setOverlay(juncTrace);

    }

    @Override
    public void mousePressed(MouseEvent e) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void mouseReleased(MouseEvent e) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void mouseEntered(MouseEvent e) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void mouseExited(MouseEvent e) {
        //To change body of implemented methods use File | Settings | File Templates.
    }


}
