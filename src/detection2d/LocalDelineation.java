package detection2d;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.*;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import tracing2d.BayesianTracer2D;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

/**
 * Created by miroslav on 11-11-14.
 */
public class LocalDelineation implements PlugIn, MouseListener {

    // input
//    ImagePlus   img;
    ImageCanvas cnv;
    String      image_path;
    float[][] likelihood_xy;

    // classes used for processing
    SemiCircle scirc;

    // params
    float radius = 3;
    int Ni = 5;
    int Ns = 20;
    float sigma_deg = 60;

    // output
    float[][][] xt = new float[Ni][Ns][4];
    float[][]   wt = new float[Ni][Ns];

    public void run(String s) {

        // load the image through the menu
        String in_folder = Prefs.get("id.folder", System.getProperty("user.home"));
        OpenDialog.setDefaultDirectory(in_folder);
        OpenDialog dc = new OpenDialog("Select file");
        in_folder = dc.getDirectory();
        image_path = dc.getPath();
        if (image_path==null) return;
        Prefs.set("id.folder", in_folder);

        ImagePlus img = new ImagePlus(image_path);
        if(img==null) return;

        // set image as float[][]
        likelihood_xy = new float[img.getWidth()][img.getHeight()]; 	// x~column, y~row
        if (img.getType()== ImagePlus.GRAY8) {
            byte[] read = (byte[]) img.getProcessor().getPixels();
            for (int idx=0; idx<read.length; idx++) {
                inimg_xy[idx%ip_load.getWidth()][idx/ip_load.getWidth()] = (float) (read[idx] & 0xff);
            }
        }
        else if (ip_load.getType()==ImagePlus.GRAY32) {
            float[] read = (float[]) ip_load.getProcessor().getPixels();
            for (int idx=0; idx<read.length; idx++) {
                inimg_xy[idx%ip_load.getWidth()][idx/ip_load.getWidth()] = read[idx];
            }
        }
        else {
            IJ.log("image type not recognized");
            return;
        }

        scirc = new SemiCircle(radius);

        // show it and get the canvas
        img.show();
        cnv = img.getCanvas();
        cnv.addMouseListener(this);
        IJ.setTool("hand");


    }

    public void mouseClicked(MouseEvent e) {
        int clickX = cnv.offScreenX(e.getX());
        int clickY = cnv.offScreenY(e.getY());

        System.out.println("(x,y) = "+clickX+" , " +clickY);

        BayesianTracer2D.run(clickX, clickY,
                likelihood_xy, scirc, sigma_deg, Ni, Ns, xt, wt);

        System.out.println("done.");

        cnv.setOverlay(viz_xt(xt));

    }

    public static Overlay viz_xt(float[][][] xt)
    {

        float rad = .5f;

        Overlay ov = new Overlay();
        for (int i = 0; i < xt.length; i++) { // iterations
            for (int j = 0; j <xt[i].length; j++) {

                OvalRoi p = new OvalRoi(xt[i][j][0]-rad+.5, xt[i][j][1]-rad+.5, 2*rad, 2*rad);
                Line l = new Line(xt[i][j][0]+.5, xt[i][j][1]+.5, xt[i][j][0]+xt[i][j][2]+.5, xt[i][j][1]+xt[i][j][3]+.5);

                ov.add(p);
                ov.add(l);

            }
        }

        return ov;
    }

    public void mousePressed(MouseEvent e) {}
    public void mouseReleased(MouseEvent e) {}
    public void mouseEntered(MouseEvent e) {}
    public void mouseExited(MouseEvent e) {}

}
