import ij.IJ;
import ij.ImagePlus;
import ij.gui.ImageCanvas;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.plugin.PlugIn;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 6/21/13
 * Time: 11:05 AM
 */
public class ShowAnnotation implements PlugIn, MouseListener {

    ImagePlus   inimg;
    String      imgPath;
    String      annFile;

    public void run(String s) {

        annFile = MyOpener.open("Open annotation file", false);

        String extension = "";
        int dotIdx = annFile.lastIndexOf('.');
        if (dotIdx > 0)
            extension = annFile.substring(dotIdx+1);

        double[][] A = null;

        if (!extension.equals(new String("pos"))) {
            IJ.log("wrong extension: "+extension);
            return;
        }
        else {
            AnalyzeCSV readCSV;
            readCSV = new AnalyzeCSV(annFile);
            A = readCSV.readLn(2);
        }

        imgPath = MyOpener.open("Open image file", false);

        inimg = new ImagePlus(imgPath);

        Overlay o = new Overlay();
        if (A!=null) {
            for (int i = 0; i < A.length; i++) {
                Roi pt = new PointRoi(A[i][0]+0.5, A[i][1]+0.5);
                pt.setName("line,"+IJ.d2s(i+1,0)+",("+A[i][0]+","+A[i][1]+")");
                o.add(pt);
            }
        }

        o.drawNames(true);

        inimg.setOverlay(o);
        inimg.show();
        inimg.getCanvas().addMouseListener(this);
        IJ.setTool("hand");

    }

    public void mouseClicked(MouseEvent e) {
        ImageCanvas srcCanv = (ImageCanvas) e.getSource();
        int atX = 	srcCanv.offScreenX(e.getX());
        int atY = 	srcCanv.offScreenY(e.getY());

        // append line
        try {
            PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(annFile, true)));
            out.println(IJ.d2s(atX, 0, 0)+","+IJ.d2s(atY, 0, 0));
            out.close();
        } catch (IOException e1) {
            //oh noes!
        }

        Roi pt = new PointRoi(atX+0.5, atY+0.5);
        pt.setName("new,("+atX+","+atY+")");

        srcCanv.getImage().getOverlay().add(pt);
        srcCanv.getImage().updateAndDraw();
    }

    public void mousePressed(MouseEvent e) {
    }

    public void mouseReleased(MouseEvent e) {
    }

    public void mouseEntered(MouseEvent e) {
    }

    public void mouseExited(MouseEvent e) {
    }
}
