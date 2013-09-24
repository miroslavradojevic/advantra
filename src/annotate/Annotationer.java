package annotate;

import aux.AnalyzeCSV;
import aux.Tools;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.*;
import ij.io.FileInfo;
import ij.plugin.PlugIn;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.*;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 9/22/13
 * Time: 10:22 PM
 */
public class Annotationer  implements PlugInFilter, MouseListener {

	/**
	 * searches folder where the file is to open .csv with annotations
	 * if there is not one - it creates new .csv annotation file
	 * with the same name as image
	 */

	ImagePlus   inimg;
	ImageCanvas canvas;
	FileInfo	inimgInfo;
	String      imgPath;
	String      annFile;
	File 		fread;

	public int setup(String s, ImagePlus imagePlus) {

		if(imagePlus==null) return DONE;

		inimg = Tools.convertToFloatImage(imagePlus);

		canvas = imagePlus.getCanvas();

		inimgInfo = imagePlus.getOriginalFileInfo();

		fread = new File(inimgInfo.directory+inimgInfo.fileName);

		if (!fread.isFile()) return DONE;

		imgPath = fread.getAbsolutePath();

		return DOES_8G+DOES_32+NO_CHANGES;
	}

	public void run(ImageProcessor imageProcessor) {

		File dir = new File(fread.getParent());
		String[] dirList = dir.list();

		String pattern = inimgInfo.fileName.substring(0,inimgInfo.fileName.length()-4)+".csv";

		boolean found = false;
		annFile = inimgInfo.directory + pattern ;

		for (int i=0; i<dirList.length; i++) {

			if (dirList[i].equals(pattern)) {

				found = true;
				break;

			}

		}

		double[][] A = null;

		if (found) {
			// read it
			AnalyzeCSV readCSV;
			readCSV = new AnalyzeCSV(annFile);
			A = readCSV.readLnDouble(2);
		}
		else {
		    // initialize new empty one
			try {
				PrintWriter writer = new PrintWriter(annFile);
				writer.println("#");
				writer.close();
			}
			catch (FileNotFoundException ex) {}

		}

		Overlay o = new Overlay();

		if (A!=null) {
			for (int i = 0; i < A.length; i++) {
				Roi pt = new PointRoi(A[i][0]+0.5, A[i][1]+0.5);
				//pt.setName("line,"+IJ.d2s(i+1,0)+",("+A[i][0]+","+A[i][1]+")");
				o.add(pt);
			}
		}

		//o.drawNames(true);

		canvas.setOverlay(o);
		//inimg.show();
		canvas.addMouseListener(this);
		IJ.setTool("hand");


	}

//		annFile = MyOpener.open("Open annotation file [ ."+ext+" extension ]", false);
//
//		if (annFile==null) return;
//
//		String extension = "";
//		int dotIdx = annFile.lastIndexOf('.');
//		if (dotIdx > 0)
//			extension = annFile.substring(dotIdx+1);
//
//		double[][] A = null;
//
//		if (!(extension.equals(ext))) {
//			IJ.log("wrong extension: " + extension);
//			return;
//		}
//		else {
//			AnalyzeCSV readCSV;
//			readCSV = new AnalyzeCSV(annFile);
//			A = readCSV.readLnDouble(2);
//		}
//
//		imgPath = MyOpener.open("Open image file", false);
//
//		if (imgPath==null) return;
//
//		inimg = new ImagePlus(imgPath);
//
//		if (inimg.getDimensions()[3]>1) {
//			IJ.log("image has to be 2d");
//			return; // there was more than 1 slice
//		}
//
//		Overlay o = new Overlay();
//		if (A!=null) {
//			for (int i = 0; i < A.length; i++) {
//				Roi pt = new PointRoi(A[i][0]+0.5, A[i][1]+0.5);
//				pt.setName("line,"+IJ.d2s(i+1,0)+",("+A[i][0]+","+A[i][1]+")");
//				o.add(pt);
//			}
//		}
//
//		o.drawNames(true);
//
//		inimg.setOverlay(o);
//		inimg.show();
//		inimg.getCanvas().addMouseListener(this);
//		IJ.setTool("hand");



	public void mouseClicked(MouseEvent e) {



//		ImageCanvas srcCanv = (ImageCanvas) e.getSource();
		int atX = 	canvas.offScreenX(e.getX());
		int atY = 	canvas.offScreenY(e.getY());

		//double atXsubp = e.getPoint().getX();
		//atXsubp = atX + atXsubp - Math.floor(atXsubp);
		//double atYsubp = e.getPoint().getY();
		//atYsubp = atY + atYsubp - Math.floor(atYsubp);

		Roi pt = new PointRoi(atX+0.5, atY+0.5);
		//pt.setName("new?("+atX+","+atY+")");

		Overlay currOverlay = canvas.getImage().getOverlay();

		currOverlay.add(pt);
		canvas.getImage().updateAndDraw();

		GenericDialog gd = new GenericDialog("OK?");
		gd.showDialog();
		if (gd.wasCanceled()) {
			currOverlay.remove(currOverlay.size()-1);
			canvas.getImage().updateAndDraw();
			return;
		}

		// append line
		try {
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(annFile, true)));
			out.println(IJ.d2s(atX, 0, 0)+","+IJ.d2s(atY, 0, 0));
			out.close();
		} catch (IOException e1) {
			//oh noes!
		}
	}

	public void mousePressed(MouseEvent e)  {}

	public void mouseReleased(MouseEvent e) {}

	public void mouseEntered(MouseEvent e)  {}

	public void mouseExited(MouseEvent e)   {}




}


