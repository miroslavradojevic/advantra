package profile;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.*;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.Arrays;
import java.util.Random;

import static ij.gui.Plot.*;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 7/31/13
 * Time: 10:43 AM
 */
public class DynamicProfileInspector implements PlugInFilter, MouseListener, MouseMotionListener, KeyListener, WindowListener {

	private ImagePlus   imp;
	private ImageCanvas canvas;
	private PlotWindow  pw;
	private int dim[];

	private double	D = 3;
	private double	s = 1.5;
	private float[] exProf = null; // profile

	// MS
	private int pointsMS = 200;
	private double[] startMS 	= new double[pointsMS];
	private double[] finishMS 	= new double[pointsMS];
	float[] peakIdxs = null;
	// will use MS parameters hardcoded in Analyzer class

	public int setup(String s, ImagePlus imp)
	{

		if(imp!=null) {

			int dim[] = imp.getDimensions();
			this.dim = new int[dim.length];
			for (int i=0; i<dim.length; i++) this.dim[i] = dim[i];
			boolean is2d = dim[0] > 0 && dim[1] > 0 && dim[2] == 1 && dim[3] == 1 && dim[4] == 1;
			if(!is2d) {
				IJ.error("This plugin only works with 2d images (one slice) without channels or frames.");
				return DONE;
			}
			else {
				this.imp = Tools.convertToFloatImage(imp);
				canvas = imp.getCanvas();
				return DOES_ALL+NO_CHANGES;
			}

		}
		else {
			return DONE;
		}

	}

	public void run(ImageProcessor imageProcessor)
	{

        IJ.log("profiler + peak detector live...");

		D       =  	Prefs.get("advantra.critpoint.D", 	3);
		s       = 	Prefs.get("advantra.critpoint.scale", 			1.5);

		GenericDialog gd = new GenericDialog("JUNCTION DET.");
		gd.addNumericField("neuron diameter ",  D, 0, 10, " pix");
		gd.addNumericField("scale ",  s, 1, 10, " x(neuron diameter)");

		gd.showDialog();
		if (gd.wasCanceled()) return;

		D       =  gd.getNextNumber();
		Prefs.set("advantra.critpoint.D",   D);

		s               = (float) gd.getNextNumber();
		Prefs.set("advantra.critpoint.scale",   s);

		turnOn();

		IJ.setTool("hand");

		IJ.log("move mouse...");

	}

	private void turnOn()
	{
		canvas.addMouseListener(this);
		canvas.addMouseMotionListener(this);
		canvas.addKeyListener(this);
	}

	private void turnOff() {
		canvas.removeMouseMotionListener(this);
		canvas.removeMouseListener(this);
		canvas.removeKeyListener(this);
	}

	public void keyPressed(KeyEvent e)
	{

		// Catch the event that enables or disables plot updating
		if (e.getKeyCode() == KeyEvent.VK_Q) {
			IJ.log("exporting ");

		}

		if (e.getKeyCode() == KeyEvent.VK_U) {
//
//			if (pw!=null) {
//
//				// export to csv
//				String fileName = "profile.csv";
//
//				// empty the file
//				PrintWriter writer = null;
//				try {
//					writer = new PrintWriter(fileName);
//					writer.print("");
//					writer.close();
//				} catch (FileNotFoundException ex) {}
//
//				try {
//					PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(fileName, true)));
//
//					out.println("Angle, Response");
//
//					int idx = 0;
//					for (int aDeg = 0; aDeg<360; aDeg+=Profiler.resolDeg) {
//
//						out.println(aDeg+", "+exProf[idx++]);
//
//					}
//
//					out.close();
//
//				} catch (IOException e1) {}
//
//				IJ.showMessage("exported to : " + new File(fileName).getAbsolutePath());
//
//					/*
//					for (int angleCnt=0; angleCnt<Profiler.offsets.size(); angleCnt++) {
//
//						out.print(""+Profiler.offsets.get(angleCnt).size()+", "+(angleCnt*Profiler.resolDeg)+", ");
//
//						for (int idx=0; idx<Profiler.offsets.get(angleCnt).size(); idx++) {
//
//							out.print(
//											 IJ.d2s(Profiler.offsets.get(angleCnt).get(idx)[0], 	3)+", "+
//													 IJ.d2s(Profiler.offsets.get(angleCnt).get(idx)[1], 	3)+", "+
//													 IJ.d2s(Profiler.weights.get(angleCnt).get(idx), 	3)
//							);
//
//							if (idx==Profiler.offsets.get(angleCnt).size()-1) {
//								out.println("");
//							}
//							else {
//								out.print(", ");
//							}
//
//						}
//
//					}
//					*/
//
//				//out.println("");
//
//			}
//			else {
//
//				IJ.showMessage("no plot window opened");
//
//			}
//
		}

	}

	public void mouseClicked(MouseEvent e)
	{
		int offscreenX = canvas.offScreenX(e.getX());
		int offscreenY = canvas.offScreenY(e.getY());

        IJ.log(offscreenX+" : "+offscreenY);

		long t1, t2;

		t1 = System.currentTimeMillis();
		exProf = Profiler.extractProfile(D, s, offscreenX, offscreenY, (FloatProcessor) imp.getProcessor());
        peakIdxs = Analyzer.extractPeakIdxs(exProf, startMS, finishMS);
		t2 = System.currentTimeMillis();
		IJ.log("found "+peakIdxs.length+" peaks, "+((t2-t1)/1000f)+" sec. elapsed");

		// Fill in X axis (frame number)
        float[] x = new float[exProf.length];
        for (int i = 1; i <= x.length; i++)
            x[i - 1] = (i-1)*Profiler.getResolDeg(s);

//		float[] mm = Tools.getMinMax(exProf);

        // Prepare plot window
        Plot chart = new Plot("", "", "", x, exProf);
		chart.setSize(800, 400);
		// Add the points for prettier plots
		//chart.addPoints(x, exProf, Plot.CIRCLE);
        //chart.setLimits(0, 360, mm[0]-1, mm[1]+1);

		float[] tempx = new float[pointsMS]  ;
		float[] tempy = new float[pointsMS]  ;

		/*t1 = System.currentTimeMillis();
		peakIdxs = Analyzer.extractPeakIdxs(exProf, startMS, finishMS); // does max-shift
		t2 = System.currentTimeMillis();
		IJ.log(""+peakIdxs.length+" peaks! "+((t2-t1)/1000f)+" sec elapsed");*/

		for (int yy=0; yy<pointsMS; yy++) {
			tempx[yy] = (float) (startMS[yy] * Profiler.getResolDeg(s));
			tempy[yy] = (float) Tools.interp1Darray((float) startMS[yy], exProf);
		}
		chart.draw(); chart.setColor(Color.BLUE);
		chart.addPoints(tempx, 	tempy, X);

		for (int yy=0; yy<pointsMS; yy++) {
			tempx[yy] = (float) (finishMS[yy] * Profiler.getResolDeg(s));
			tempy[yy] = (float) Tools.interp1Darray((float) finishMS[yy], exProf);
		}
		chart.draw(); chart.setColor(Color.RED);
		chart.addPoints(tempx, 	tempy, PlotWindow.CIRCLE);

        if (pw == null) {
            pw = chart.show();
            pw.addWindowListener(this);
        } else
            pw.setTitle("");

        pw.drawPlot(chart);
		pw.setTitle("Profile, x = " + offscreenX + ", y = " + offscreenY + " : " );

	}

	private static Color getColor(int i)
	{
		if (i<=0) {
			return Color.RED;
		}
		else if (i==1) {
			return Color.GREEN;
		}
		else {
			return Color.BLUE;
		}
	}

	public void mouseMoved(MouseEvent e)
	{
		mouseClicked(e);
	}

	public void windowClosed(WindowEvent e)
	{
		turnOff();
	}

	/*
	unused
	 */
	public void mousePressed(MouseEvent e) {}
	public void mouseReleased(MouseEvent e) {}
	public void mouseEntered(MouseEvent e) {}
	public void mouseExited(MouseEvent e) {}
	public void mouseDragged(MouseEvent e) {}
	public void keyTyped(KeyEvent e) {}
	public void keyReleased(KeyEvent e) {}

	@Override
	public void windowOpened(WindowEvent e) {
		//To change body of implemented methods use File | Settings | File Templates.
	}

	@Override
	public void windowClosing(WindowEvent e) {
		//To change body of implemented methods use File | Settings | File Templates.
	}

	@Override
	public void windowIconified(WindowEvent e) {
		//To change body of implemented methods use File | Settings | File Templates.
	}

	@Override
	public void windowDeiconified(WindowEvent e) {
		//To change body of implemented methods use File | Settings | File Templates.
	}

	@Override
	public void windowActivated(WindowEvent e) {
		//To change body of implemented methods use File | Settings | File Templates.
	}

	@Override
	public void windowDeactivated(WindowEvent e) {
		//To change body of implemented methods use File | Settings | File Templates.
	}
}
