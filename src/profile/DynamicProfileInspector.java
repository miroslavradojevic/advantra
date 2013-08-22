package profile;

import ij.IJ;
import ij.ImagePlus;
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

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 7/31/13
 * Time: 10:43 AM
 */
public class DynamicProfileInspector implements PlugInFilter, MouseListener, MouseMotionListener,
														KeyListener {

	private ImagePlus   imp;
	private ImageCanvas canvas;
	private PlotWindow  pw;
	private int dim[];


//	private JFrame jf;
//	private JLabel enabled;
//	private JTextField D_JInput;
//	private JTextField s_JInput;


    private boolean moveFlag = true;
	private double	D = 4;
	private double	s = 1.0;
	private float[] exProf = null;

    private Overlay ov = new Overlay();
//    private float[] colinP = new float[50];
//    private float[] colinN = new float[50];
//    private float[] iJumpP = new float[50];
//    private float[] iJumpN = new float[50];
//    private int cntP = 0, cntN = 0;
    private PlotWindow fw;

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
				return DOES_ALL;
			}

		}
		else {
			return DONE;
		}

	}

	public void run(ImageProcessor imageProcessor)
	{

        IJ.log("test profiler live...");



//		enabled = new JLabel();
//		enabled.setBounds(10, 10, 100, 25);
//		changeEnabledLabel();
//
//		D_JInput = new JTextField(Double.toString(D), 15);
//		D_JInput.setName("D");
//		D_JInput.setBounds(120, 10, 100, 25);
//		D_JInput.addActionListener(this);
//
//		s_JInput = new JTextField(Double.toString(s), 15);
//		s_JInput.setName("s");
//		s_JInput.setBounds(230, 10, 100, 25);
//		s_JInput.addActionListener(this);
//
//		jf = createFrame();
//		jf.addWindowListener(this);
//		jf.setVisible(true);
		turnOn();

		// profiler ready
		ByteProcessor mask = new ByteProcessor(dim[0], dim[1]);
		for (int i=0; i<dim[0]*dim[1]; i++) mask.set(i, (byte)255);

        Profiler2.loadTemplate(imp.getProcessor(), mask);
        Profiler2.loadParams(D, s, true);

		//Profiler.loadTemplate(imp.getProcessor(), mask);
		//Profiler.loadParams(D, s, true);   // true to see how profiles are created (for figures mainly)

//        Arrays.fill(colinN, Float.NaN);
//        Arrays.fill(colinP, Float.NaN);
//        Arrays.fill(iJumpN, Float.NaN);
//        Arrays.fill(iJumpP, Float.NaN);
//        cntP = 0;
//        cntN = 0;

		IJ.setTool("hand");
        IJ.log("DONE");

	}

//	private void changeEnabledLabel()
//	{
//		if (moveFlag) {
//			enabled.setText("Enabled");
//			enabled.setForeground(Color.green);
//		} else {
//			enabled.setText("Disabled");
//			enabled.setForeground(Color.red);
//		}
//	}

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

//	private JFrame createFrame()
//	{
//		JFrame frame = new JFrame();
//
//		frame.setTitle("Intensity Profile Inspector");
//		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
//		frame.setBounds(100, 100, 400, 150);
//
//		JPanel contentPane = new JPanel();
//		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
//		frame.setContentPane(contentPane);
//		contentPane.setLayout(null);
//
//
//
//		contentPane.add(enabled);
//		contentPane.add(D_JInput);
//		contentPane.add(s_JInput);
//
//		JLabel enabledLabel = new JLabel("LIVE MODE");
//		enabledLabel.setBounds(10, 40, 100, 25);
//		contentPane.add(enabledLabel);
//
//		JLabel D_JInputLabel = new JLabel("NEURON D.");
//		D_JInputLabel.setBounds(120, 40, 100, 25);
//		contentPane.add(D_JInputLabel);
//
//		JLabel s_JInputLabel = new JLabel("RANGE(xNEURON D.)");
//		s_JInputLabel.setBounds(230, 40, 150, 25);
//		contentPane.add(s_JInputLabel);
//
////		JCheckBox chckboxInvertValues = new JCheckBox("Invert values");
////		chckboxInvertValues.setBounds(10, 50, 97, 23);
////		chckboxInvertValues.addActionListener(this);
////		contentPane.add(chckboxInvertValues);
//
//		JLabel lblPressCtrl = new JLabel("Press 'q' for enabling or disabling" +
//												 " continous plotting");
//		lblPressCtrl.setFont(new Font("Tahoma", Font.PLAIN, 8));
//		lblPressCtrl.setBounds(10, 80, 216, 15);
//		contentPane.add(lblPressCtrl);
//
//		JLabel lblExport = new JLabel("Press 'u' to export" +
//												 " filters and profile");
//		lblExport.setFont(new Font("Tahoma", Font.PLAIN, 8));
//		lblExport.setBounds(10, 80+15, 216, 14);
//		contentPane.add(lblExport);
//
//		return frame;
//	}

	public void keyPressed(KeyEvent e)
	{
//
//		// Catch the event that enables or disables plot updating
//		if (e.getKeyCode() == KeyEvent.VK_Q) {
//			moveFlag = !moveFlag;
////			changeEnabledLabel();
//            if (moveFlag) {
//                ov.clear();
//                canvas.getImage().setOverlay(ov);
////                Arrays.fill(colinN, Float.NaN);
////                Arrays.fill(colinP, Float.NaN);
////                Arrays.fill(iJumpN, Float.NaN);
////                Arrays.fill(iJumpP, Float.NaN);
////                cntP = 0;
////                cntN = 0;
//            }
//		}
//
//		if (e.getKeyCode() == KeyEvent.VK_U) {
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
//		}
//
	}

	public void mouseClicked(MouseEvent e)
	{
		int offscreenX = canvas.offScreenX(e.getX());
		int offscreenY = canvas.offScreenY(e.getY());

        IJ.log(offscreenX+" : "+offscreenY);
		//exProf = Profiler2.extractProfile(D, s, offscreenX, offscreenY, (FloatProcessor) imp.getProcessor());

//        Analyzer.extractPeakIdxs(exProf); // does mean-shift,
        // this one returns the array now take it and see what you can print,
        // don't use Analyzer.peakIdx

        // Fill in X axis (frame number)
        float[] x = new float[exProf.length];
        for (int i = 1; i <= x.length; i++)
            x[i - 1] = (i-1)*Profiler.resolDeg;

        // Prepare plot window
        Plot chart = new Plot("", "", "", x, exProf);
        //float[] mm = Tools.getMinMax(exProf);
        //chart.setLimits(0, 360, mm[0]-1, mm[1]+1);
        chart.setSize(900, 450);

        if (pw == null) {
            pw = chart.show();
            //pw.addWindowListener(this);
        } else
            pw.setTitle("Profile, x = " + offscreenX + ", y = " + offscreenY);

        // Add the points for prettier plots
        chart.addPoints(x, exProf, PlotWindow.CIRCLE);

        /*
        // Add MS convergence points (use static method for this, as everything is now in Analyzer.extractPeakIdxs, needs to return)
        float[] xMS = new float[Analyzer.nrPoints];
        float[] yMS = new float[Analyzer.nrPoints];
        float[] yLW;
        for (int i=0; i<Analyzer.nrPoints; i++) {
            xMS[i] = Analyzer.convIdx.get(0).get(0)[i] * Profiler.resolDeg;
            yMS[i] = (float) Tools.interp1Darray(Analyzer.convIdx.get(0).get(0)[i], exProf);
        }
        chart.draw();
        chart.setColor(Color.GREEN);
        chart.addPoints(xMS, yMS, Plot.X);
*/

        //Add MS plot
//        String domes = "";
//        if (Analyzer.peakIdx.get(0).get(0).size() != 0) { // >=3 peaks
//            chart.draw();
//            chart.setColor(Color.RED);
//            chart.setLineWidth(5);
//            xMS = new float[3];
//            yMS = new float[3];
//            yLW = new float[3];
//			int[] roundedPeaks = new int[]{
//												  Math.round(Analyzer.peakIdx[0][0][0]),
//												  Math.round(Analyzer.peakIdx[0][0][1]),
//												  Math.round(Analyzer.peakIdx[0][0][2])
//			};
//            for (int i=0; i<3; i++) {
//                xMS[i] = Analyzer.peakIdx[0][0][i] * Profiler.resolDeg;
//                yMS[i] = (float) Tools.interp1Darray(Analyzer.peakIdx[0][0][i], exProf);
//                yLW[i] = Tools.findNextStationaryValue(Analyzer.peakIdx[0][0][i], roundedPeaks, exProf);
//                chart.drawLine(xMS[i], yLW[i], xMS[i], yMS[i]);
//                domes += "   " + IJ.d2s((yMS[i]>yLW[i])?(yMS[i]-yLW[i]):0, 2) + "  ";
//            }
            //chart.addPoints(xMS, yMS, Plot.BOX);

//        }

        pw.drawPlot(chart);
		pw.setTitle("Profile, x = " + offscreenX + ", y = " + offscreenY + " : " );

	}

    private  float colin(float[] center, float[][] pos) {


/*        // max area of the chosen triangle  3 out of 4 combinations, pos has to have length 3
        float maxTriangleArea = Float.MIN_VALUE;
        float currTriangleArea;

        // 0 : 0 1
        currTriangleArea = Tools.triangleArea(center[0], center[1], pos[0][0], pos[0][1], pos[1][0], pos[1][1]); // 0, 1
        if (currTriangleArea>maxTriangleArea) maxTriangleArea = currTriangleArea;
        // 0 : 0 2
        currTriangleArea = Tools.triangleArea(center[0], center[1], pos[0][0], pos[0][1], pos[2][0], pos[2][1]); // 0, 2
        if (currTriangleArea>maxTriangleArea) maxTriangleArea = currTriangleArea;
        // 0 : 1 2
        currTriangleArea = Tools.triangleArea(center[0], center[1], pos[1][0], pos[1][1], pos[2][0], pos[2][1]); // 1, 2
        if (currTriangleArea>maxTriangleArea) maxTriangleArea = currTriangleArea;
        //   : 0 1 2
        currTriangleArea = Tools.triangleArea(pos[0][0], pos[0][1], pos[1][0], pos[1][1], pos[2][0], pos[2][1]); // 0, 1, 2
        if (currTriangleArea>maxTriangleArea) maxTriangleArea = currTriangleArea;
        return maxTriangleArea;*/

        // new version
        float x10 = pos[0][0] - center[0];
        float y10 = pos[0][1] - center[1];

        float x31 = pos[2][0] - pos[0][0];
        float y31 = pos[2][1] - pos[0][1];

        float x21 = pos[1][0] - pos[0][0];
        float y21 = pos[1][1] - pos[1][1];

        float x32 = pos[2][0] - pos[1][0];
        float y32 = pos[2][1] - pos[1][1];

        return (float) (0.5*(1+(x21*x32+y21*y32)/(Math.sqrt(x21*x21+y21*y21)*Math.sqrt(x32*x32+y32*y32)))*0.5*(1+(x10*x31+y10*y31)/(Math.sqrt(x10*x10+y10*y10)*Math.sqrt(x31*x31+y31*y31))));

    }

    private  float ijump(float[] peakH, float[] peakL) {
        float out = 0;
        for (int i=0; i<peakH.length; i++) {
            out += (peakH[i]>peakL[i])? peakH[i]-peakL[i] : 0 ;
        }
        return out/peakH.length;
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
		//if (moveFlag)
			mouseClicked(e);
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

}
