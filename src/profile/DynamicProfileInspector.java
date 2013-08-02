package profile;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.ImageCanvas;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.*;
import java.io.*;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 7/31/13
 * Time: 10:43 AM
 */
public class DynamicProfileInspector implements PlugInFilter, ActionListener,
														WindowListener, MouseListener, MouseMotionListener,
														KeyListener {

	private ImagePlus imp;
	private ImageCanvas canvas;
	private PlotWindow pw;
	private int dim[];
	private JFrame jf;
	private JLabel enabled;
	private JTextField D_JInput;
	private JTextField s_JInput;
	private boolean moveFlag = true;
	private double	D = 3;
	private double	s = 1.5;
	private float[] exProf = null;

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
//				cal = imp.getCalibration();
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

		enabled = new JLabel();
		enabled.setBounds(10, 10, 100, 25);
		changeEnabledLabel();

		D_JInput = new JTextField(Double.toString(D), 15);
		D_JInput.setName("D");
		D_JInput.setBounds(120, 10, 100, 25);
		D_JInput.addActionListener(this);

		s_JInput = new JTextField(Double.toString(s), 15);
		s_JInput.setName("s");
		s_JInput.setBounds(230, 10, 100, 25);
		s_JInput.addActionListener(this);

		jf = createFrame();
		jf.addWindowListener(this);
		jf.setVisible(true);
		turnOn();

		// profiler ready
		ByteProcessor mask = new ByteProcessor(dim[0], dim[1]);
		for (int i=0; i<dim[0]*dim[1]; i++) mask.set(i, (byte)255);
		Profiler.loadTemplate(imp.getProcessor(), mask);
		Profiler.loadParams(D, s, false);   // true to see how profiles are created (for figures mainly)

		IJ.setTool("hand");

	}

	private void changeEnabledLabel()
	{
		if (moveFlag) {
			enabled.setText("Enabled");
			enabled.setForeground(Color.green);
		} else {
			enabled.setText("Disabled");
			enabled.setForeground(Color.red);
		}
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

	private JFrame createFrame()
	{
		JFrame frame = new JFrame();

		frame.setTitle("Intensity Profile Inspector");
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		frame.setBounds(100, 100, 400, 150);

		JPanel contentPane = new JPanel();
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		frame.setContentPane(contentPane);
		contentPane.setLayout(null);



		contentPane.add(enabled);
		contentPane.add(D_JInput);
		contentPane.add(s_JInput);

		JLabel enabledLabel = new JLabel("LIVE MODE");
		enabledLabel.setBounds(10, 40, 100, 25);
		contentPane.add(enabledLabel);

		JLabel D_JInputLabel = new JLabel("NEURON D.");
		D_JInputLabel.setBounds(120, 40, 100, 25);
		contentPane.add(D_JInputLabel);

		JLabel s_JInputLabel = new JLabel("RANGE(xNEURON D.)");
		s_JInputLabel.setBounds(230, 40, 150, 25);
		contentPane.add(s_JInputLabel);

//		JCheckBox chckboxInvertValues = new JCheckBox("Invert values");
//		chckboxInvertValues.setBounds(10, 50, 97, 23);
//		chckboxInvertValues.addActionListener(this);
//		contentPane.add(chckboxInvertValues);

		JLabel lblPressCtrl = new JLabel("Press 'q' for enabling or disabling" +
												 " continous plotting");
		lblPressCtrl.setFont(new Font("Tahoma", Font.PLAIN, 8));
		lblPressCtrl.setBounds(10, 80, 216, 15);
		contentPane.add(lblPressCtrl);

		JLabel lblExport = new JLabel("Press 'u' to export" +
												 " filters and profile");
		lblExport.setFont(new Font("Tahoma", Font.PLAIN, 8));
		lblExport.setBounds(10, 80+15, 216, 14);
		contentPane.add(lblExport);

		return frame;
	}

	public void actionPerformed(ActionEvent e)
	{

		String type = e.getSource().getClass().getName();

		if (type.equals("javax.swing.JCheckBox")) {
			JCheckBox jcb = (JCheckBox)e.getSource();
//			out = jcb.isSelected();
		}

		if (type.equals("javax.swing.JTextField")) {

			JTextField jtf = (JTextField)e.getSource();

			if (jtf.getName()=="D") {
				D = Double.valueOf(jtf.getText());
			}

			if (jtf.getName()=="s") {
				s = Double.valueOf(jtf.getText());
			}

			Profiler.loadParams(D, s, false);

		}

	}

	public void keyPressed(KeyEvent e)
	{

		// Catch the event that enables or disables plot updating
		if (e.getKeyCode() == KeyEvent.VK_Q) {
			moveFlag = !moveFlag;
			changeEnabledLabel();
		}

		if (e.getKeyCode() == KeyEvent.VK_U) {

			if (pw!=null) {

				// export to csv
				String fileName = "profile.csv";

				// empty the file
				PrintWriter writer = null;
				try {
					writer = new PrintWriter(fileName);
					writer.print("");
					writer.close();
				} catch (FileNotFoundException ex) {}

				try {
					PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(fileName, true)));

					out.println("Angle, Response");

					int idx = 0;
					for (int aDeg = 0; aDeg<360; aDeg+=Profiler.resolDeg) {

						out.println(aDeg+", "+exProf[idx++]);

					}

					out.close();

				} catch (IOException e1) {}

				IJ.showMessage("exported to : " + new File(fileName).getAbsolutePath());

					/*
					for (int angleCnt=0; angleCnt<Profiler.offsets.size(); angleCnt++) {

						out.print(""+Profiler.offsets.get(angleCnt).size()+", "+(angleCnt*Profiler.resolDeg)+", ");

						for (int idx=0; idx<Profiler.offsets.get(angleCnt).size(); idx++) {

							out.print(
											 IJ.d2s(Profiler.offsets.get(angleCnt).get(idx)[0], 	3)+", "+
													 IJ.d2s(Profiler.offsets.get(angleCnt).get(idx)[1], 	3)+", "+
													 IJ.d2s(Profiler.weights.get(angleCnt).get(idx), 	3)
							);

							if (idx==Profiler.offsets.get(angleCnt).size()-1) {
								out.println("");
							}
							else {
								out.print(", ");
							}

						}

					}
					*/

				//out.println("");

			}
			else {

				IJ.showMessage("no plot window opened");

			}

		}

	}

	public void mouseClicked(MouseEvent e)
	{
		int offscreenX = canvas.offScreenX(e.getX());
		int offscreenY = canvas.offScreenY(e.getY());

		exProf = Profiler.extractProfile(offscreenX, offscreenY);

		if (exProf != null) {

			// Fill in X axis (frame number)
			float[] x = new float[exProf.length];
			for (int i = 1; i <= x.length; i++)
				x[i - 1] = (i-1)*Profiler.resolDeg;

			// Prepare plot window
			Plot chart = new Plot("Profile,x = "
										  + offscreenX + ", y = " + offscreenY,
										 "Frame number", "Intensity", x, exProf);
			if (pw == null) {
				pw = chart.show();
				//pw.addWindowListener(this);
			} else
				pw.setTitle("Profile, x = " + offscreenX + ", y = " + offscreenY);

			// Add the points for prettier plots
			chart.addPoints(x, exProf, PlotWindow.CIRCLE);
			pw.drawPlot(chart);

		}

	}

	public void mouseMoved(MouseEvent e)
	{
		if (moveFlag)
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
	public void windowOpened(WindowEvent e) {}
	public void windowClosing(WindowEvent e) {}
	public void windowIconified(WindowEvent e) {}
	public void windowDeiconified(WindowEvent e) {}
	public void windowActivated(WindowEvent e) {}
	public void windowDeactivated(WindowEvent e) {}
	public void keyTyped(KeyEvent e) {}
	public void keyReleased(KeyEvent e) {}

}
