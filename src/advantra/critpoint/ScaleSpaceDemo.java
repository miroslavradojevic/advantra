package advantra.critpoint;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import advantra.feature.ScaleSpace;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.measure.Calibration;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import imagescience.image.Axes;
import imagescience.image.Coordinates;
import imagescience.image.FloatImage;
import imagescience.image.Image;

public class ScaleSpaceDemo implements PlugInFilter, MouseListener {

	Image img;
	ImagePlus scale_space;
	double[] sigma = null;
	double[] profile = null;
	
	PlotWindow pw;
	
	public void run(ImageProcessor arg0) {
		
		pw = new PlotWindow("L(sigma)", "sigma", "Intensity", new double[]{1,2}, new double[]{1,2});
		
		GenericDialog gd = new GenericDialog("Scale space");
		//gd.addStringField("folder with annotations", dir_path, 50);
		gd.addMessage("Choose scales (gaussian std.)");
		gd.addNumericField( "sigma start:", 			2, 	 0, 5, "" );
		gd.addNumericField( "sigma end  :", 			6, 	 0, 5, "" );	
		gd.addNumericField( "number of scales : ", 		10,  0, 5, "");
		
		gd.showDialog();
		if (gd.wasCanceled()) return;
		
		double 		sigma_1 		=  (double)gd.getNextNumber();
		double 		sigma_2 		=  (double)gd.getNextNumber();
		int 		nr				=  (int)gd.getNextNumber();
		
		//double[] scales = new double[];
		
		// scales
		sigma = new double[nr];
		for (int i = 0; i < nr; i++) {
			sigma[i] = (i==0)? sigma_1 : sigma_1+i*((sigma_2-sigma_1)/(nr-1));
		}
		
		
		
//		Vector<Image> imgs_at_scales = scsp_extractor.extract(img, sigma);
//		for (int i = 0; i < imgs_at_scales.size(); i++) {
//			imgs_at_scales.get(i).imageplus().show();
//		}
		
		ScaleSpace scsp_extractor = new ScaleSpace();
		scale_space = scsp_extractor.extractAsStack(img, sigma).imageplus();
		scale_space.show();
		scale_space.getWindow().getCanvas().addMouseListener(this);
		IJ.showMessage("click on the point ...");
		
	}

	public int setup(String arg0, ImagePlus imp) {
		if(imp!=null) {
					// reset calibration before going further 
					Calibration cal = new Calibration(imp);
					cal.pixelWidth = cal.pixelHeight = cal.pixelDepth = 1;
					cal.setUnit("pixel");
					imp.setCalibration(cal);
					img = new FloatImage(Image.wrap(imp));
		}
		return DOES_8G+DOES_16+DOES_32+NO_CHANGES;
	}

	public void mouseClicked(MouseEvent e) {
		
		int mouseX = scale_space.getWindow().getCanvas().offScreenX(e.getX());
		int mouseY = scale_space.getWindow().getCanvas().offScreenX(e.getY());
//		IJ.showMessage(
//				"X: "++
//				"," +
//				"Y: "+);
		
		// take the data out
		
		
		Image sc_sp = Image.wrap(scale_space);
		profile = new double[sc_sp.dimensions().z];
		sc_sp.axes(Axes.Z);
		Coordinates coord = new Coordinates();
		coord.y = mouseY;
		coord.x = mouseX;
		sc_sp.get(coord, profile);
		// show it in plot window
		//PlotWindow pw = new PlotWindow("L(sigma)", "sigma", "Intensity", sigma, L);
		//pw.addPoints(sigma, profile, PlotWindow.LINE);
		pw.drawPlot(new Plot("", "", "", sigma, profile));
		//pw.draw();
		
	}

	public void mousePressed(MouseEvent e) {}

	public void mouseReleased(MouseEvent e) {}

	public void mouseEntered(MouseEvent e) {}

	public void mouseExited(MouseEvent e) {}

}
