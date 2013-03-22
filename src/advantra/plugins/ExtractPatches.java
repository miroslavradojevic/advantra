package advantra.plugins;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.Vector;

import advantra.feature.MyHessian;

import weka.core.Utils;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.measure.Calibration;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import imagescience.image.Axes;
import imagescience.image.Coordinates;
import imagescience.image.Dimensions;
import imagescience.image.FloatImage;
import imagescience.image.Image;

public class ExtractPatches implements PlugInFilter, MouseListener {

	ImagePlus 	img;
	Image 		L1scales, L2scales, v1scales, v2scales; // each will be stack with each layer corresponding to a scale
	int 		ptch_size;
	double 		s_start;
	double 		s_end;
	int 		s_number;
	double[]	s;
	
	public void run(ImageProcessor arg0) {
		
		// extract hessian eigen values and eigen vectors for different scales
		Image im = new FloatImage(Image.wrap(img));
		
		Dimensions dims 		= im.dimensions();
		
		System.out.println("dims(x,y,z) = "+dims.x+" , "+dims.y+" , "+dims.z);
		
		Dimensions new_dims 	= new Dimensions(img.getWidth(), img.getHeight(), s.length);
		
		L1scales = new FloatImage(new_dims); L1scales.axes(Axes.X);
		L2scales = new FloatImage(new_dims); L2scales.axes(Axes.X);
		v1scales = new FloatImage(new_dims); v1scales.axes(Axes.X);
		v2scales = new FloatImage(new_dims); v2scales.axes(Axes.X);
		
		MyHessian my_hess = new MyHessian();
		
		double[] aL1 = new double[dims.x];
		double[] aL2 = new double[dims.x];
		double[] aV11 = new double[dims.x];
		double[] aV12 = new double[dims.x];
		
		
		Coordinates coords 	= new Coordinates();
		
		for (int i = 0; i < s.length; i++) {
			
			Vector<Image> hess = my_hess.eigs(im.duplicate(), s[i], true);
			
			// assign values to layers of L1, L2, v1, v2
			Image L2 	= hess.get(0); // L2
			L2.axes(Axes.X);
			Image L1 	= hess.get(1); // L1
			L1.axes(Axes.X);
			Image V11 	= hess.get(2); // V11
			V11.axes(Axes.X);
			Image V12 	= hess.get(3); // V12
			V12.axes(Axes.X);
			
			for (coords.y=0; coords.y<dims.y; ++coords.y) {
				
				//System.out.println(""+coords.y+" of "+dims.y);
				
				coords.z = 0;
				L2.get(coords,aL2);
				L1.get(coords,aL1);
				V11.get(coords,aV11);
				V12.get(coords,aV12);
				
				coords.z = i;
				L2scales.set(coords, aL2);
				L1scales.set(coords, aL1);
				v1scales.set(coords, aV11);
				v2scales.set(coords, aV12);
				
			}
			
		}
		
		L1scales.imageplus().show();
		L2scales.imageplus().show();
		v1scales.imageplus().show();
		v2scales.imageplus().show();
		
		img.getWindow().getCanvas().addMouseListener(this);
		System.out.format("Click!\n");
		
	}

	public int setup(String arg0, ImagePlus arg1) {
		if(arg1!=null){
			img = arg1;
		}
		else{
			IJ.showMessage("No image!");
			return DONE;
		}
		
		//modify calibration 
		Calibration cal = img.getCalibration();
		cal.pixelWidth 	= cal.pixelHeight = cal.pixelDepth = 1;
		cal.setUnit("pixels");
		img.setCalibration(cal);
		
		GenericDialog gd = new GenericDialog("EXTRACT PATCHES");
		gd.addNumericField("patch  size", 10, 0);
		gd.addNumericField("start scale", 1, 0);
		gd.addNumericField("end   scale", 4, 0);
		gd.addNumericField("nr   scales", 4, 0);
		gd.showDialog();
		if (gd.wasCanceled()) return DONE;
		ptch_size 	= (int)	gd.getNextNumber();
		s_start 	= 		gd.getNextNumber();
		s_end		= 		gd.getNextNumber();
		s_number	= (int)	gd.getNextNumber();
		
		s = new double[s_number];
		for (int i = 0; i < s_number; i++) {
			s[i] = (i==0)? s_start : s_start+i*((s_end-s_start)/(s_number-1));
		}
		
		IJ.log("scales: "+Utils.arrayToString(s));
		
		return DOES_8G+NO_CHANGES;
	}

	public void mouseClicked(MouseEvent e) {
		
		int mouseX = 	img.getWindow().getCanvas().offScreenX(e.getX());
		int mouseY = 	img.getWindow().getCanvas().offScreenY(e.getY()); 
		System.out.format("mouse clicked at: x: %d y: %d \n", mouseX, mouseY);
		
		// extract direction
		
		
		// extract patch
		
		
	}

	private double[] extractDirection(){
//		for (int i = 0; i < s.length; i++) {
//			
//		}
		return new double[2];
	}
	
	public void mousePressed(MouseEvent e) {
		
	}

	public void mouseReleased(MouseEvent e) {
		
	}

	@Override
	public void mouseEntered(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void mouseExited(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}

}
