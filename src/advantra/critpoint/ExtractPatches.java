package advantra.critpoint;

import java.awt.Button;
import java.awt.Color;
import java.awt.Panel;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.File;
import java.util.Vector;

import javax.swing.JFileChooser;

import advantra.feature.MyHessian;
import advantra.file.AnalyzeCSV;
import advantra.processing.IntensityCalc;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.Line;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import imagescience.image.Axes;
import imagescience.image.Coordinates;
import imagescience.image.Dimensions;
import imagescience.image.FloatImage;
import imagescience.image.Image;

public class ExtractPatches implements PlugInFilter, MouseListener, ActionListener {

	ImagePlus 	img;
	Image 		L1scales, L2scales, v1scales, v2scales; // each will be stack with each layer corresponding to a scale
	int 		ptch_size;
	double 		s_start;
	double 		s_end;
	int 		s_number;
	double[]	s;
	Button 		loadPoints;
	String 		load_points_file;
	
	ImageStack 	trn;
	ImagePlus	viz;
	
	private int half_diag;
	
	public void 	run(ImageProcessor arg0) {
		
		//modify calibration of the input image
		Calibration cal = img.getCalibration();
		cal.pixelWidth 	= cal.pixelHeight = cal.pixelDepth = 1;
		cal.setUnit("pixel");
		img.setCalibration(cal);
		
		s_start  = Prefs.get("advantra.critpoint.start_scale", 0.5);
		s_end    = Prefs.get("advantra.critpoint.end_scale", 5.0);
		s_number = Prefs.getInt("advantra.critpoint.nr_scales", 10);
		
		GenericDialog gd = new GenericDialog("Extract Patches");
		
		gd.addNumericField("patch  size", 15, 0);
		gd.addNumericField("start scale", s_start, 1);
		gd.addNumericField("end   scale", s_end, 1);
		gd.addNumericField("nr   scales", s_number, 0);
		
		Panel buttons_panel = new Panel();
		loadPoints = new Button("Load points");
		loadPoints.addActionListener(this);
		
		buttons_panel.add(loadPoints);
		gd.add(buttons_panel);
		
		gd.showDialog();
		if (gd.wasCanceled()) return;
		ptch_size 	= (int)	gd.getNextNumber();
		s_start 	= 		gd.getNextNumber();
		s_end		= 		gd.getNextNumber();
		s_number	= (int)	gd.getNextNumber();
		
		Prefs.set("advantra.critpoint.start_scale", s_start);
		Prefs.set("advantra.critpoint.end_scale", s_end);
		Prefs.set("advantra.critpoint.nr_scales", s_number);
		
		half_diag = (int)Math.ceil(ptch_size/Math.sqrt(2));
		
		s = new double[s_number];
		for (int i = 0; i < s_number; i++) {
			s[i] = (i==0)? s_start : s_start+i*((s_end-s_start)/(s_number-1));
		}

		trn = new ImageStack(ptch_size, ptch_size); 
		
		viz = new ImagePlus();
		
		// extract hessian eigen values and eigen vectors for different scales
		Image im = new FloatImage(Image.wrap(img));
		
		Dimensions dims 		= im.dimensions();
		
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
			
			System.out.println("extracting hessian for scale "+IJ.d2s(s[i]));
			
			Vector<Image> hess = my_hess.eigs(im.duplicate(), s[i], true);
			
			// assign values to layers of L1, L2, v1, v2
			Image L2 	= hess.get(0); L2.axes(Axes.X);
			
			Image L1 	= hess.get(1); L1.axes(Axes.X);
			
			Image V11 	= hess.get(2); V11.axes(Axes.X);
			
			Image V12 	= hess.get(3); V12.axes(Axes.X);
			
			
			for (coords.y=0; coords.y<dims.y; ++coords.y) {
				
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
		
//		L1scales.imageplus().show();
//		L2scales.imageplus().show();
//		v1scales.imageplus().show();
//		v2scales.imageplus().show();
		
		img.getWindow().getCanvas().addMouseListener(this);
		
		IJ.showMessage("Hessian calculations done.\nClick now!");
		
		// read locations if button selected the locations
		
		if(load_points_file!=null){
			
			AnalyzeCSV analyze_csv = new AnalyzeCSV(load_points_file);
			System.out.println("read "+analyze_csv.getLinesNr()+" locations");
			double[][] poss = analyze_csv.readLn(2);// col, row
			System.out.println("loaded "+poss.length+" locations");
			IntensityCalc im_calc = new IntensityCalc(img.getStack());
			
			for (int i = 0; i < poss.length; i++) {
				
				int x = (int)poss[i][0];
				int y = (int)poss[i][1];
				addPatchToStack(x, y, im_calc);
				
			}
			
			viz.setStack(trn);
			viz.setTitle("train");
			viz.show();
			viz.getCanvas().zoomIn(0,0);
			
		}
		
	}

	public int 		setup(String arg0, ImagePlus arg1) {
		if(arg1!=null){
			img = arg1;
		}
		else{
			IJ.showMessage("Open image first!");
			return DONE;
		}
		
		return DOES_8G+NO_CHANGES;
	}

	public void 		mouseClicked(MouseEvent e) {
		
		if(e.getClickCount()==2){
			// save the stack and get out
			FileSaver fs = new FileSaver(viz);
			
			GenericDialog gd = new GenericDialog("Output stack file name");
			gd.addStringField("TIF file name: ", "...");
			gd.showDialog();
			if (gd.wasCanceled()) return;
			String name = gd.getNextString();
			
			ptch_size 	= (int)	gd.getNextNumber();
			fs.saveAsTiffStack(System.getProperty("user.home")+File.separator+name+".tif");
			System.out.println(System.getProperty("user.home")+File.separator+name+".tif"+"   SAVED.");
			return;
		}
		
		int mouseX = 	img.getWindow().getCanvas().offScreenX(e.getX());
		int mouseY = 	img.getWindow().getCanvas().offScreenY(e.getY());
		
		// check whether patch can be extracted due to borders
		if(
				mouseX<=half_diag || mouseX>=img.getWidth()-half_diag ||
				mouseY<=half_diag || mouseX>=img.getHeight()-half_diag
				){
			System.out.println("won't extract on this point");
			return;
		}
		
		// extract direction (only at selected point)
		double[] vec = extractDirection(mouseX, mouseY);
		
		double theta = Math.atan(vec[1]/vec[0]);
		
		// extract locs
		float[][] locs = extractLocs(mouseX, mouseY, theta);
		
		// overlay locs - just for vizualization
		Roi roi = new PointRoi(locs[0], locs[1], ptch_size*ptch_size);
		Roi roi_line = new Line(mouseX, mouseY, (mouseX+vec[0]*ptch_size), (mouseY+vec[1]*ptch_size));
		roi.setStrokeColor(Color.red);
		Overlay ovly = new Overlay(roi);
		ovly.add(roi_line);
		img.setOverlay(ovly);

		// extract patch
		IntensityCalc im_calc = new IntensityCalc(img.getStack());
		byte[] ptch = extractPatch(mouseX, mouseY, theta, im_calc);
		
		trn.addSlice(new ByteProcessor(ptch_size, ptch_size, ptch));
		
		viz.setStack(trn);
		viz.setTitle("train");
		viz.show();
		viz.getCanvas().zoomIn(0,0);
		
		System.out.println("layer: "+trn.getHeight()+" x "+trn.getWidth()+" , curr. stack size: "+trn.getSize());
		
	}
	
	private void addPatchToStack(int x, int y, IntensityCalc im_calc){
		
		// check whether patch can be extracted due to borders
		if(x<=half_diag || x>=img.getWidth()-half_diag ||
			y<=half_diag || y>=img.getHeight()-half_diag){
			System.out.println("won't extract on this point");
			return;
		}
		
		// extract direction (only at selected point)
		double[] vec = extractDirection(x, y);
		double theta = Math.atan(vec[1]/vec[0]);
		byte[] ptch = extractPatch(x, y, theta, im_calc);
		trn.addSlice(new ByteProcessor(ptch_size, ptch_size, ptch));
		
	}

	private double[] 	extractDirection(int x, int y){

		double[] v = new double[2];
		
		Coordinates coord = new Coordinates();
		coord.x = x;	coord.y = y;
		
		L1scales.axes(Axes.Z);
		L2scales.axes(Axes.Z);
		
		double[] aL1 = new double[s.length];
		double[] aL2 = new double[s.length];
		
		L1scales.get(coord, aL1);
		L2scales.get(coord, aL2);
		
		double min_ratio = (aL2[0]<0)?Math.abs(aL1[0]/aL2[0]):Double.MAX_VALUE;
		int scale_idx = 0;
		
		for (int i = 0; i < s.length; i++) {
			if(aL2[i]<0){
				double ratio = Math.abs(aL1[i]/aL2[i]);
				if(ratio<min_ratio){
					min_ratio = ratio;
					scale_idx = i;
				}
			}
		}
		
		coord.z = scale_idx;
		v[0] = v1scales.get(coord);
		v[1] = v2scales.get(coord);
		
		return v;
	}
	
	private byte[] 		extractPatch(int x, int y, double theta, IntensityCalc img_calc){
		
		byte[] 	patch 		= new byte		[ptch_size*ptch_size];
		float[][] 	locs 	= new float		[2][ptch_size*ptch_size];
		
		int N1 = ptch_size/2;
		
		float x_beg = x-N1;
		float y_beg = y-N1;
		float x_cur, y_cur;
		
		int cnt = 0;
		for (int i = 0; i < ptch_size; i++) {
			
			x_cur = x_beg + i;
			
			for (int j = 0; j < ptch_size; j++) {
				
				y_cur = y_beg + j;
				
				locs[0][cnt] = x_cur-x;
				locs[1][cnt] = y_cur-y;
				cnt++;
				
			}
		}
		
		// rotate
		for (int i = 0; i < locs[0].length; i++) {
			float x_rot = locs[0][i]*(float)Math.cos(theta)-locs[1][i]*(float)Math.sin(theta);
			float y_rot = locs[0][i]*(float)Math.sin(theta)+locs[1][i]*(float)Math.cos(theta);
			locs[0][i] = x_rot;
			locs[1][i] = y_rot;
		}
		
		// return back
		for (int i = 0; i < locs[0].length; i++) {
			locs[0][i] += x;
			locs[1][i] += y;
		}
		
		cnt = 0;
		for (int i = 0; i < ptch_size; i++) {
			for (int j = 0; j < ptch_size; j++) {
				patch[cnt] = (byte)Math.round(img_calc.interpolateAt(locs[1][cnt], locs[0][cnt]));
				cnt++;
			}
		}
		
		return patch;
		
	}
	
	private float[][] 	extractLocs(int x, int y, double theta){
		
		float[][] locs = new float[2][ptch_size*ptch_size]; 
		
		int N1 = ptch_size/2;
		
		float x_beg = x-N1;
		float y_beg = y-N1;
		float x_cur, y_cur;
		
		int cnt = 0;
		
		for (int i = 0; i < ptch_size; i++) {
			
			x_cur = x_beg + i;
			
			for (int j = 0; j < ptch_size; j++) {
				
				y_cur = y_beg + j;
				
				locs[0][cnt] = x_cur-x;
				locs[1][cnt] = y_cur-y;
				cnt++;
				
			}
		}
		
		// rotate
		for (int i = 0; i < locs[0].length; i++) {
			float x_rot = locs[0][i]*(float)Math.cos(theta)-locs[1][i]*(float)Math.sin(theta);
			float y_rot = locs[0][i]*(float)Math.sin(theta)+locs[1][i]*(float)Math.cos(theta);
			locs[0][i] = x_rot;
			locs[1][i] = y_rot;
		}
		
		// return back
		for (int i = 0; i < locs[0].length; i++) {
			locs[0][i] += x;
			locs[1][i] += y;
		}
		
		return locs;
	}
	
	public void mousePressed(MouseEvent e) {
		
	}

	public void mouseReleased(MouseEvent e) {
		
	}

	public void mouseEntered(MouseEvent e) {
		
	}

	public void mouseExited(MouseEvent e) {
		
	}

	public void actionPerformed(ActionEvent e) {
		if(e.getSource()==loadPoints){
			JFileChooser fc = new JFileChooser();
			int returnVal = fc.showOpenDialog(null);
	        if (returnVal != JFileChooser.APPROVE_OPTION) {
	            return;
	        }
	        
	        load_points_file = fc.getSelectedFile().getAbsolutePath();
	        IJ.showStatus("Opened " + load_points_file);
		}
	}

}
