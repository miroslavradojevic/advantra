package advantra.plugins;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.File;
import java.text.SimpleDateFormat;
import java.util.Calendar;

import flanagan.analysis.Stat;

import advantra.general.CreateDirectory;
import advantra.tools.OtsuBinarisation;
import advantra.trace.NeuronTrace;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.ImageCanvas;
import ij.measure.Calibration;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

public class NeuronTraceRecursive implements PlugInFilter, MouseListener {

	ImagePlus img;
	ImageCanvas canvas;
	int MAX_BRANCHES;
	
	public int setup (String arg, ImagePlus img){
		this.img = img;
		if(img==null){ //.getStack().getSize()<=1
			IJ.error("The image was NULL... Select new image...");
			return DONE;
		}
		if(img.getStack().getSize()<=1) {
			IJ.error("The image was not a stack... Select new image...");
			return DONE;
		}
		if(img.getType()!=ImagePlus.GRAY8) {
			IJ.error("The image was not GRAY8... Select new image...");
			return DONE;
		}
		
		//modify calibration 
		Calibration cal = img.getCalibration();
		cal.pixelWidth 	= 1;
		cal.pixelHeight	= 1;
		cal.pixelDepth	= 1;
		cal.setUnit("pixels");
		img.setCalibration(cal);
		
		MAX_BRANCHES = 4;
		return DOES_8G+NO_CHANGES;
	}
	
	public void run(ImageProcessor ip) {
		img.getWindow().getCanvas().addMouseListener(this);
		System.out.println("click on the point to start tracing from there...");
		IJ.log("click on the point to start tracing from there...");
	}
	
	public void mouseClicked(MouseEvent e) {
		
		int mouseX = 	img.getWindow().getCanvas().offScreenX(e.getX());
		int mouseY = 	img.getWindow().getCanvas().offScreenY(e.getY()); 
		int mouseZ =  	img.getCurrentSlice()-1;
		
		double value = img.getStack().getVoxel(mouseX, mouseY, mouseZ);
		System.out.format("value at (%d, %d, %d) is %f \n", mouseX, mouseY, mouseZ, value);
		
		// check if the point was roughly on neuron
		OtsuBinarisation otsu = new OtsuBinarisation(img); // binarise
		ImagePlus out = otsu.run();
		
		value = out.getStack().getVoxel(mouseX, mouseY, mouseZ);
		System.out.format("binary value at (%d, %d, %d) is %s \n", mouseX, mouseY, mouseZ, (value==255)?"ON":"OFF");
		
		// median filtering in the neighborhood extract 3x3x3 neighborhood around and see if it's still "ON"
		double[] neighborhood = new double[27];
		int idx = 0;
		int r, c, l;
		for (int i = -1; i <= 1; i++) {
			for (int j = -1; j <= 1; j++) {
				for (int k = -1; k <= 1; k++) {
					r = mouseY+i;	c = mouseX+j;	l = mouseZ+k;
					if(r>=0 && r<img.getHeight() && c>=0 && c<img.getWidth() && l>=0 && l<img.getStackSize()){
						neighborhood[idx] = out.getStack().getVoxel(c, r, l); 
						idx++;
					}
				}
			}
		}
		
		if((int)Math.round(Stat.median(neighborhood))==255){
			IJ.log("start trace from "+mouseX+", "+mouseY+", "+mouseZ);
		}
		else{
			IJ.log("couldn't start the trace from this point... click again...");
			return;
		}
		
		IJ.log("trace()...");
		trace(mouseY, mouseX, mouseZ); // trace from the selected point mouseX==col, mouseY==row
		
	}
	
	public void mousePressed(MouseEvent e) {}
	public void mouseReleased(MouseEvent e) {}
	public void mouseEntered(MouseEvent e) {}
	public void mouseExited(MouseEvent e) {}
	
	private void trace(double mouseX, double mouseY, double mouseZ){ // trace(row, col, lay)

		System.out.format("starting from: %f, %f, %f ...\n", mouseX, mouseY, mouseZ);

		String log_dir_name 		= 
				System.getProperty("user.home")																	+File.separator+
				"TraceLog_"+(new SimpleDateFormat("dd-MM-yyyy-HH-mm-ss")).format(Calendar.getInstance().getTime())	+File.separator;
		
		CreateDirectory.createOneDir(log_dir_name);
		
		NeuronTrace neuron_tr = new NeuronTrace(img, MAX_BRANCHES); 
		
		long start_time = System.currentTimeMillis();
		
		neuron_tr.trace(mouseX, mouseY, mouseZ);
		
		System.out.format("total elasped time: %f sec.\n", (double)(System.currentTimeMillis()-start_time)/1000 );
		
		neuron_tr.export_swc(log_dir_name+"recon.swc");
		
	}
	
}