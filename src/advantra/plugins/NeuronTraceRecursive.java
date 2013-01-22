package advantra.plugins;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.File;
import java.text.SimpleDateFormat;
import java.util.Calendar;

import advantra.general.CreateDirectory;
import advantra.trace.NeuronTrace;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.ImageCanvas;
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
		MAX_BRANCHES = 4;
		return DOES_8G+NO_CHANGES;
	}
	
	public void run(ImageProcessor ip) {
		img.getWindow().getCanvas().addMouseListener(this);
		System.out.println("click on the point to start tracing from there...");
	}
	
	public void mouseClicked(MouseEvent e) {
		
		int mouseY = 	img.getWindow().getCanvas().offScreenX(e.getX());
		int mouseX = 	img.getWindow().getCanvas().offScreenY(e.getY()); 
		int mouseZ =  	img.getCurrentSlice()-1;
		
		trace(mouseX, mouseY, mouseZ); // trace from the selected point
		
	}
	
	public void mousePressed(MouseEvent e) {}
	public void mouseReleased(MouseEvent e) {}
	public void mouseEntered(MouseEvent e) {}
	public void mouseExited(MouseEvent e) {}
	
	private void trace(double mouseX, double mouseY, double mouseZ){

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
