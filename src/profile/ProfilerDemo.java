package profile;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.ImageCanvas;
import ij.gui.Plot;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 7/7/13
 * Time: 8:13 PM
 */
public class ProfilerDemo implements PlugInFilter, MouseListener {

	ImagePlus 	inimg;
	String 		inimgPath;
	ImagePlus   inmask;

    // visualization
    ImagePlus       vizProfileImage;
    ImageStack      vizProfileStack;

    // to store profiles, list of values per location
    ArrayList<ArrayList<float[]>>   profilesPerLocation;
    ArrayList<ArrayList<Integer>>   angResDegPerLocation;
    ArrayList<ArrayList<String>>    profileNamePerLocation;

    int CPU_NR;

	public int setup(String s, ImagePlus imagePlus) {

		if(imagePlus==null) {
			IJ.showMessage("needs image opened"); return DONE; }
		inimg = Tools.convertToFloatImage(imagePlus);
		inimg.setTitle("inimg");
		inimgPath = imagePlus.getOriginalFileInfo().directory+imagePlus.getOriginalFileInfo().fileName;
		return DOES_8G+DOES_32+NO_CHANGES;

	}

	public void run(ImageProcessor imageProcessor) {

		String inmaskPath = Tools.removeExtension(inimgPath)+".mask";

		if (new File(inmaskPath).exists()) {
			inmask = new ImagePlus(inmaskPath);
			inmask.setTitle("inmask");
		}
		else {
			byte[] array = new byte[inimg.getWidth() * inimg.getHeight()];
			for (int a = 0; a < array.length; a++) array[a] = (byte) 255;
			inmask = new ImagePlus("inmask", new ByteProcessor(inimg.getWidth(), inimg.getHeight(), array));
		}

        CPU_NR = 5;

		Profiler.loadTemplate(inimg.getProcessor(), (ByteProcessor) inmask.getProcessor());

        int totalLocations = Profiler.locations.length;
        profilesPerLocation     = new ArrayList<ArrayList<float[]>>(totalLocations);
        angResDegPerLocation    = new ArrayList<ArrayList<Integer>>(totalLocations);
        profileNamePerLocation  = new ArrayList<ArrayList<String>>(totalLocations);

        //loop parameters
        for (double neuronDiam = 3; neuronDiam<=3; neuronDiam++) {
            for (double scale=2; scale<=3; scale+=1) {

                Profiler.loadParams(neuronDiam, scale);
                IJ.log("calculating profiles... neuronDiam="+neuronDiam+", scale="+scale);
                long t1 = System.currentTimeMillis();
                int totalProfiles = Profiler.offsets.size();
                Profiler ms_jobs[] = new Profiler[CPU_NR];
                for (int i = 0; i < ms_jobs.length; i++) {
                    ms_jobs[i] = new Profiler(i*totalProfiles/CPU_NR,  (i+1)*totalProfiles/CPU_NR);
                    ms_jobs[i].start();
                }
                for (int i = 0; i < ms_jobs.length; i++) {
                    try {
                        ms_jobs[i].join();
                    } catch (InterruptedException e) {
                        e.printStackTrace();
                    }
                }
                long t2 = System.currentTimeMillis();
                IJ.log("done extracting profiles "+((t2-t1)/1000f)+" sec.");

                updateList();

            }
        }

        // ms, use profilesPerLocation
        IJ.log("calculating local peaks... ");
        long t1 = System.currentTimeMillis();

        Analyzer.loadProfiles(profilesPerLocation);
        Analyzer ms_jobs[] = new Analyzer[CPU_NR];
        for (int i = 0; i < ms_jobs.length; i++) {
            ms_jobs[i] = new Analyzer(i*totalLocations/CPU_NR,  (i+1)*totalLocations/CPU_NR);
            ms_jobs[i].start();
        }
        for (int i = 0; i < ms_jobs.length; i++) {
            try {
                ms_jobs[i].join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }

        long t2 = System.currentTimeMillis();
        IJ.log("done extracting peaks "+((t2-t1)/1000f)+" sec.");

        inimg.show();
        inimg.getCanvas().addMouseListener(this);

        vizProfileImage = new ImagePlus();
        IJ.setTool("hand");

        inmask.show();
		IJ.selectWindow("inimg");
		IJ.run("Add Image...", "image=inmask x="+0+" y="+0+" opacity=20");
		inmask.close();

	}

    private void updateList() // uses Profiler
    {
        // store it in the list for each loc
        for (int loopLocations=0; loopLocations<Profiler.locations.length; loopLocations++) {

            // profileToAdd from Profiler.profiles
            float[] profileToAdd = new float[Profiler.profiles[loopLocations].length];
            for (int k=0; k<profileToAdd.length; k++) {
                profileToAdd[k] = Profiler.profiles[loopLocations][k];
            }

            // stringToAdd from Profiler.neuronDiam, Profiler.scale, Profiler.resolDeg
            String stringToAdd = "nD_"+Profiler.neuronDiam+"_s_"+Profiler.scale+"_a_"+Profiler.resolDeg;

            // intToAdd from Profiler.resolDeg
            int intToAdd = Profiler.resolDeg;

            if (profilesPerLocation.size()<Profiler.locations.length) {

                ArrayList<float[]> A = new ArrayList<float[]>();
                A.add(profileToAdd);
                profilesPerLocation.add(A);

                ArrayList<String> B = new ArrayList<String>();
                B.add(stringToAdd);
                profileNamePerLocation.add(B);

                ArrayList<Integer> C = new ArrayList<Integer>();
                C.add(intToAdd);
                angResDegPerLocation.add(C);

            }
            else {
                profilesPerLocation.get(loopLocations).add(profileToAdd);
                profileNamePerLocation.get(loopLocations).add(stringToAdd);
                angResDegPerLocation.get(loopLocations).add(intToAdd);
            }

        }
    }

    public void mouseClicked(MouseEvent e) {

        ImageCanvas srcCanv = (ImageCanvas) e.getSource();
        int atX = 	srcCanv.offScreenX(e.getX());
        int atY = 	srcCanv.offScreenY(e.getY());

        vizProfileStack = new ImageStack(600, 300);

        String fileName = "plotProfiles.r";
        // empty the file
        PrintWriter writer = null;
        try {
            writer = new PrintWriter(fileName);
        } catch (FileNotFoundException ex) {
            //ex.printStackTrace();
        }
        writer.print("");
        writer.close();

        for (int i=0; i<Profiler.locations.length; i++) {
            if (Profiler.locations[i][0]==atX && Profiler.locations[i][1]==atY) {

                // pick all the profiles from the list for that location

                int nrProfiles = profilesPerLocation.get(i).size();

                for (int g=0; g<nrProfiles; g++) {

                    int currProfileLen = profilesPerLocation.get(i).get(g).length;

                    /*
                    plot ij
                     */
                    float[] angIdx = new float[currProfileLen];

                    for (int q=0; q<angIdx.length; q++) {
                        angIdx[q] = q * angResDegPerLocation.get(i).get(g);
                    }

                    Plot p = new Plot("", "orient.[deg]", profileNamePerLocation.get(i).get(g), angIdx, profilesPerLocation.get(i).get(g));
                    p.setSize(600, 300);
                    vizProfileStack.addSlice("", p.getProcessor());

                    /*
                    export R
                     */

                    String printProfile = "";
                    String printAngle = "";
                    String xName = "ang_"+profileNamePerLocation.get(i).get(g);
                    String yName = profileNamePerLocation.get(i).get(g);

                    printProfile    += yName+" <- c(";
                    printAngle      += xName+" <- c(";

                    for (int i1=0; i1<currProfileLen; i1++) {

                        printProfile+=profilesPerLocation.get(i).get(g)[i1]+"";
                        printAngle+=angIdx[i1]+"";

                        if(i1<profilesPerLocation.get(i).get(g).length-1) {
                            printProfile+=", ";
                            printAngle+=", ";
                        }
                    }

                    printProfile+=")\n";
                    printAngle+=")\n";

                    try {
                        PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(fileName, true)));
                        out.println(printAngle);
                        out.println(printProfile);
                        if (g==0)
                            out.println("plot("+xName+","+yName+",type=\"l\")\n grid()");
                        else
                            out.println("lines("+xName+","+yName+",type=\"l\")");
                        out.close();

                    } catch (IOException e1) {}

                }

            }
        }

        if (vizProfileStack.getSize()>0) {
            vizProfileImage.setStack(vizProfileStack);
            vizProfileImage.updateAndDraw();
            vizProfileImage.setTitle("profiles");
            vizProfileImage.show();
        }

    }

    public void mousePressed(MouseEvent e) {}

    public void mouseReleased(MouseEvent e) {}

    public void mouseEntered(MouseEvent e) {}

    public void mouseExited(MouseEvent e) {}
}
