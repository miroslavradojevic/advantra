package advantra.critpoint;

import java.awt.Color;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

import advantra.feature.CircHaarFeat;
import advantra.feature.FilterSet;
import advantra.general.Sort;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.plugin.PlugIn;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

public class DemoFeatures implements PlugInFilter, MouseListener {

    // is meant to test the features at particular locations
    // defined with the mouse click


    ImagePlus img; // this is the main input image

    // filter set parameters
    int     angScale1, angScale2;
    double  ring1, ring2;
    int     nr_rings;
    int     p_radius;

    int[]   angScale;
    double[] rings;

    float[] vals;
    float[] rads;
    float[] angs;
    FilterSet fs;

	public void run(ImageProcessor arg0) {


        img.show();

        for (int i = 0; i < 5; i++) { img.getCanvas().zoomIn(0, 0); }

        img.getCanvas().addMouseListener(this);

		/*
		 * test profile features
		 */
        GenericDialog gd = new GenericDialog("Demo Features");

        gd.addChoice("alfa_1",      new String[]{"20", "40", "60"}, "20");
        gd.addChoice("alfa_2",      new String[]{"20", "40", "60"}, "40");

        gd.addNumericField("ring_1:",       0.4,    1); //ring1
        gd.addNumericField("ring_2:",       0.7,    1); //ring2
        gd.addNumericField("nr_rings:",     2, 	    0,  5, "");//nr_rings

        gd.addNumericField("patch_size: ",    3.0, 	1);

        gd.showDialog();
        if (gd.wasCanceled()) return;

        switch ((int) gd.getNextChoiceIndex()) {
            case 0: angScale1 = 20;
                break;
            case 1: angScale1 = 40;
                break;
            case 2: angScale1 = 60;
                break;
            default: angScale1 = 40;
                break;
        }
        System.out.println("angScale1 = "+angScale1);

        switch ((int) gd.getNextChoiceIndex()) {
            case 0: angScale2 = 20;
                break;
            case 1: angScale2 = 40;
                break;
            case 2: angScale2 = 60;
                break;
            default: angScale2 = 40;
                break;
        }
        System.out.println("angScale2 = "+angScale2);

        ring1 =  gd.getNextNumber();            System.out.println("ring1 = "+ring1);
        ring2 =  gd.getNextNumber();            System.out.println("ring2 = "+ring2);
        nr_rings = (int) gd.getNextNumber();    System.out.println("nr_rings = "+nr_rings);
        p_radius =  (int)gd.getNextNumber();             System.out.println("patch_radius = "+p_radius);


        angScale = new int[(angScale2-angScale1)/20+1];
        int cnt = 0;
        for (int i = angScale1; i <= angScale2; i+=20){
            angScale[cnt] = i;
            System.out.println("angScale["+cnt+"] = "+angScale[cnt]);
            cnt++;
        }

        rings = new double[nr_rings];
        for (int i = 0; i <nr_rings; i++){
            rings[i] = (i==0)? ring1 : ring1+i*((ring2-ring1)/(nr_rings-1)) ;
            System.out.println("rings["+i+"] = "+rings[i]);
        }

        // create filter set
        fs = new FilterSet(angScale, rings, new double[1]);
        int nrFilters = fs.circConfs.size()+fs.radlConfs.size();
        System.out.println(nrFilters+" filters formed!");
        /*
        show them
        */
        ImagePlus all_feats = new ImagePlus("FEATURES", fs.plot(65));
        all_feats.setTitle("All_Features");
        all_feats.show();

        int toAlloc = 0;
        for (int x = -p_radius; x <= p_radius; x++){
            for (int y = -p_radius; y <= p_radius; y++){
                if(x*x+y*y<=p_radius*p_radius)
                    toAlloc++;
            }

        }

        System.out.println("to allocate: "+toAlloc+" with "+p_radius+" patch radius");

        vals = new float[toAlloc];
        angs = new float[toAlloc];
        rads = new float[toAlloc];

//        gd.addNumericField("min scale", a1, 0, 5, "deg");
//        gd.addNumericField("max scale", a2, 0, 5, "deg");
//        gd.addNumericField("ang. step", as, 0, 5, "deg");
//
//        gd.addMessage("radial scale:");
//        gd.addNumericField("min scale", r1, 1);
//        gd.addNumericField("max scale", r2, 1);
//        gd.addNumericField("ang. step", rn, 0, 5, "#");
//

//
//        a1 	= 		    (int)gd.getNextNumber();
//        a2	= 		    (int)gd.getNextNumber();
//        as	= 		    (int)gd.getNextNumber();
//
//        r1 	= 		    gd.getNextNumber();
//        r2	= 		    gd.getNextNumber();
//        rn	= 		    (int)gd.getNextNumber();
//
//        int aSclNr = 0;
//        for (int i = a1; i <= a2; i+=as) aSclNr++;
//
//		int[] aScl = new int[aSclNr];
//        for (int i = 0; i < aSclNr; i++) aScl[i] = a1+i*as;
//
//        double[] rScl = new double[rn];
//        for (int i = 0; i < rn; i ++){
//            rScl[i] = (i==0)? r1 : r1+i*((r2-r1)/(rn-1));
//        }
//
////        for (int i = 0; i < aSclNr; i++) System.out.println(i + " : " + aScl[i]);
////        for (int i = 0; i < rn; i++) System.out.println(i + " : " + rScl[i]);
//
//		FilterSet fs = new FilterSet(aScl, new double[]{0.5, 1.0}, rScl);
//        fs.print();
//		new ImagePlus("features", fs.plot(101)).show();
//		System.out.println("total nr. configurations: "+ (fs.circConfs.size()+fs.radlConfs.size()));
		
	}

    public int setup(String arg0, ImagePlus arg1) {

        if(arg1!=null){
            img = arg1;
        }
        else{
            IJ.showMessage("Open image to 'demo features' first!");
            return DONE;
        }

        return DOES_8G+NO_CHANGES;

    }

    public void mouseClicked(MouseEvent e) {



        int atX = 	img.getWindow().getCanvas().offScreenX(e.getX());
        int atY = 	img.getWindow().getCanvas().offScreenY(e.getY());

        System.out.println("click "+atX+" , "+atY);

        profile(img, atX, atY, p_radius);

        int choose_ft;
        GenericDialog gd = new GenericDialog("Choose Feature");
        gd.addNumericField("choose feature:" ,	1, 0);
        gd.showDialog();
        if (gd.wasCanceled()) return;
        choose_ft = (int)       gd.getNextNumber();

        System.out.print("chosen ft: "+choose_ft);
        float[] ft_id = new float[fs.circConfs.get(choose_ft).nrRot];
        for (int i = 0; i < ft_id.length; i++) ft_id[i] = i;

        new ImagePlus("chosenfeature", plotProfile()).show();

//        for (int i = 0; i < vals.length; i++){
//            System.out.println(""+i+" -> "+vals[i]);
//        }

        new ImagePlus("chosen_ft"+choose_ft+"/"+(fs.circConfs.size()-1), fs.circConfs.get(choose_ft).plotAllRotations(65)).show();
        Plot p = new Plot("ft"+choose_ft+"/"+(fs.circConfs.size()-1), "rot#", "scores per rot.", ft_id, fs.circConfs.get(choose_ft).scoreAllRot(vals, angs, rads));
        p.show();


    }

    public void profile(
            ImagePlus img,
            int xin,
            int yin,
            int rin
    )
    {

        // will extract out vals[], rads[], and angs[]
        // those will be used to score() on filterSet

        int xs = xin-rin;
        int xe = xin+rin;

        int ys = yin-rin;
        int ye = yin+rin;

        int cnt = 0;
        for (int x = xs; x <= xe; x++){
            for (int y = ys; y <= ye; y++){

//                float d = (float) (Math.pow(, 2) + Math.pow(y-yin, 2));
                int d =  (x-xin)*(x-xin)+(y-yin)*(y-yin);

                if(d <= rin*rin){

                    vals[cnt] = img.getProcessor().getPixelValue(x, y);
                    rads[cnt] = (float) (Math.sqrt(d) / rin);
                    angs[cnt] = (float) (Math.atan2(y-yin, x-xin) + Math.PI);
                    angs[cnt] = (angs[cnt]>=(float)(2*Math.PI))? 0 : angs[cnt];
                    angs[cnt] = (angs[cnt]<0)? 0 : angs[cnt];
                    cnt++;

                }

            }
        }

        System.out.println("cnt: "+cnt);
    }


    public ImageStack plotProfile()
    {
        ImageStack viz = new ImageStack(400, 200);

        // find max for plotting
        float max_val = vals[0];
        for (int i = 1; i < vals.length; i++){
            if (vals[i] > max_val) max_val = vals[i];
        }

        Plot p = new Plot("circular_profile", "angle[rad]", "value");
        p.setLimits(0, 2*Math.PI, 0, max_val);
        p.setSize(400, 200);
        p.addPoints(angs, vals, Plot.BOX);
        viz.addSlice("circular_profile", p.getProcessor());

        Plot p1 = new Plot("radial_profile", "radius", "value");
        p1.setLimits(0, 1, 0, max_val);
        p1.setSize(400, 200);
        p1.addPoints(rads, vals, Plot.CIRCLE);
        viz.addSlice("radial_profile", p1.getProcessor());

        // check circular score
        float[] c_idx = new float[fs.circConfs.size()];
        float[] c_sco = new float[fs.circConfs.size()];
        for (int i = 0; i < fs.circConfs.size(); i++){
            c_idx[i] = i;
            c_sco[i] = fs.score[i];
        }
        Plot p3 = new Plot("filtering_score", "circular_configuration", "score", c_idx, c_sco);
        p3.setSize(400, 200);
        viz.addSlice("scores", p3.getProcessor());

        // check radial score
        float[] r_idx = new float[fs.radlConfs.size()];
        float[] r_sco = new float[fs.radlConfs.size()];
        for (int i = 0; i < fs.radlConfs.size(); i++){
            r_idx[i] = i;
            r_sco[i] = fs.score[fs.circConfs.size()+i];
        }
//        Plot p4 = new Plot("filtering_score", "radial_configuration", "score", r_idx, r_sco);
//        p4.setSize(400, 200);
//        viz.addSlice("scores", p4.getProcessor());

        return viz;
    }

    @Override
    public void mousePressed(MouseEvent e) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void mouseReleased(MouseEvent e) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void mouseEntered(MouseEvent e) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void mouseExited(MouseEvent e) {
        //To change body of implemented methods use File | Settings | File Templates.
    }
}
