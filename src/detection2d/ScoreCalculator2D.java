package detection2d;

/**
 * Created by miroslav on 1/9/14.
 * takes the robustly estimated peaks
 * loops the locations' peak skeletons extracted using PeakAnalyzer2D (the usage of that class draws th usage of the other classes as well)
 * detector checks if the robust, consistent peak recursion evolved till the proposed end, along at least three branches
 * the "spatial consistency" is checked to account for the peak robustness (it's neighbours have to point to the spatially close location)
 * to avoid having some outlier peaks included but have groups of peaks pointing together, agreeing on the same information
 * gives out the filtered version of delin2 where only the stable ones are taken
 */

public class ScoreCalculator2D extends Thread {

    private int begN, endN;

    // VARIABLES (mainly used as a table)
    public static int[][] 	    i2xy;                   // selected locations
    public static int[][]     	xy2i;                   // need for recursion

    // INPUT: list of extracted peaks
    public static int[][][]     peaks_xy;             	// N x 4(max. threads) x 2   every FG location with 4 selected peaks in XY format

    // PARAMETERS


    // OUTPUT



}
