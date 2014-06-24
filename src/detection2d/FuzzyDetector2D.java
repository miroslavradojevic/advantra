package detection2d;

/**
 * Created by miroslav on 22-6-14.
 */
public class FuzzyDetector2D extends Thread {

    // loop through the features extracted at each location and makes the fuzzy detection reading from the extracted features

    private int begN, endN;

    // INPUT
    public static float[][]     ncc;                    // N(foreground locs.) x 4(max. threads) ncc of the patch (from Ncc2D)
    public static float[][]		lhood; 					// N(foreground locs.) x 4(max. threads) normalized 0-1 likelihoods (from PeakExtractor2D)
    public static float[][]		smthness;				// N(foreground locs.) x 4(max. threads) smoothness score (from Delineator2D)

    private static  int nr_points = Integer.MIN_VALUE;

    private static float        ncc_high            = Float.NaN;
    private static float        ncc_low             = Float.NaN;
    private static float        likelihood_high     = Float.NaN;
    private static float        likelihood_low      = Float.NaN;
    private static float        smoothness_high     = Float.NaN;
    private static float        smoothness_low      = Float.NaN;
    private static float        output_sigma        = Float.NaN;

    // OUTPUT
    public static float[]   	endpoint_score;         // N(foreground locs.) x 1(endpoint) given by fuzzy logic
    public static float[]   	junction_score;         // N(foreground locs.) x 1(bifpoint) given by fuzzy logic

    public FuzzyDetector2D(int n0, int n1)
    {
        this.begN = n0;
        this.endN = n1;
    }

    public static void loadTemplate(
            float[][] _ncc,
            float[][] _lhood,
            float[][] _smthness,
            float       _ncc_high,
            float       _ncc_low,
            float       _likelihood_high,
            float       _likelihood_low,
            float       _smoothness_high,
            float       _smoothness_low,
            float       _output_sigma
    )
    {

        ncc = _ncc;
        lhood = _lhood;
        smthness = _smthness;

        nr_points = ncc.length; // total points to process (other features should have the same length)

        ncc_high = _ncc_high;
        ncc_low = _ncc_low;

        likelihood_high = _likelihood_high;
        likelihood_low = _likelihood_low;

        smoothness_high = _smoothness_high;
        smoothness_low = _smoothness_low;

        output_sigma = _output_sigma;

        // allocate outputs
        endpoint_score = new float[nr_points];
        junction_score = new float[nr_points];

    }

    public void run()
    {
        // aux
        Fuzzy2D fls = new Fuzzy2D(// allocate one for each thread, important not to make it static so that all the calculations are independent for every location
                ncc_high, // ncc
                ncc_low,

                likelihood_high, // lhood
                likelihood_low,

                smoothness_high, // smoothness
                smoothness_low,

                output_sigma  // std output membership functions - defines separation
        );
        fls.verbose = false; // make it sure that there is no logging...
        float[] tmp = new float[3];// will take end,non,jun score using fls, will be sequaentially used
        float ncc_1=0, lhood_1=0, smthness_1=0,
                ncc_2=0, lhood_2=0, smthness_2=0,
                    ncc_3=0, lhood_3=0, smthness_3=0,
                        ncc_4=0, lhood_4=0, smthness_4=0; // allocate potential inputs (up to 4 scores for three categories)
        int cnt_valid;
        // aux

        for (int locationIdx = begN; locationIdx < endN; locationIdx++) {

            if (ncc[locationIdx]!=null && lhood[locationIdx]!=null && smthness[locationIdx]!=null) {

                // take those that are !=NaN and supply them to the decision scheme
                cnt_valid = 0;

                for (int i = 0; i < 4; i++) {

                    boolean score_exists =
                            Float.isNaN(ncc[locationIdx][i]) &&
                            Float.isNaN(lhood[locationIdx][i]) &&
                            Float.isNaN(smthness[locationIdx][i]);

                    if (score_exists) {

                        cnt_valid++;
                        // assign to dummy variables
                        switch (cnt_valid) {
                            case 1:
                                ncc_1       = ncc[locationIdx][i];
                                lhood_1     = lhood[locationIdx][i];
                                smthness_1  = smthness[locationIdx][i];
                                break;
                            case 2:
                                ncc_2       = ncc[locationIdx][i];
                                lhood_2     = lhood[locationIdx][i];
                                smthness_2  = smthness[locationIdx][i];
                                break;
                            case 3:
                                ncc_3       = ncc[locationIdx][i];
                                lhood_3     = lhood[locationIdx][i];
                                smthness_3  = smthness[locationIdx][i];
                            case 4:
                                ncc_4       = ncc[locationIdx][i];
                                lhood_4     = lhood[locationIdx][i];
                                smthness_4  = smthness[locationIdx][i];
                        }

                    }

                }

                switch (cnt_valid) {
                    case 0: break; // all were NaN do nothing, they're already initialized to zero
                    case 1: fls.critpointScore(ncc_1, lhood_1, smthness_1, tmp); break;
                    case 2: fls.critpointScore(ncc_1, lhood_1, smthness_1, ncc_2, lhood_2, smthness_2, tmp); break;
                    case 3: fls.critpointScore(ncc_1, lhood_1, smthness_1, ncc_2, lhood_2, smthness_2, ncc_3, lhood_3, smthness_3, tmp); break;
                    case 4: fls.critpointScore(ncc_1, lhood_1, smthness_1, ncc_2, lhood_2, smthness_2, ncc_3, lhood_3, smthness_3, ncc_4, lhood_4, smthness_4, tmp); break;
                }

                // store the critpoint scores
                endpoint_score[locationIdx] = tmp[0];
                junction_score[locationIdx] = tmp[2];


            } // else it stays zero... do nothing, they're already initialized to zero

        }

    }

}
