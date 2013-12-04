package detection3d;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 12/4/13
 * Time: 2:34 PM
 */
public class Stat {

    /*
        collection of statistical functions used in detection
     */

    public static float median(float[] a)
    {
        int n = a.length;
        int i, j, l, m, k;
        double x;
        if (n % 2 == 0) k = (n/2)-1;
        else k = (n/2);
        l=0 ; m=n-1 ;
        while (l < m)
        {
            x=a[k] ;
            i = l ;
            j = m ;
            do
            {
                while (a[i] < x) i++ ;
                while (x < a[j]) j-- ;
                if (i <= j) {
                    float temp = a[i];
                    a[i] = a[j];
                    a[j] = temp;
                    i++ ; j-- ;
                }
            } while (i <= j) ;
            if (j < k) l = i ;
            if (k < i) m = j ;
        }
        return a[k] ;
    }

    public static float quantile(float[] a, int ratioNum, int ratioDen) // ratioNum/ratioDen first ones
    {
        int n = a.length;
        int i, j, l, m, k;
        double x;

        if ((ratioNum*n) % ratioDen == 0) k = ((ratioNum*n)/ratioDen)-1;
        else k = (ratioNum*n)/ratioDen;

        l=0 ; m=n-1 ;
        while (l < m)
        {
            x=a[k] ;
            i = l ;
            j = m ;
            do
            {
                while (a[i] < x) i++ ;
                while (x < a[j]) j-- ;
                if (i <= j) {
                    float temp = a[i];
                    a[i] = a[j];
                    a[j] = temp;
                    i++ ; j-- ;
                }
            } while (i <= j) ;
            if (j < k) l = i ;
            if (k < i) m = j ;
        }
        return a[k] ;
    }


    public static float average(float[] a)
    {
        float meanVal = 0;
        for (int i=0; i<a.length; i++)
            meanVal += a[i];
        return meanVal/a.length;

    }

    public static float std(float[] in, float avg)
    {
        float std = 0;
        for (int i=0; i<in.length; i++) {
            std += (in[i]-avg)*(in[i]-avg);
        }
        std /= in.length;
        std = (float) Math.sqrt(std);
        return std;
    }


}
