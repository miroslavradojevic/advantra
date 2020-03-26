package advantra.filter;

public class FastMedian {

    /**
     * background median estimation ( Torben's algorithm )
     *
     * The following code is public domain.
     * Algorithm by Torben Mogensen, implementation by N. Devillard.
     * This code in public domain.
     */
    public static double median_Torben(float[] m)
    {
        int n = m.length;
        int    less, greater, equal;
        double min, max, guess, maxltguess, mingtguess;
        min = max = m[0] ;
        for (int i=1 ; i<n ; i++)
        {
            if (m[i]<min) min=m[i];
            if (m[i]>max) max=m[i];
        }
        while (true)
        {
            guess = (min+max)/2;
            less = 0; greater = 0; equal = 0;
            maxltguess = min ;
            mingtguess = max ;
            for (int i=0; i<n; i++)
            {
                if (m[i]<guess)
                {
                    less++;
                    if (m[i]>maxltguess) maxltguess = m[i] ;
                } else if (m[i]>guess) {
                    greater++;
                    if (m[i]<mingtguess) mingtguess = m[i] ;
                } else equal++;
            }
            if (less <= (n+1)/2 && greater <= (n+1)/2) break ;
            else if (less>greater) max = maxltguess ;
            else min = mingtguess;
        }
        if (less >= (n+1)/2) return maxltguess;
        else if (less+equal >= (n+1)/2) return guess;
        else return mingtguess;
    }

    /**
     * background median estimation ( Wirth's algorithm )
     * Title: Algorithms + data structures = programs
     * Publisher: Englewood Cliffs: Prentice-Hall, 1976
     * Physical description: 366 p.
     * Series: Prentice-Hall Series in Automatic Computation
     */
    public static double median_Wirth(float[] a)
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

	
}
