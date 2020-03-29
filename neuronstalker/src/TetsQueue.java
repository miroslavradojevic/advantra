import java.util.*;

/**
 * Created by miroslav on 3-4-15.
 */
public class TetsQueue {

    public static void main(String[] args) {

        System.out.println("testing...");

        BfsQueue bfs_queue = new BfsQueue();

        System.out.println("enqueue:");
        int[] e = new int[]{2,3};
        bfs_queue.enqueue(e);           System.out.println(Arrays.toString(e));

        e = new int[]{3,4};
        bfs_queue.enqueue(e);           System.out.println(Arrays.toString(e));

        e = new int[]{99,99};
        bfs_queue.enqueue(e);           System.out.println(Arrays.toString(e));

        e = new int[]{88, 88};
        bfs_queue.enqueue(e);           System.out.println(Arrays.toString(e));

        System.out.println(bfs_queue.size());
        System.out.println("----------------------------");

        while (bfs_queue.hasItems()) System.out.println(Arrays.toString((int[]) bfs_queue.dequeue()));

        System.out.println("****************************");

        ArrayList<Boolean> bb= new ArrayList<Boolean>();
        bb.add(true);
        bb.add(true);
        bb.add(false);

        System.out.println(Collections.frequency(bb, false) + " elements");

    }

}
