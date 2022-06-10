public class Recursion
{
    public static void Print_Increasing(int a, int b)
    {
        if(a == b+1)
        {
            return;
        }
        System.out.print(a + " ");
        Print_Increasing(a+1, b);
    }
    public static void Print_Decreasing(int a, int b)
    {
        if(a == b+1)
        {
            return;
        }
        
        Print_Decreasing(a+1, b);
        System.out.print(a + " ");
    }
    public static void Recursion_Basics()
    {
        //Print_Increasing(0, 10);
        //Print_Decreasing(0, 10);
    }
    public static void solve()
    {
        Recursion_Basics();
    }
    public static void main(String[] args)
    {
        solve();
    }
}