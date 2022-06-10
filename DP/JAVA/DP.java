import java.util.ArrayList;
import java.util.Arrays;

public class DP {
    public static void display_1D(int[] arr) 
    {
        for (int ele : arr) 
        {
            System.out.print(ele + "\t");
        }
        System.out.println();
    }
    public static void display_2D(int[][] arr) 
    {
        for (int[] ar : arr) 
        {
            display_1D(ar);
        }
        System.out.println();
    }

    public static int Fib_Simple(int n) 
    {
        if (n <= 1) 
        {
            return n;
        }

        int ans = Fib_Simple(n - 1) + Fib_Simple(n - 2);
        return ans;
    }
    public static int Fib_Memoization(int n, int[] qb) 
    {
        if (n <= 1) 
        {
            return qb[n] = n;
        }

        if (qb[n] != 0) 
        {
            return qb[n];
        }
        int ans = Fib_Memoization(n - 1, qb) + Fib_Memoization(n - 2, qb);
        return qb[n] = ans;
    }
    public static int Fib_Tabulation(int N, int[] dp) 
    {
        for (int n = 0; n <= N; n++) 
        {
            if (n <= 1) 
            {
                dp[n] = n;
                continue;
            }

            int ans = dp[n - 1] + dp[n - 2];

            dp[n] = ans;
        }
        return dp[N];
    }
    public static void Fibanocci_Series()
    {
        int n = 10;
        int res = Fib_Simple(n);
        System.out.println("Fibannocci of " + n + " is " + res);

        int[] qb = new int[n + 1];
        System.out.println("Fibannocci_Memoization of " + n + " is " + Fib_Memoization(n, qb));
        display_1D(qb);

        int[] dp = new int[n + 1];
        System.out.println("Fibannocci_Tabulation of " + n + " is " + Fib_Tabulation(n, dp));
        display_1D(dp);
    }

    public static int MazePath_Memoization(int sr, int sc, int dr, int dc, int[][] dp) 
    {
        if (sr == dr && sc == dc) 
        {
            return dp[sr][sc] = 1;
        }

        if (dp[sr][sc] != 0) 
        {
            return dp[sr][sc];
        }

        int count = 0;
        if (sc + 1 <= dc) 
        {
            count += MazePath_Memoization(sr, sc + 1, dr, dc, dp);
        }

        if (sr + 1 <= dr) 
        {
            count += MazePath_Memoization(sr + 1, sc, dr, dc, dp);
        }

        if (sr + 1 <= dr && sc + 1 <= dc) 
        {
            count += MazePath_Memoization(sr + 1, sc + 1, dr, dc, dp);
        }

        return dp[sr][sc] = count;
    }
    public static int MazePath_Tabulation(int sr, int sc, int dr, int dc, int[][] dp) 
    {
        for (sr = dr; sr >= 0; sr--) 
        {
            for (sc = dc; sc >= 0; sc--) 
            {
                if (sr == dr && sc == dc) 
                {
                    dp[sr][sc] = 1;
                    continue;
                }

                int count = 0;
                if (sc + 1 <= dc) 
                {
                    count += dp[sr][sc + 1];
                }

                if (sr + 1 <= dr) 
                {
                    count += dp[sr + 1][sc];
                }

                if (sr + 1 <= dr && sc + 1 <= dc) 
                {
                    count += dp[sr + 1][sc + 1];
                }

                dp[sr][sc] = count;
            }
        }
        return dp[0][0];
    }
    public static void MazePath_H_V_D()
    {
        int n=5, m=5;
        int[][] dp = new int[n][m];
        System.out.println("This is Memoization Path :- " + MazePath_Memoization(0, 0, n-1, m-1, dp));
        display_2D(dp);
        System.out.println("This is Tabulation Path :- " + MazePath_Tabulation(0, 0, n-1, m-1, dp));
        display_2D(dp);
    }

    public static int MazePath_MultiMove_Memoization(int sr, int sc, int dr, int dc, int[][] dp) 
    {
        if (sr == dr && sc == dc) 
        {
            return dp[sr][sc] = 1;
        }
        if (dp[sr][sc] != 0) 
        {
            return dp[sr][sc];
        }

        int count = 0;
        for (int ms = 1; ms <= dc - sc; ms++) 
        {
            count += MazePath_MultiMove_Memoization(sr, sc + ms, dr, dc, dp);
        }
        for (int ms = 1; ms <= dr - sr; ms++) 
        {
            count += MazePath_MultiMove_Memoization(sr + ms, sc, dr, dc, dp);
        }
        for (int ms = 1; ms <= dr - sr && ms <= dc - sc; ms++) 
        {
            count += MazePath_MultiMove_Memoization(sr + ms, sc + ms, dr, dc, dp);
        }

        return dp[sr][sc] = count;
    }
    public static int MazePath_MultiMove_Tabulation(int sr, int sc, int dr, int dc, int[][] dp) 
    {
        for (sr = dr; sr >= 0; sr--) 
        {
            for (sc = dc; sc >= 0; sc--) 
            {
                if (sr == dr && sc == dc) 
                {
                    dp[sr][sc] = 1;
                    continue;
                }

                int count = 0;
                for (int ms = 1; ms <= dc - sc; ms++)
                {
                    count += dp[sr][sc + ms];
                }
                for (int ms = 1; ms <= dr - sr; ms++) 
                {
                    count += dp[sr + ms][sc];
                }
                for (int ms = 1; ms <= dr - sr && ms <= dc - sc; ms++) 
                {
                    count += dp[sr + ms][sc + ms];
                }

                dp[sr][sc] = count;
            }
        }

        return dp[0][0];
    }
    public static void MazePath_MultiMove()
    {
        int n=4, m=4;
        int[][] dp = new int[n][m];
        System.out.println("This is Memoization Path of MazePath_MultiMove :- " + MazePath_MultiMove_Memoization(0, 0, n-1, m-1, dp));
        display_2D(dp);
        System.out.println("This is Tabulation Path of MazePath_MultiMove :- " + MazePath_MultiMove_Tabulation(0, 0, n-1, m-1, dp));
        display_2D(dp);
    }

    // sp--> starting point
    // ep--> Ending point
    public static int Dice_Path_Memoization(int sp, int ep, int[] dp) 
    {
        if (sp == ep) 
        {
            return dp[sp] = 1;
        }

        if (dp[sp] != 0) 
        {
            return dp[sp];
        }

        int count = 0;
        for (int jump = 1; jump <= ep - sp && jump <= 6; jump++) 
        {
            count += Dice_Path_Memoization(sp + jump, ep, dp);
        }

        return dp[sp] = count;
    }
    public static int Dice_Path_Tabulation(int sp, int ep, int[] dp) 
    {
        for (sp = ep; sp >= 0; sp--) 
        {
            if (sp == ep) 
            {
                dp[sp] = 1;
                continue;
            }

            int count = 0;
            for (int jump = 1; jump <= ep - sp && jump <= 6; jump++) 
            {
                count += dp[sp + jump];
            }

            dp[sp] = count;
        }

        return dp[0];
    }
    public static void Dice_Throw_Path()
    {
        int sp=0, ep=10;
        int[] dp = new int[ep+1];
        System.out.println("Dice path Memoization :- " + Dice_Path_Memoization(sp, ep, dp));
        display_1D(dp);
        System.out.println("Dice path Tabulation :- " + Dice_Path_Tabulation(sp, ep, dp));
        display_1D(dp);
    }

    public static int Dice_Path_Array_Memoization(int sp, int ep, int[] dp, int[] diceArray) 
    {
        if (sp == ep) 
        {
            return dp[sp] = 1;
        }

        if (dp[sp] != 0) 
        {
            return dp[sp];
        }

        int count = 0;
        for (int jump = 0; jump < diceArray.length; jump++) 
        {
            if (diceArray[jump] <= ep - sp) 
            {
                count += Dice_Path_Array_Memoization(sp + diceArray[jump], ep, dp, diceArray);
            }
        }

        return dp[sp] = count;
    }
    public static int Dice_Path_Array_Tabulation(int sp, int ep, int[] dp, int[] diceArray) 
    {
        for (sp = ep; sp >= 0; sp--) 
        {
            if (sp == ep) 
            {
                dp[sp] = 1;
                continue;
            }

            int count = 0;
            for (int dice = 0; dice < diceArray.length; dice++) 
            {
                if (sp + diceArray[dice] <= ep) 
                {
                    count += dp[sp + diceArray[dice]];
                }
            }

            dp[sp] = count;
        }

        return dp[0];
    }
    public static void Dice_with_Array()
    {
        int sp=0, ep=10;
        int[] diceArray = {1, 2, 3, 4, 5, 6};
        int[] dp = new int[ep+1];
        System.out.println("Dice Array Path Memoization :- " + Dice_Path_Array_Memoization(sp, ep, dp, diceArray));
        display_1D(dp);
        System.out.println("Dice Array Path Tabulation :- " + Dice_Path_Array_Tabulation(sp, ep, dp, diceArray));
        display_1D(dp);
    }

    // **********************************************************************************************
    // */
    // https://www.geeksforgeeks.org/friends-pairing-problem/
    public static int Friends_Pairing_Problem_Memoization(int n, int[] dp) 
    {
        if (n <= 1) 
        {
            return dp[n] = 1;
        }

        if (dp[n] != 0) 
        {
            return dp[n];
        }

        int single = Friends_Pairing_Problem_Memoization(n - 1, dp);
        int pairup = Friends_Pairing_Problem_Memoization(n - 2, dp) * (n - 1);

        int ans = single + pairup;
        return dp[n] = ans;
    }
    public static int Friends_Pairing_Problem_Tabulation(int n, int[] dp) 
    {
        int N = n;
        for (n = 0; n <= N; n++) 
        {
            if (n <= 1) 
            {
                dp[n] = 1;
                continue;
            }

            int single = dp[n - 1];
            int pairup = dp[n - 2] * (n - 1);

            int ans = single + pairup;
            dp[n] = ans;
        }
        return dp[N];
    }
    public static void Friends_Pairing_Problem()
    {
        int n = 10;
        int[] dp = new int[n+1];
        System.out.println("Memoization :- " + Friends_Pairing_Problem_Memoization(n,dp));
        display_1D(dp);
        System.out.println("Tabulation :- " + Friends_Pairing_Problem_Tabulation(n,dp));
        display_1D(dp);
    }

    // https://www.geeksforgeeks.org/gold-mine-problem/
    public static int GoldMine_Memoization(int sr, int sc, int[][] grid, int[][] dp) 
    {
        if (sc == grid[0].length - 1) 
        {
            return dp[sr][sc] = grid[sr][sc];
        }

        if (dp[sr][sc] != 0) 
        {
            return dp[sr][sc];
        }

        int[][] dir = { { -1, 1 }, { 0, 1 }, { 1, 1 } };
        int maxCoins = 0;
        for (int d = 0; d < 3; d++) 
        {
            int x = sr + dir[d][0];
            int y = sc + dir[d][1];

            if (x >= 0 && y >= 0 && x < grid.length && y < grid[0].length) 
            {
                maxCoins = Math.max(maxCoins, GoldMine_Memoization(x, y, grid, dp));
            }
        }

        return dp[sr][sc] = maxCoins + grid[sr][sc];
    }
    public static int GoldMine_Tabulation(int sr, int sc, int[][] grid, int[][] dp) 
    {
        for (sc = grid[0].length - 1; sc >= 0; sc--) 
        {
            for (sr = grid.length - 1; sr >= 0; sr--) 
            {
                if (sc == grid[0].length - 1) 
                {
                    dp[sr][sc] = grid[sr][sc];
                    continue;
                }

                int[][] dir = { { -1, 1 }, { 0, 1 }, { 1, 1 } };
                int maxCoins = 0;
                for (int d = 0; d < 3; d++) 
                {
                    int x = sr + dir[d][0];
                    int y = sc + dir[d][1];

                    if (x >= 0 && y >= 0 && x < grid.length && y < grid[0].length) 
                    {
                        maxCoins = Math.max(maxCoins, dp[x][y]);
                    }
                }

                dp[sr][sc] = maxCoins + grid[sr][sc];
            }
        }

        int max = 0;
        for (int i = 0; i < grid.length; i++) 
        {
            max = Math.max(max, dp[i][0]);
        }

        return max;
    }
    public static void GoldMine_Problem()
    { 
        int[][] grid = { 
                        {1, 3, 1, 5}, 
                        {2, 2, 4, 1}, 
                        {5, 0, 2, 3},
                        {0, 6, 1, 2} 
        }; 
        int n = grid.length; int m = grid[0].length; 
        int[][] dp = new int[n][m]; 
        int maxCoin = 0; 
        for(int i = 0; i < n; i++) 
        { 
            maxCoin = Math.max(maxCoin, GoldMine_Memoization(i, 0, grid, dp)); 
        } 
        System.out.println("Memoization :- " + maxCoin); 
        display_2D(dp);
        System.out.println("Tabulation :- " + GoldMine_Tabulation(0, 0, grid, dp));
        display_2D(dp);
    }

    // https://www.geeksforgeeks.org/count-number-of-ways-to-partition-a-set-into-k-subsets/
    public static int Count_of_Ways_Memoization(int n, int k, int[][] dp) 
    {
        if (n < k) 
        {
            return dp[k][n] = 0;
        }

        if (n == k || k == 1) 
        {
            return dp[k][n] = 1;
        }

        if (dp[k][n] != 0) 
        {
            return dp[k][n];
        }

        int newGroup = Count_of_Ways_Memoization(n - 1, k - 1, dp);
        int existingGroup = Count_of_Ways_Memoization(n - 1, k, dp) * k;
        int ans = newGroup + existingGroup;
        return dp[k][n] = ans;
    }
    public static int Count_of_Ways_Tabulation(int n, int k, int[][] dp) 
    {
        int N = n, K = k;
        for (k = 1; k <= K; k++) 
        {
            for (n = 0; n <= N; n++) 
            {
                if (n < k) 
                {
                    dp[k][n] = 0;
                    continue;
                }

                if (n == k || k == 1) 
                {
                    dp[k][n] = 1;
                    continue;
                }

                int newGroup = dp[k - 1][n - 1];
                int existingGroup = dp[k][n - 1] * k;
                int ans = newGroup + existingGroup;

                dp[k][n] = ans;
            }
        }

        return dp[K][N];
    }
    public static void Count_of_Ways(int n, int k) 
    {
        if (n < k) 
        {
            return;
        }

        int[][] dp = new int[k + 1][n + 1];
        System.out.println("Memoization :- " + Count_of_Ways_Memoization(n, k, dp));
        display_2D(dp);
        System.out.println("Tabulation :- " + Count_of_Ways_Tabulation(n, k, dp));
        display_2D(dp);
    }
    
    //****************************************************************************** */

    public static boolean[][] is_Pallindrome_SubString(String str)
    {
        int n = str.length();
        boolean[][] dp = new boolean[n][n];

        for(int gap = 0; gap < n; gap++)
        {
            for(int i = 0; i < n - gap; i++)
            {
                int j = i + gap;
                if(gap == 0)
                {
                    dp[i][j] = true;
                }
                else if(gap==1 && str.charAt(i) == str.charAt(j))
                {
                    dp[i][j] = true;
                }
                else 
                {
                    dp[i][j] = str.charAt(i) == str.charAt(j) && dp[i+1][j-1]; 
                }
            }
        }

        return dp;
    }

    public static String Leetcode_005_Longest_Pallindromic_SubString(String str)
    {
        int n = str.length();
        int[][] dp = new int[n][n];

        int maxLen = 0;
        int si=0, ei=0;
        for(int gap = 0; gap < n; gap++)
        {
            for(int i = 0; i < n - gap; i++)
            {
                int j = i + gap;
                if(gap == 0)
                {
                    dp[i][j] = 1;
                }
                else if(gap==1 && str.charAt(i) == str.charAt(j))
                {
                    dp[i][j] = 2;
                }
                else if(str.charAt(i) == str.charAt(j) && dp[i+1][j-1]!=0)
                {
                    dp[i][j] = gap + 1; 
                }

                if(dp[i][j] > maxLen)
                {
                    maxLen = dp[i][j];
                    si = i;
                    ei = j;
                }
            }
        }

        System.out.println("Maximum Length :- " + maxLen);

        for(int i = 0; i < dp.length; i++)
        {
            for(int j = 0; j < dp[0].length; j++)
            {
                System.out.print(dp[i][j] + "\t");
            }
            System.out.println();
        }

        return str.substring(si, ei+1);
    }
    public static void Leetcode_005_Longest_Pallindromic_SubString()
    {
        System.out.println(Leetcode_005_Longest_Pallindromic_SubString("abcaacbefgpgf"));
    }

    public static int Leetcode_647_Count_All_Pallindromic_SubString(String str)
    {
        int n = str.length();
        int[][] dp = new int[n][n];

        int count = 0;
        for(int gap = 0; gap < n; gap++)
        {
            for(int i = 0; i < n - gap; i++)
            {
                int j = i + gap;
                if(gap == 0)
                {
                    dp[i][j] = 1;
                }
                else if(gap==1 && str.charAt(i) == str.charAt(j))
                {
                    dp[i][j] = 2;
                }
                else if(str.charAt(i) == str.charAt(j) && dp[i+1][j-1]!=0)
                {
                    dp[i][j] = gap + 1; 
                }

                count += dp[i][j]!=0 ? 1 : 0;
            }
        }


        display_2D(dp);

        return count;  
    }
    public static void Leetcode_647_Count_All_Pallindromic_SubString()
    {
        System.out.println(Leetcode_647_Count_All_Pallindromic_SubString("aaa"));        // abccbc---> 9
    }

    public static int Leetcode_516_Longest_Pallindromic_SubSequence_Memoization(String str, int si, int ei, int[][] dp, boolean[][] isPalindrome)
    {
        if(isPalindrome[si][ei])
        {
            return dp[si][ei] = ei-si+1;
        }

        if(dp[si][ei] != 0)
        {
            return dp[si][ei];
        }

        int len=0;
        if(str.charAt(si) == str.charAt(ei))
        {
            len = Leetcode_516_Longest_Pallindromic_SubSequence_Memoization(str, si+1, ei-1, dp, isPalindrome) + 2;
        }
        else
        {
            len = Math.max(Leetcode_516_Longest_Pallindromic_SubSequence_Memoization(str, si+1, ei, dp, isPalindrome),
                           Leetcode_516_Longest_Pallindromic_SubSequence_Memoization(str, si, ei-1, dp, isPalindrome));
        }

        return dp[si][ei] = len;
    }
    public static int Leetcode_516_Longest_Pallindromic_SubSequence_Tabulation(String str, int si, int ei, int[][] dp, boolean[][] isPalindrome)
    {
        int n=str.length();
        for(int gap = 0; gap < n; gap++)
        {
            for(si = 0; si < n - gap; si++)
            {
                ei = si + gap;

                if(isPalindrome[si][ei])
                {
                    dp[si][ei] = ei-si+1;
                    continue;
                }

                int len=0;
                if(str.charAt(si) == str.charAt(ei))
                {
                    len = dp[si+1][ei-1] + 2;
                }
                else
                {
                    len = Math.max(dp[si+1][ei], dp[si][ei-1]);
                }

                dp[si][ei] = len;
            }
        }

        return dp[0][str.length()-1];
    }
    public static void Leetcode_516_Longest_Pallindromic_SubSequence()
    {
        String str = "geeksforgeeks";
        int n = str.length();
        int si = 0, ei = n - 1;
        int[][] dp = new int[n][n];

        boolean[][] isPalindrome = is_Pallindrome_SubString(str);
        System.out.println("Memoization :- " + Leetcode_516_Longest_Pallindromic_SubSequence_Memoization(str, si, ei, dp, isPalindrome));
        display_2D(dp);
        System.out.println("Tabulation :- " + Leetcode_516_Longest_Pallindromic_SubSequence_Tabulation(str, si, ei, dp, isPalindrome));
        display_2D(dp);
    }

    public static int LeetCode_115_Distinct_SubSequence_Memoization(String S, String T, int n, int m, int[][] dp)
    {
        if(m==0)
        {
            return dp[n][m] = 1;
        }
        if(m>n)
        {
            return dp[n][m] = 0;
        }

        if(dp[n][m] != 0)
        {
            return dp[n][m];
        }

        if(S.charAt(n-1) == T.charAt(m-1))
        {
            return dp[n][m] = LeetCode_115_Distinct_SubSequence_Memoization(S, T, n-1, m-1, dp)
                              + LeetCode_115_Distinct_SubSequence_Memoization(S, T, n-1, m, dp);
        }

        return dp[n][m] = LeetCode_115_Distinct_SubSequence_Memoization(S, T, n-1, m, dp);
    }
    public static int LeetCode_115_Distinct_SubSequence_Tabulation(String S, String T, int n, int m, int[][] dp)
    {
        int N=n, M=m;
        for(n=0; n<=N; n++)
        {
            for(m=0; m<=M; m++)
            {
                if(m==0)
                {
                    dp[n][m] = 1;
                    continue;
                }
                if(m>n)
                {
                    dp[n][m] = 0;
                    continue;
                }

                if(S.charAt(n-1) == T.charAt(m-1))
                {
                    dp[n][m] = dp[n-1][m-1] + dp[n-1][m];
                }
                else
                {
                    dp[n][m] = dp[n-1][m];
                }
            }
        }

        return dp[N][M];
    }
    public static int LeetCode_115_Distinct_SubSequence_02_Memoization(String S, String T, int i, int j, int[][] dp)
    {
        if(T.length()-j == 0)
        {
            return dp[i][j] = 1;
        }
        if(T.length()-j > S.length()-i)
        {
            return dp[i][j] = 0;
        }

        if(dp[i][j] != 0)
        {
            return dp[i][j];
        }

        if(S.charAt(i) == T.charAt(j))
        {
            return dp[i][j] = LeetCode_115_Distinct_SubSequence_02_Memoization(S, T, i+1, j+1, dp)
                              + LeetCode_115_Distinct_SubSequence_02_Memoization(S, T, i+1, j, dp);
        }

        return dp[i][j] = LeetCode_115_Distinct_SubSequence_02_Memoization(S, T, i+1, j, dp);
    }
    public static int LeetCode_115_Distinct_SubSequence_02_Tabulation(String S, String T, int i, int j, int[][] dp)
    {
        int N=i, M=j;
        for(i=N; i>=0; i--)
        {
            for(j=M; j>=0; j--)
            {
                if(T.length()-j == 0)
                {
                    dp[i][j] = 1;
                    continue;
                }
                if(T.length()-j > S.length()-i)
                {
                    dp[i][j] = 0;
                    continue;
                }

                if(S.charAt(i) == T.charAt(j))
                {
                    dp[i][j] = dp[i+1][j+1] + dp[i+1][j];
                }
                else
                {
                    dp[i][j] = dp[i+1][j];
                }
            }
        }

        return dp[N][M];
    }
    public static void LeetCode_115_Distinct_SubSequence()
    {
        String S = "babgbag";
        String T = "bag";

        int n = S.length();
        int m= T.length();

        int[][] dp = new int[n+1][m+1];
        
        System.out.println("Memoization :- " + LeetCode_115_Distinct_SubSequence_Memoization(S, T, n, m, dp));
        display_2D(dp);
        System.out.println("Tabulation :- " + LeetCode_115_Distinct_SubSequence_Tabulation(S, T, n, m, dp));
        display_2D(dp);
        System.out.println("Memoization_02 :- " + LeetCode_115_Distinct_SubSequence_02_Memoization(S, T, 0, 0, dp));
        display_2D(dp);
        System.out.println("Tabulation_02 :- " + LeetCode_115_Distinct_SubSequence_02_Tabulation(S, T, 0, 0, dp));
        display_2D(dp);
    }

    //Geeks: https://practice.geeksforgeeks.org/problems/count-palindromic-subsequences/1
    public static int Count_Palindromic_SubSequence_Memoization(String s, int i, int j, int[][] dp)
    {
        if(i>j)
        {
            return dp[i][j] = 0;
        }

        if(i==j)
        {
            return dp[i][j] = 1;
        }

        if(dp[i][j] != 0)
        {
            return dp[i][j];
        }

        int middleString = Count_Palindromic_SubSequence_Memoization(s, i+1, j-1, dp);
        int excludingLast = Count_Palindromic_SubSequence_Memoization(s, i, j-1, dp);
        int excludingFirst = Count_Palindromic_SubSequence_Memoization(s, i+1, j, dp);

        int ans = excludingFirst + excludingLast;

        return dp[i][j] = (s.charAt(i) == s.charAt(j)) ? ans + 1 : ans - middleString;
    }
    public static int Count_Palindromic_SubSequence_Tabulation(String s, int[][] dp)
    {
        int n = s.length();
        for(int gap = 0; gap < n; gap++)
        {
            for(int i = 0; i < n - gap; i++)
            {
                int j = i + gap;

                if(i==j)
                {
                    dp[i][j] = 1;
                    continue;
                }

                int middleString = dp[i+1][j-1];
                int excludingLast = dp[i][j-1];
                int excludingFirst = dp[i+1][j];

                int ans = excludingFirst + excludingLast;

                dp[i][j] = (s.charAt(i) == s.charAt(j)) ? ans + 1 : ans - middleString;
            }
        }

        return dp[0][s.length()-1];
    }
    public static void Count_Palindromic_SubSequence()
    {
        String s = "abckycbc";
        int n=s.length();
        int[][] dp = new int[n][n];
        int i = 0, j = n-1;
        System.out.println("Memoization :- " + Count_Palindromic_SubSequence_Memoization(s, i, j, dp));
        display_2D(dp);
        System.out.println("Tabulation :- " + Count_Palindromic_SubSequence_Tabulation(s, dp));
        display_2D(dp);
    }

    // *****************************************************************************
    // ******************************************************************************/

    public static int LeetCode_70_Climbing_Stairs_Memoization(int n, int[] dp) 
    {
        if (n <= 1) 
        {
            return dp[n] = 1;
        }

        if (dp[n] != 0) 
        {
            return dp[n];
        }

        int ans = LeetCode_70_Climbing_Stairs_Memoization(n - 1, dp)
                + LeetCode_70_Climbing_Stairs_Memoization(n - 2, dp);

        return dp[n] = ans;
    }
    public static int LeetCode_70_Climbing_Stairs_Tabulation(int n, int[] dp) 
    {
        int N = n;
        for (n = 0; n <= N; n++) 
        {
            if (n <= 1) 
            {
                dp[n] = 1;
                continue;
            }

            int ans = dp[n - 1] + dp[n - 2];

            dp[n] = ans;
        }

        return dp[N];
    }
    public static int LeetCode_70_Climbing_Stairs_Better(int n) 
    {
        int a = 1;
        int b = 1;
        int sum = 0;

        for (int i = 2; i <= n; i++) 
        {
            sum = a + b;
            a = b;
            b = sum;
        }

        return sum;
    }
    public static void LeetCode_70_Climbing_Stairs()
    {
        int n=7;
        int[] dp = new int[n+1];
        System.out.println("Memoization :- " + LeetCode_70_Climbing_Stairs_Memoization(n, dp));
        display_1D(dp);
        System.out.println("Tabulation :- " + LeetCode_70_Climbing_Stairs_Tabulation(n, dp));
        display_1D(dp);
        System.out.println("Better :- " + LeetCode_70_Climbing_Stairs_Better(n));
    }

    public static int LeetCode_746_Min_Cost_Climbing_Stairs_Memoization(int n, int[] dp, ArrayList<Integer> cost) 
    {
        if (n <= 1) 
        {
            return dp[n] = cost.get(n);
        }

        if (dp[n] != 0) 
        {
            return dp[n];
        }

        int ans = Math.min(LeetCode_746_Min_Cost_Climbing_Stairs_Memoization(n - 1, dp, cost),
                LeetCode_746_Min_Cost_Climbing_Stairs_Memoization(n - 2, dp, cost));

        return dp[n] = ans + cost.get(n);
    }
    public static int LeetCode_746_Min_Cost_Climbing_Stairs_Tabulation(int n, int[] dp, ArrayList<Integer> cost) 
    {
        int N = n;
        for (n = 0; n <= N; n++) 
        {
            if (n <= 1) 
            {
                dp[n] = cost.get(n);
                continue;
            }

            int ans = Math.min(dp[n - 1], dp[n - 2]);

            dp[n] = ans + cost.get(n);
        }
        return dp[N];
    }
    public static void LeetCode_746_Min_Cost_Climbing_Stairs()
    {
        ArrayList<Integer> cost = new ArrayList<>();
        cost.addAll(Arrays.asList(1, 100, 1, 1, 1, 100, 1, 1, 100, 1));
        int n = cost.size();
        cost.add(0);
        int[] dp = new int[n+1];
        System.out.println("Memoization :- " + LeetCode_746_Min_Cost_Climbing_Stairs_Memoization(n, dp, cost));
        display_1D(dp);
        System.out.println("Tabulation :- " + LeetCode_746_Min_Cost_Climbing_Stairs_Tabulation(n, dp, cost));
        display_1D(dp);
    }

    public static int LeetCode_64_min_Path_Sum_Memoization(int sr, int sc, int[][] grid, int[][] dp) 
    {
        if (sr == grid.length - 1 && sc == grid[0].length - 1) 
        {
            return dp[sr][sc] = grid[sr][sc];
        }

        if (dp[sr][sc] != 0) 
        {
            return dp[sr][sc];
        }

        int minCost = Integer.MAX_VALUE;
        if (sr + 1 <= grid.length - 1) 
        {
            minCost = Math.min(minCost, LeetCode_64_min_Path_Sum_Memoization(sr + 1, sc, grid, dp));
        }

        if (sc + 1 <= grid[0].length - 1) 
        {
            minCost = Math.min(minCost, LeetCode_64_min_Path_Sum_Memoization(sr, sc + 1, grid, dp));
        }

        return dp[sr][sc] = minCost + grid[sr][sc];
    }
    public static int LeetCode_64_min_Path_Sum_Tabulation(int sr, int sc, int[][] grid, int[][] dp) 
    {
        for (sr = grid.length - 1; sr >= 0; sr--) 
        {
            for (sc = grid[0].length - 1; sc >= 0; sc--) 
            {
                if (sr == grid.length - 1 && sc == grid[0].length - 1) 
                {
                    dp[sr][sc] = grid[sr][sc];
                    continue;
                }

                int minCost = Integer.MAX_VALUE;
                if (sr + 1 <= grid.length - 1) 
                {
                    minCost = Math.min(minCost, dp[sr + 1][sc]);
                }

                if (sc + 1 <= grid[0].length - 1) 
                {
                    minCost = Math.min(minCost, dp[sr][sc + 1]);
                }

                dp[sr][sc] = minCost + grid[sr][sc];
            }
        }
        return dp[0][0];
    }
    public static void LeetCode_64_min_Path_Sum()
    {
        int[][] grid = {
                       {1,3,1},
                       {1,5,1},
                       {4,2,1}
        };
        int n = grid.length;
        int[][] dp = new int[n][grid[0].length];
        System.out.println("Memoization :- " + LeetCode_64_min_Path_Sum_Memoization(0, 0, grid, dp));
        display_2D(dp);
        System.out.println("Tabulation :- " + LeetCode_64_min_Path_Sum_Tabulation(0, 0, grid, dp));
        display_2D(dp);
    }

    // ****************************************************************************/
    // ------------------------------LeetCode Problems-----------------------------
    // ****************************************************************************/
    public static void LeetCode_Problems() 
    {
        //LeetCode_70_Climbing_Stairs();

        //LeetCode_746_Min_Cost_Climbing_Stairs();

        //LeetCode_64_min_Path_Sum();

    }
    // ****************************************************************************/
    // -----------------------------------END-------------------------------------
    // ****************************************************************************/

    public static void DP_SubString_SubSequence()
    {
        //Leetcode_005_Longest_Pallindromic_SubString();
        
        //Leetcode_647_Count_All_Pallindromic_SubString();
        
        //Leetcode_516_Longest_Pallindromic_SubSequence();

        LeetCode_115_Distinct_SubSequence();

        //Count_Palindromic_SubSequence();
    }

    public static void Geeks_For_Geeks() 
    {
        //Friends_Pairing_Problem();

        //GoldMine_Problem();

        //Count_of_Ways(7, 3);      // int n=7, k=3;

    }

    public static void DP_Basics_1() 
    {
        //MazePath_H_V_D();

        //MazePath_MultiMove();

        //Dice_Throw_Path();

        //Dice_with_Array();

    }

    public static void DP_Basics_0() 
    {
        //Fibanocci_Series(); 
    }

    public static void solve() {
        //DP_Basics_0();

        //DP_Basics_1();

        //LeetCode_Problems();

        //Geeks_For_Geeks();

        DP_SubString_SubSequence();
    }

    public static void main(String[] args) 
    {
        solve();
    }
}