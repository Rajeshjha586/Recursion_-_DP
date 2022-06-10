#include <iostream>
#include <vector>
#include <list>
#include <algorithm>

using namespace std;

void display_1D(vector<int> &arr)
{
    for (int ele : arr)
    {
        cout << ele << "\t";
    }
    cout << endl;
}

void display_2D(vector<vector<int>> &arr)
{
    for (vector<int> ar : arr)
    {
        display_1D(ar);
    }
    cout << endl;
}

int Fib_Memoization(int n, vector<int> &dp)
{
    if (n <= 1)
    {
        return dp[n] = n;
    }

    if (dp[n] != 0)
    {
        return dp[n];
    }

    int ans = Fib_Memoization(n - 1, dp) + Fib_Memoization(n - 2, dp);
    return dp[n] = ans;
}
int Fib_Tabulation(int n, vector<int> &dp)
{
    int N = n;
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
int Fib_Bttr(int n)
{
    int a = 0;
    int b = 1;
    int sum = 0;
    for (int i = 0; i < n; i++)
    {
        sum = a + b;
        a = b;
        b = sum;
    }

    return a;
}
void Fibanocci_Series()
{
    int n = 5;
    vector<int> dp(n + 1, 0);
    cout << "Memoization of Fibanocci :- " << Fib_Memoization(n, dp) << endl;
    display_1D(dp);
    cout << "Tabulation of Fibanocci :- " << Fib_Tabulation(n, dp) << endl;
    display_1D(dp);
    cout << "Fibanocci better :- " << Fib_Bttr(n) << " " << endl;
}

int MazePath_Memoization(int sr, int sc, int dr, int dc, vector<vector<int>> &dp)
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
int MazePath_Tabulation(int sr, int sc, int dr, int dc, vector<vector<int>> &dp)
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
void MazePath_H_V_D()
{
    int n=3, m=3;
    vector<vector<int>> dp(n, vector<int>(m, 0));
    cout << "Memoization of MazePath --> Path :- " << MazePath_Memoization(0, 0, n-1, m-1, dp) << endl;
    display_2D(dp);
    cout << "Tabulation of MazePath --> Path :- " << MazePath_Tabulation(0, 0, n-1, m-1, dp) << endl;
    display_2D(dp);

}

int MazePath_MultiMove_Memoization(int sr, int sc, int dr, int dc, vector<vector<int>> &dp)
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
int MazePath_MultiMove_Tabulation(int sr, int sc, int dr, int dc, vector<vector<int>> &dp)
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
void MazePath_MultiMove()
{
    int n=3, m=3;
    vector<vector<int>> dp(n, vector<int>(m, 0));
    cout << "Memoization of MazePath_MultiMove --> Path :- " << MazePath_MultiMove_Memoization(0, 0, n-1, m-1, dp) << endl;
    display_2D(dp);
    cout << "Tabulation of MazePath_MultiMove --> Path :- " << MazePath_MultiMove_Tabulation(0, 0, n-1, m-1, dp) << endl;
    display_2D(dp);
}

int Dice_Memoization(int sp, int ep, vector<int> &dp)
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
    for (int dice = 1; dice <= ep - sp && dice <= 6; dice++)
    {
        count += Dice_Memoization(sp + dice, ep, dp);
    }

    return dp[sp] = count;
}
int Dice_Tabulation(int sp, int ep, vector<int> &dp)
{
    for (sp = ep; sp >= 0; sp--)
    {
        if (sp == ep)
        {
            dp[sp] = 1;
            continue;
        }

        int count = 0;
        for (int dice = 1; dice <= ep - sp && dice <= 6; dice++)
        {
            count += dp[sp + dice];
        }

        dp[sp] = count;
    }

    return dp[0];
}
int Dice_best(int sp, int ep)
{
    list<int> ll;
    for (sp = ep; sp >= 0; sp--)
    {
        if (sp > ep - 2)
        {
            ll.push_front(1);
            continue;
        }

        if (ll.size() <= 6)
            ll.push_front(2 * ll.front());
        else
        {
            ll.push_front(2 * ll.front() - ll.back());
            ll.pop_back();
        }
    }

    return ll.front();
}
void Dice_Throw()
{
    int sp=0, ep=10;
    vector<int> dp(ep+1, 0);
    cout << "Dice Memoization --> Path :- " << Dice_Memoization(sp, ep, dp) << endl;
    display_1D(dp);
    cout << "Dice Tabulation --> Path :- " << Dice_Tabulation(sp, ep, dp) << endl;
    display_1D(dp);
    cout << "Dice Best --> Path :- " << Dice_best(sp, ep) << endl;
}

int Dice_Array_Memoization(int sp, int ep, vector<int> &dp, vector<int> &diceArray)
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
    for (int dice = 0; diceArray[dice] <= ep - sp && dice < diceArray.size(); dice++)
    {
        count += Dice_Array_Memoization(sp + diceArray[dice], ep, dp, diceArray);
    }

    return dp[sp] = count;
}
int Dice_Array_Tabulation(int sp, int ep, vector<int> &dp, vector<int> &diceArray)
{
    for (sp = ep; sp >= 0; sp--)
    {
        if (sp == ep)
        {
            dp[sp] = 1;
            continue;
        }

        int count = 0;
        for (int dice = 0; diceArray[dice] <= ep - sp && dice < diceArray.size(); dice++)
        {
            count += dp[sp + diceArray[dice]];
        }

        dp[sp] = count;
    }
}
void Dice_With_Array()
{
    int sp=0, ep=10;
    vector<int> diceArray{1, 2, 3, 4, 5, 6};
    vector<int> dp(ep+1, 0);
    cout << "Dice_Array Memoization --> Path :- " << Dice_Array_Memoization(sp, ep, dp, diceArray) << endl;
    display_1D(dp);
    cout << "Dice_Array Tabulation --> Path :- " << Dice_Array_Tabulation(sp, ep, dp, diceArray) << endl;
    display_1D(dp);
}


//******************************************************************************************

// https://www.geeksforgeeks.org/friends-pairing-problem/
int Friends_Pairing_Problem_Memoization(int n, vector<int>& dp)
{
    if(n<=1)
    {
        return dp[n] = 1;
    }

    if(dp[n] != 0)
    {
        return dp[n];
    }

    int Single = Friends_Pairing_Problem_Memoization(n-1, dp);
    int PairUp = Friends_Pairing_Problem_Memoization(n-2, dp) * (n-1);
    int ans = Single + PairUp;
    return dp[n] = ans;
}
int Friends_Pairing_Problem_Tabulation(int n, vector<int>& dp)
{
    int N=n;
    for(n=0; n<=N; n++)
    {
        if(n<=1)
        {
            dp[n] = 1;
            continue;
        }

        int Single = dp[n-1];
        int PairUp = dp[n-2] * (n-1);
        int ans = Single + PairUp;

        dp[n] = ans;
    }

    return dp[N];
}
void Friends_Pairing_Problem()
{
    int n=4;
    vector<int> dp(n+1, 0);
    cout << "Memoization :- " << Friends_Pairing_Problem_Memoization(n, dp) << endl;
    display_1D(dp);
    cout << "Tabulation :- " << Friends_Pairing_Problem_Tabulation(n, dp) << endl;
    display_1D(dp);
}

// https://www.geeksforgeeks.org/gold-mine-problem/
int GoldMine_Memoization(int sr, int sc, vector<vector<int>>& grid, vector<vector<int>>& dp)
{
    if(sc == grid[0].size()-1)
    {
        return dp[sr][sc] = grid[sr][sc];
    }

    if(dp[sr][sc] != 0)
    {
        return dp[sr][sc];
    }

    int dir[3][2] = {{-1, 1}, {0, 1}, {1, 1}};
    int maxCoins = 0;
    for(int d=0; d<3; d++)
    {
        int x = sr + dir[d][0];
        int y = sc + dir[d][1];

        if(x>=0 && y>=0 && x<grid.size() && y<grid[0].size())
        {
            maxCoins = max(maxCoins, GoldMine_Memoization(x, y, grid, dp));
        }
    }

    return dp[sr][sc] = maxCoins + grid[sr][sc];
}
int GoldMine_Tabulation(int sr, int sc, vector<vector<int>>& grid, vector<vector<int>>& dp)
{
    for(sr=grid.size()-1; sr>=0; sr--)
    {
        for(sc=grid[0].size()-1; sc>=0; sc--)
        {
            if(sc == grid[0].size()-1)
            {
                dp[sr][sc] = grid[sr][sc];
                continue;
            }

            int dir[3][2] = {{-1, 1}, {0, 1}, {1, 1}};
            int maxCoins = 0;
            for(int d=0; d<3; d++)
            {
                int x = sr + dir[d][0];
                int y = sc + dir[d][1];

                if(x>=0 && y>=0 && x<grid.size() && y<grid[0].size())
                {
                    maxCoins = max(maxCoins, dp[x][y]);
                }
            }

            dp[sr][sc] = maxCoins + grid[sr][sc];
        }
    }

    int maxCoin = 0;
    for(int i=0; i<grid.size(); i++)
    {
        maxCoin = max(maxCoin, dp[i][0]);
    }
    return maxCoin;
}
void GoldMine_Problem()
{
    vector<vector<int>> grid = {
                               {1, 3, 1, 5}, 
                               {2, 2, 4, 1}, 
                               {5, 0, 2, 3},
                               {0, 6, 1, 2} 
    };

    int n=grid.size(), m=grid[0].size();
    vector<vector<int>> dp(n, vector<int>(m, 0));

    int maxCoin=0;
    for(int i=0; i<n; i++)
    {
        maxCoin = max(maxCoin, GoldMine_Memoization(i, 0, grid, dp));
    }
    cout << "Memoization :- " << maxCoin << endl; 
    display_2D(dp);
    cout << "Tabulation :- " << GoldMine_Tabulation(0, 0, grid, dp) << endl; 
    display_2D(dp);
}

// https://www.geeksforgeeks.org/count-number-of-ways-to-partition-a-set-into-k-subsets/
int Count_Number_of_Ways_Memoization(int n, int k, vector<vector<int>>& dp)
{
    if(n<k)
    {
        return dp[k][n] = 0;
    }

    if(n==k || k==1)
    {
        return dp[k][n] = 1; 
    }

    if(dp[k][n] != 0)
    {
        return dp[k][n];
    }

    int newGrp = Count_Number_of_Ways_Memoization(n-1, k-1, dp);
    int ExistingGrp = Count_Number_of_Ways_Memoization(n-1, k, dp) * k;
    int ans = newGrp + ExistingGrp;

    return dp[k][n] = ans;
}
int Count_Number_of_Ways_Tabulation(int n, int k, vector<vector<int>>& dp)
{
    int N = n, K = k;
    for (k = 1; k <= K; k++) 
    {
        for (n = 0; n <= N; n++) 
        {
            if(n<k)
            {
                dp[k][n] = 0;
                continue;
            }

            if(n==k || k==1)
            {
                dp[k][n] = 1;
                continue; 
            }

            int newGrp = dp[k-1][n-1];
            int ExistingGrp = dp[k][n-1] * k;
            int ans = newGrp + ExistingGrp;

            dp[k][n] = ans;
        }
    }

    return dp[K][N];
}
void Count_Number_of_Ways(int n, int k)
{
    if(n<k)
    {
        return;
    }

    vector<vector<int>> dp(k+1, vector<int>(n+1, 0));

    cout << "Memoization :- " << Count_Number_of_Ways_Memoization(n, k, dp) << endl;
    display_2D(dp);
    cout << "Tabulation :- " << Count_Number_of_Ways_Tabulation(n, k, dp) << endl;
    display_2D(dp);
}

//******************************************************************************************
vector<vector<bool>> is_Palindrome_SubString(string str)
{
    int n = str.length();
    vector<vector<bool>> dp(n, vector<bool>(n, 0));
    for(int gap=0; gap<n; gap++)
    {
        for(int i=0; i<n-gap; i++)
        {
            int j = i + gap;
            
            if(gap == 0)
            {
                dp[i][j] = true;
            }
            else if(gap == 1 && str[i] == str[j])
            {
                dp[i][j] = true;
            }
            else
            {
                dp[i][j] = str[i]==str[j] && dp[i+1][j-1];
            }
            
        }
    }

    return dp;
}
string Leetcode_005_Longest_Pallindromic_SubString_Tabulation(string str)
{
    int n = str.length();
    vector<vector<int>> dp(n, vector<int>(n, 0));
    
    int maxLen = 0;
    int si=0, ei=0; 

    for(int gap=0; gap<n; gap++)
    {
        for(int i=0; i<n-gap; i++)
        {
            int j = i + gap; 
            if(gap == 0)
            {
                dp[i][j] = 1;
            }
            else if(gap == 1 && str[i] == str[j])
            {
                dp[i][j] = 2;
            }
            else if(str[i]==str[j] && dp[i+1][j-1]!=0)
            {
                dp[i][j] = gap+1;
            }

            if(dp[i][j] > maxLen)
            {
                maxLen = dp[i][j];
                si = i;
                ei = j;
            }
            
        }
    }

    cout << "Maximum Length :- " << maxLen << endl;

    display_2D(dp);

    return str.substr(si, (ei-si+1));   
}
void Leetcode_005_Longest_Pallindromic_SubString()
{
    cout << "Tabulation :- " << Leetcode_005_Longest_Pallindromic_SubString_Tabulation("abcaacbefgpgf") << endl;
}

int Leetcode_647_Count_All_Pallindromic_SubString_Tabulation(string str)
{
    int n = str.length();
    vector<vector<int>> dp(n, vector<int>(n, 0));
    
    int count=0;

    for(int gap=0; gap<n; gap++)
    {
        for(int i=0; i<n-gap; i++)
        {
            int j = i + gap; 
            if(gap == 0)
            {
                dp[i][j] = 1;
            }
            else if(gap == 1 && str[i] == str[j])
            {
                dp[i][j] = 2;
            }
            else if(str[i]==str[j] && dp[i+1][j-1]!=0)
            {
                dp[i][j] = gap+1;
            }

            count += dp[i][j] !=0 ? 1 : 0;
            
        }
    }

    display_2D(dp);

    return count;   
}
void Leetcode_647_Count_All_Pallindromic_SubString()
{
    cout << "Tabulation :- " << Leetcode_647_Count_All_Pallindromic_SubString_Tabulation("aaa") << endl;
}

int Leetcode_516_Longest_Pallindromic_SubSequence_Memoization(string str, int si, int ei, vector<vector<int>>& dp, vector<vector<bool>>& is_Palindrome)
{
    if(is_Palindrome[si][ei])
    {
        return dp[si][ei] = ei-si+1;
    }

    if(si==ei)
    {
        return dp[si][ei] = 1;
    }

    if(dp[si][ei] != 0)
    {
        return dp[si][ei];
    }

    int len = 0;
    if(str[si] == str[ei])
    {
        len = Leetcode_516_Longest_Pallindromic_SubSequence_Memoization(str, si+1, ei-1, dp, is_Palindrome) + 2;
    }
    else
    {
        len = max(Leetcode_516_Longest_Pallindromic_SubSequence_Memoization(str, si+1, ei, dp, is_Palindrome), 
              Leetcode_516_Longest_Pallindromic_SubSequence_Memoization(str, si, ei-1, dp, is_Palindrome));
    }
    
    return dp[si][ei] = len;
}
int Leetcode_516_Longest_Pallindromic_SubSequence_Tabulation(string str, int si, int ei, vector<vector<int>>& dp, vector<vector<bool>>& is_Palindrome)
{
    int n = str.length();
    for(int gap=0; gap<n; gap++)
    {
        for(si=0; si<n-gap; si++)
        {
            ei = si + gap;
            if(is_Palindrome[si][ei])
            {
                dp[si][ei] = ei-si+1;
                continue;
            }

            if(si==ei)
            {
                dp[si][ei] = 1;
                continue;
            }

            int len = 0;
            if(str[si] == str[ei])
            {
                len = dp[si+1][ei-1] + 2;
            }
            else
            {
                len = max(dp[si+1][ei], dp[si][ei-1]);
            }
            
            dp[si][ei] = len;
        }
    }
    
    return dp[0][n-1];
}
void Leetcode_516_Longest_Pallindromic_SubSequence()
{
    string str = "geeksforgeeks";
    int n = str.length();
    int si=0, ei=n-1;

    vector<vector<int>> dp(n, vector<int>(n, 0));

    vector<vector<bool>> is_Palindrome = is_Palindrome_SubString(str);
    cout << "Memoization :- " << Leetcode_516_Longest_Pallindromic_SubSequence_Memoization(str, si, ei, dp, is_Palindrome) << endl;
    display_2D(dp);
    cout << "Tabulation :- " << Leetcode_516_Longest_Pallindromic_SubSequence_Tabulation(str, si, ei, dp, is_Palindrome) << endl;
    display_2D(dp);
}

int LeetCode_115_Distinct_SubSequence_Memoization(string S, string T, int n, int m, vector<vector<int>> &dp)
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

    if(S[n-1] == T[m-1])
    {
        return dp[n][m] = LeetCode_115_Distinct_SubSequence_Memoization(S, T, n-1, m-1, dp)
                              + LeetCode_115_Distinct_SubSequence_Memoization(S, T, n-1, m, dp);
    }

    return dp[n][m] = LeetCode_115_Distinct_SubSequence_Memoization(S, T, n-1, m, dp);
}
int LeetCode_115_Distinct_SubSequence_Tabulation(string S, string T, int n, int m, vector<vector<int>> &dp)
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

            if(S[n-1] == T[m-1])
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
int LeetCode_115_Distinct_SubSequence_Memoization_02(string S, string T, int i, int j, vector<vector<int>> &dp)
{
    if (T.length() - j == 0)
    {
        return dp[i][j] = 1;
    }
    if(T.length() - j > S.length() - i)
    {
        return dp[i][j] = 0;
    }

    if (dp[i][j] != 0)
    {
        return dp[i][j];
    }

    if (S[i] == T[j])
    {
         return dp[i][j] = LeetCode_115_Distinct_SubSequence_Memoization_02(S, T, i + 1, j + 1, dp) 
                          + LeetCode_115_Distinct_SubSequence_Memoization_02(S, T, i + 1, j, dp);
    }

    return dp[i][j] = LeetCode_115_Distinct_SubSequence_Memoization_02(S, T, i + 1, j, dp);
}
int LeetCode_115_Distinct_SubSequence_Tabulation_02(string S, string T, int i, int j, vector<vector<int>> &dp)
{
    int N=i, M=j;
    for(i=N; i>=0; i--)
    {
        for(j=M; j>=0; j--)
        {
            if (T.length() - j == 0)
            {
                dp[i][j] = 1;
                continue;
            }
            if(T.length() - j > S.length() - i)
            {
                dp[i][j] = 0;
                continue;
            }

            if (S[i] == T[j])
            {
                dp[i][j] = dp[i+1][j+1] + dp[i+1][j];
            }
            else
            {
                dp[i][j] = dp[i+1][j];
            }
        }
    }
    return dp[0][0];
}
void LeetCode_115_Distinct_SubSequence()
{
    string S = "babgbag";
    string T = "bag";

    int n = S.length();
    int m= T.length();

    vector<vector<int>> dp(n+1, vector<int>(m+1, 0));    //vector<vector<int>> dp(n+1, vector<int>(m+1, -1)); check call unwanted

    //cout << "Memoization :- " << LeetCode_115_Distinct_SubSequence_Memoization(S, T, n, m, dp) << endl;
    //display_2D(dp);
    //cout << "Tabulation :- " << LeetCode_115_Distinct_SubSequence_Tabulation(S, T, n, m, dp) << endl;
    //display_2D(dp);
    cout << "Memoization_02 :- " << LeetCode_115_Distinct_SubSequence_Memoization_02(S, T, 0, 0, dp) << endl;
    display_2D(dp);
    cout << "Tabulation_02 :- " << LeetCode_115_Distinct_SubSequence_Tabulation_02(S, T, 0, 0, dp) << endl;
    display_2D(dp);
}

//**************************************************************
//---------------------------LeetCode---------------------------

// 1 or 2 jumps are allowed
int LeetCode_70_ClimbStairs_Memoization(int n, vector<int> &dp)
{
    if (n <= 1)
    {
        return dp[n] = 1;
    }

    if (dp[n] != 0)
    {
        return dp[n];
    }

    int ans = LeetCode_70_ClimbStairs_Memoization(n - 1, dp) + LeetCode_70_ClimbStairs_Memoization(n - 2, dp);

    return dp[n] = ans;
}
int LeetCode_70_ClimbStairs_Tabulation(int n, vector<int> &dp)
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
int LeetCode_70_ClimbStairs_Bttr(int n)
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
void LeetCode_70_ClimbStairs()
{
    int n=7;
    vector<int> dp(n+1, 0);
    cout <<"Memoization :- Path " << LeetCode_70_ClimbStairs_Memoization(n, dp) << endl;
    display_1D(dp);
    cout <<"Tabulation :- Path " << LeetCode_70_ClimbStairs_Tabulation(n, dp) << endl;
    display_1D(dp);
    cout <<"Better Sol :- Path " << LeetCode_70_ClimbStairs_Bttr(n) << endl;
}


int LeetCode_746_Min_Cost_Climbing_Stairs_Memoization(int n, vector<int> &dp, vector<int> &cost)
{
    if (n <= 1)
    {
        return dp[n] = cost[n];
    }

    if (dp[n] != 0)
    {
        return dp[n];
    }

    int ans = min(LeetCode_746_Min_Cost_Climbing_Stairs_Memoization(n - 1, dp, cost),
                  LeetCode_746_Min_Cost_Climbing_Stairs_Memoization(n - 2, dp, cost));

    return dp[n] = ans + cost[n];
}
int LeetCode_746_Min_Cost_Climbing_Stairs_Tabulation(int n, vector<int> &dp, vector<int> &cost)
{
    int N = n;
    for (n = 0; n <= N; n++)
    {
        if (n <= 1)
        {
            dp[n] = cost[n];
            continue;
        }

        int ans = min(dp[n - 1], dp[n - 2]);

        dp[n] = ans + cost[n];
    }

    return dp[N];
}
void LeetCode_746_Min_Cost_Climbing_Stairs()
{
    vector<int> cost{1, 100, 1, 1, 1, 100, 1, 1, 100, 1};
    int n = cost.size();
    cost.push_back(0);
    vector<int> dp(n + 1, 0);
    cout << "Memoization :- " << LeetCode_746_Min_Cost_Climbing_Stairs_Memoization(n, dp, cost) << endl;
    display_1D(dp);
    cout << "Tabulation :- " << LeetCode_746_Min_Cost_Climbing_Stairs_Tabulation(n, dp, cost) << endl;
    display_1D(dp);
}


int LeetCode_64_Min_Path_Sum_Memoization(int sr, int sc, vector<vector<int>> &grid, vector<vector<int>> &dp)
{
    if (sr == grid.size() - 1 && sc == grid[0].size() - 1)
    {
        return dp[sr][sc] = grid[sr][sc];
    }

    if (dp[sr][sc] != 0)
    {
        return dp[sr][sc];
    }

    int minCost = 1e8;
    if (sr + 1 < grid.size())
    {
        minCost = min(minCost, LeetCode_64_Min_Path_Sum_Memoization(sr + 1, sc, grid, dp));
    }

    if (sc + 1 < grid[0].size())
    {
        minCost = min(minCost, LeetCode_64_Min_Path_Sum_Memoization(sr, sc + 1, grid, dp));
    }

    return dp[sr][sc] = minCost + grid[sr][sc];
}
int LeetCode_64_Min_Path_Sum_Tabulation(int sr, int sc, vector<vector<int>> &grid, vector<vector<int>> &dp)
{
    for(sr=grid.size()-1; sr>=0; sr--)
    {
        for(sc=grid[0].size()-1; sc>=0; sc--)
        {
            if (sr == grid.size() - 1 && sc == grid[0].size() - 1)
            {
                dp[sr][sc] = grid[sr][sc];
                continue;
            }

            int minCost = 1e8;
            if (sr + 1 < grid.size())
            {
                minCost = min(minCost, dp[sr + 1][sc]);
            }

            if (sc + 1 < grid[0].size())
            {
                minCost = min(minCost, dp[sr][sc+1]);
            }

            dp[sr][sc] = minCost + grid[sr][sc];
        }
    }

    return dp[0][0];
}
void LeetCode_64_Min_Path_Sum()
{
    vector<vector<int>> grid = {
        {1, 3, 1},
        {1, 5, 1},
        {4, 2, 1}};
    int n = grid.size(), m = grid[0].size();
    vector<vector<int>> dp(n, vector<int>(m, 0));
    cout << "Memoization :- " << LeetCode_64_Min_Path_Sum_Memoization(0, 0, grid, dp) << endl;
    display_2D(dp);
    cout << "Tabulation :- " << LeetCode_64_Min_Path_Sum_Tabulation(0, 0, grid, dp) << endl;
    display_2D(dp);
}

//***************************************************************
                  //LeetCode Problems//
//***************************************************************
void LeetCode_Problems()
{
    //LeetCode_70_ClimbStairs();

    //LeetCode_746_Min_Cost_Climbing_Stairs();

    //LeetCode_64_Min_Path_Sum();
}
//***************************************************************
                         //END//
//***************************************************************

void Geeks_For_Geeks()
{
    //Friends_Pairing_Problem();

    //GoldMine_Problem();

    //Count_Number_of_Ways(7, 3);   // int n=7, k=3
}

void DP_Basics_1()
{
    //MazePath_H_V_D();

    //MazePath_MultiMove();

    //Dice_Throw();

    //Dice_With_Array();
}

void DP_Basics_0()
{
    //Fibanocci_Series();
}

void DP_SubString_SubSequence_Set()
{
    //Leetcode_005_Longest_Pallindromic_SubString();

    //Leetcode_647_Count_All_Pallindromic_SubString();

    //Leetcode_516_Longest_Pallindromic_SubSequence();  

    LeetCode_115_Distinct_SubSequence();  
}
void solve()
{
    //DP_Basics_0();

    //DP_Basics_1();

    //LeetCode_Problems();

    //Geeks_For_Geeks();

    DP_SubString_SubSequence_Set();
}

int main(int argc, char **argv)
{
    solve();
    return 0;
}