#include<iostream>

using namespace std;

void Print_Increasing(int a, int b)
{
    if(a == b+1)
    {
        return;
    }
    cout << a << " ";
    Print_Increasing(a+1, b);
}
void Print_Decreasing(int a, int b)
{
    if(a == b+1)
    {
        return;
    }
    
    Print_Decreasing(a+1, b);
    cout << a << " ";
}
void Print_Increasing_Decreasing_Even_Odd(int a, int b)
{
    if(a == b+1)
    {
        return;
    }
    
    if(a%2 == 0)
    {
        cout << a << " Hi " << endl;
    }
    Print_Increasing_Decreasing_Even_Odd(a+1, b);
    if(a%2 != 0)
    {
        cout << a << " Bye " << endl;
    }
}
void Recursion_Basic()
{
    //Print_Increasing(0, 10);
    //Print_Decreasing(0, 10);

    //Print_Increasing_Decreasing(0, 10);
    //Print_Increasing_Decreasing_Even_Odd(0, 10);
}
void solve()
{
    Recursion_Basic();
}
int main(int argc, char** argv)
{
    solve();
    return 0;
}