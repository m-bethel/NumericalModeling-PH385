/*

Creating a 3-dimensional plot of charges and a disk held to 0.25V in a grid using the Jacobi method

Created by: Miles Bethel
02/02/2026
miles.d.bethel@gmail.com

*/


#include <iostream>
#include <chrono>


using namespace std;

int main()
{

    double c1;
    double c2;
    double c3;
    double c4;
    double V;
    char choice;
    // User prompts for initial conditions
    // This asks user if they would like to use default conditons provided
    cout << "Would you like to accept default initial conditions? [y,n]\n";
    cout << endl << "C1 as +1 uC at (25,0,0) \n"
                 << "C2 as +1 uC at (0,25,0) \n"
                 << "C3 as -1 uC at (-25,0,0) \n"
                 << "C4 as -1 uC at (0,-25,0) \n"
                 << "V = 0.25V" << endl;

    cin >> choice;

    cin.ignore(numeric_limits<streamsize>::max(), '\n');
    
    // If statement to check if you would like to use default
    // conditions or input your own
    if (choice  == 'y' || choice == 'Y')
    {
        cout << "Using default conditions" << endl;
        c1=1e-6,c2=1e-6,c3=-1e-6,c4=-1e-6,V=0.25;
        
    } else

    return 0;
}   