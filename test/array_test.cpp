#include "array.hpp"
#include <iostream>
#include <cassert>

using namespace std;


void test_2D()
{
    int n1=2, n2=3;
    Array2D<double> a(n1, n2);

    a(0, 1) = 3.1;
    a(1, 2) = -9.2;

    double* b = a.data();

    cout << "i  j  a(i, j)  b[count]" << endl;
    {
        int count = 0;
        for(int i=0; i<n1; i++){
            for(int j=0; j<n2; j++){
                cout << i << " " << j << " " << a(i, j) << " " << b[count] << endl;
                assert(a(i, j) == b[count]);
                count++;
            }
        }
    }

    a.zeros();

    cout << "\nset zeros" << endl;
    {
        int count = 0;
        for(int i=0; i<n1; i++){
            for(int j=0; j<n2; j++){
                cout << i << " " << j << "  " << a(i, j) << " " << b[count] << endl;
                assert(a(i, j) == 0);
                count++;
            }
        }
    }

}

void test_3D()
{
    int n1=2, n2=3, n3=4;
    Array3D<double> a(n1, n2, n3);

    a(0, 1, 0) = 3.1;
    a(1, 2, 3) = -9.2;

    double* b = a.data();

    cout << "i  j  k  a(i, j, k)  b[count]" << endl;
    {
        int count = 0;
        for(int i=0; i<n1; i++){
            for(int j=0; j<n2; j++){
                for(int k=0; k<n3; k++){
                    cout << i << " " << j << " " << k << "  " << a(i, j, k) << " " << b[count] << endl;
                    assert(a(i, j, k) == b[count]);
                    count++;
                }
            }
        }
    }
}

void test_4D()
{
    int n1=2, n2=3, n3=4, n4=2;
    Array4D<double> a(n1, n2, n3, n4);

    a(0, 1, 0, 1) = 3.1;
    a(1, 2, 3, 0) = -9.2;
    a(0, 0, 1, 1) = 2.2;

    double* b = a.data();

    cout << "i  j  k  l  a(i, j, k, l)  b[count]" << endl;
    {
        int count = 0;
        for(int i=0; i<n1; i++){
            for(int j=0; j<n2; j++){
                for(int k=0; k<n3; k++){
                    for(int l=0; l<n4; l++){
                        cout << i << " " << j << " " << k << " " << l << "  " << a(i, j, k, l) << " " << b[count] << endl;
                        assert(a(i, j, k, l) == b[count]);
                        count++;
                    }
                }
            }
        }
    }
}

void test_5D()
{
    int n1=2, n2=3, n3=4, n4=2, n5=3;
    Array5D<double> a(n1, n2, n3, n4, n5);

    a(0, 1, 0, 1, 0) = 3.1;
    a(1, 2, 3, 0, 2) = -9.2;
    a(0, 0, 1, 1, 1) = 2.2;
    a(1, 0, 2, 1, 0) = -4.4;

    double* b = a.data();

    cout << "i  j  k  l  m  a(i, j, k, l, m)  b[count]" << endl;
    {
        int count = 0;
        for(int i=0; i<n1; i++){
            for(int j=0; j<n2; j++){
                for(int k=0; k<n3; k++){
                    for(int l=0; l<n4; l++){
                        for(int m=0; m<n5; m++){
                            cout << i << " " << j << " " << k << " " << l << " " << m << "  " << a(i, j, k, l, m) << " " << b[count] << endl;
                            assert(a(i, j, k, l, m) == b[count]);
                            count++;
                        }
                    }
                }
            }
        }
    }
}

main()
{
    cout << "\nTest for Array2D" << endl;
    test_2D();

    cout << "\nTest for Array3D" << endl;
    test_3D();

    cout << "\nTest for Array4D" << endl;
    test_4D();

    cout << "\nTest for Array5D" << endl;
    test_5D();
}
