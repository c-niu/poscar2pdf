//
//  poscar2pdf.cpp
//  
//  Created by Changning Niu on 2/9/16.
//  Copyright © 2016 Changning Niu. All rights reserved.
//  For more information, go to https://github.com/changning/poscar2pdf

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <limits>

using namespace std;

// read_poscar:
//		return 1 if reading successful
//		return 2 if POSCAR has no line of elements
//		return 0 if reading failed
int read_poscar (char *file, string& head, double& dLatConst,
                 double dArrVec[3][3], string strArrElem[10], int nArrElem[10],
                 char& cCoord, double dArrAtom[][3], int& nElem, int& nAtom) {
    istringstream iss;
    ifstream poscar;
    string line;
    double x, y, z;
    int num;
    
    poscar.open(file);
    if (!poscar) return 0;
    getline(poscar,head);	// line 1
    getline(poscar,line);
    iss.str(line);			// line 2
    if (iss >> dLatConst) iss.clear();	// read lattice constant
    else return 0;
    for (int i = 0; i < 3; i++) {	// line 3-5: lattice vectors
        getline(poscar,line);
        iss.str(line);
        if (iss >> x >> y >> z) {
            iss.clear();
            dArrVec[i][0] = x;
            dArrVec[i][1] = y;
            dArrVec[i][2] = z;
        }
        else return 0;
    }
    getline(poscar,line);	// line 6
    iss.str(line);
    nAtom = 0;
    nElem = 0;
    if (iss >> num) {		// if line 6 has integers
        nAtom += num;
        while (iss >> num) {	// count numberf in line 6
            nElem++;
            nAtom += num;
        }
        string strArrMark[10] = {"A","B","C","D","E","F","G","H","I","J"};
        for (int i = 0; i < nElem; i++)     // Mark artificial element types
            strArrElem[i] = strArrMark[i];
    }
    else {                  // if line 6 has element types
        iss.clear();
        iss.str(line);
        num = 0;
        while (iss >> strArrElem[num]) {    // read elements
            num++;
            nElem++;
        }
        if (nElem == 0) return 0;       // if no element read, return 0
        getline(poscar,line);
        iss.clear();
        iss.str(line);
        for (int i = 0; i < nElem; i++) {   // read line 7: number for elements
            if (iss >> num) {
                nArrElem[i] = num;
                nAtom += num;
            }
            else return 0;
        }
    }
    getline(poscar,line);
    cCoord = line[0];               // read line 8 first character
    if (cCoord == 'S' || cCoord == 's') {   // if it's selective relaxation
        getline(poscar,line);               // ignore this line
        cCoord  = line[0];
    }
    for (int i = 0; i < nAtom; i++) {
        iss.clear();
        getline(poscar,line);
        iss.str(line);
        if (iss >> x >> y >> z) {
            dArrAtom[i][0] = x;
            dArrAtom[i][1] = y;
            dArrAtom[i][2] = z;
        }
        else return 0;
    }
    return 1;
}

// cart_to_direct:
// Converts atomic coordinates in Cartesian format to Directional format.
//    return 1 if successful, 0 if failed.
int cart_to_direct(double dArrVec[3][3], double dArrAtom[][3], int nAtom) {
    double det;
    double newAtom[nAtom][3];
    double rev[3][3];
    det = dArrVec[0][0] * dArrVec[1][1] * dArrVec[2][2]
        - dArrVec[0][0] * dArrVec[1][2] * dArrVec[2][1]
        - dArrVec[0][1] * dArrVec[1][0] * dArrVec[2][2]
        + dArrVec[0][1] * dArrVec[1][2] * dArrVec[2][0]
        + dArrVec[0][2] * dArrVec[1][0] * dArrVec[2][1]
        - dArrVec[0][2] * dArrVec[1][1] * dArrVec[2][0];
    rev[0][0] = (dArrVec[1][1] * dArrVec[2][2] - dArrVec[1][2] * dArrVec[2][1]);
    rev[0][1] = (dArrVec[0][2] * dArrVec[2][1] - dArrVec[0][1] * dArrVec[2][2]);
    rev[0][2] = (dArrVec[0][1] * dArrVec[1][2] - dArrVec[0][2] * dArrVec[1][1]);
    rev[1][0] = (dArrVec[1][2] * dArrVec[2][0] - dArrVec[1][0] * dArrVec[2][2]);
    rev[1][1] = (dArrVec[0][0] * dArrVec[2][2] - dArrVec[0][2] * dArrVec[2][0]);
    rev[1][2] = (dArrVec[0][2] * dArrVec[1][0] - dArrVec[0][0] * dArrVec[1][2]);
    rev[2][0] = (dArrVec[1][0] * dArrVec[2][1] - dArrVec[1][1] * dArrVec[2][0]);
    rev[2][1] = (dArrVec[0][1] * dArrVec[2][0] - dArrVec[0][0] * dArrVec[2][1]);
    rev[2][2] = (dArrVec[0][0] * dArrVec[1][1] - dArrVec[0][1] * dArrVec[1][0]);
    for (int i = 0; i < nAtom; i++) {
        newAtom[i][0] = dArrAtom[i][0] * rev[0][0] + dArrAtom[i][1] * rev[1][0] + dArrAtom[i][2] * rev[2][0];
        newAtom[i][1] = dArrAtom[i][0] * rev[0][1] + dArrAtom[i][1] * rev[1][1] + dArrAtom[i][2] * rev[2][1];
        newAtom[i][2] = dArrAtom[i][0] * rev[0][2] + dArrAtom[i][1] * rev[1][2] + dArrAtom[i][2] * rev[2][2];
    }
    for (int i = 0; i < nAtom; i++) {
        dArrAtom[i][0] = newAtom[i][0] / det;
        dArrAtom[i][1] = newAtom[i][1] / det;
        dArrAtom[i][2] = newAtom[i][2] / det;
    }
    return 1;
}

// Bubble sort
// Based on codes from http://www.algolist.net/Algorithms/Sorting/Bubble_sort
void bubbleSort(double arr[], int n) {
    bool swapped = true;
    int j = 0;
    double tmp;
    while (swapped) {
        swapped = false;
        j++;
        for (int i = 0; i < n - j; i++) {
            if (arr[i] > arr[i + 1]) {
                tmp = arr[i];
                arr[i] = arr[i + 1];
                arr[i + 1] = tmp;
                swapped = true;
            }
        }
    }
}

int main(int argc, char* argv[]) {
    double dMax, dMin, dStep;   // min/max and step of desired PDF radius
    double dMax2, dMin2;        // squared min/max for acceleration
    char *file;                 // pointer to POSCAR
    int nElem;                  // number of elements
    int nAtom;                  // number of atoms
    string head;                // line 1 of POSCAR
    double dLatConst;           // line 2 of POSCAR
    double dArrVec[3][3];       // line 3-5: lattice vectors
    string strArrElem[10];      // element names
    int nArrElem[10];           // number of atoms of each element
    char cCoord;                // line 8: directional or cartesian
    double dArrAtom[500][3];    // atomic coordinates
    double dArrAtom2[500][3];   // atomic coordinates in Cartesian
    string strArrAtom[500];     // element of each atom
    int count;
    double A1, B1, C1, A2, B2, C2, A3, B3, C3;// coefficients A_i, B_i, C_i for finding min d
    double d1, d2, d3, dMinCell;// distances between two opposite sides of cell
    int nCell;                  // number of cells outside each side to build supercell
    int nPairs;                 // number of pairs of elements
    int nData;                  // number of data points of PDF
    int nSCatom;                // number of atoms in the supercell
    double dTemp;               // temperary usage
    
    // Check input
    if (argc != 5) {
        cout << "Usage: poscar2pdf POSCAR 0.1 0.1 4.0\n";
        exit(1);
    }
    file = argv[1];
    dMin = atof(argv[2]);
    dStep = atof(argv[3]);
    dMax = atof(argv[4]);
    if (dMin == 0.0) {
        cout << "Min value must be larger than 0.\n";
        exit(1);
    }
    // Read POSCAR
    if (read_poscar(file, head, dLatConst, dArrVec, strArrElem, nArrElem, cCoord, dArrAtom, nElem, nAtom) == 0) {
        cout << "Errors occurred when reading POSCAR.\n";
        exit(1);
    }
    else if (cCoord == 'C' && cCoord == 'c') {
        cart_to_direct(dArrVec, dArrAtom, nAtom);
    }
    
    // Get atomic species for each atom
    count = 0;
    for (int i = 0; i < nElem; i++) {
        for (int j = 0; j < nArrElem[i]; j++) {
            strArrAtom[count++] = strArrElem[i];
        }
    }
    
    // If (dMax - dMin) / dStep is not integer, decrease dMax slightly
    nData = (int)((dMax - dMin) / dStep) + 1;
    dMax = (nData - 1) * dStep + dMin;
    dMin2 = dMin * dMin / dLatConst / dLatConst;
    dMax2 = dMax * dMax / dLatConst / dLatConst;
    if (nData > 10000) {
        cout << "Number of PDF data points exceeds 10000. Change the code.\n";
        exit(1);
    }
    
    // Determine the size of supercell and expand by Cartesian
    A1 = dArrVec[1][1] * dArrVec[2][2] - dArrVec[2][1] * dArrVec[1][2];
    B1 = dArrVec[2][1] * dArrVec[0][2] - dArrVec[0][1] * dArrVec[2][2];
    C1 = dArrVec[0][1] * dArrVec[1][2] - dArrVec[1][1] * dArrVec[0][2];
    A2 = dArrVec[1][0] * dArrVec[2][2] - dArrVec[2][0] * dArrVec[1][2];
    B2 = dArrVec[2][0] * dArrVec[0][2] - dArrVec[0][0] * dArrVec[2][2];
    C2 = dArrVec[0][0] * dArrVec[1][2] - dArrVec[1][0] * dArrVec[0][2];
    A3 = dArrVec[1][0] * dArrVec[2][1] - dArrVec[2][0] * dArrVec[1][1];
    B3 = dArrVec[2][0] * dArrVec[0][1] - dArrVec[0][0] * dArrVec[2][1];
    C3 = dArrVec[0][0] * dArrVec[1][1] - dArrVec[1][0] * dArrVec[0][1];
    d1 = abs(A1 * dArrVec[0][0] + B1 * dArrVec[0][1] + C1 * dArrVec[0][2]);
    d2 = abs(A2 * dArrVec[1][0] + B2 * dArrVec[1][1] + C2 * dArrVec[1][2]);
    d3 = abs(A3 * dArrVec[2][0] + B3 * dArrVec[2][1] + C3 * dArrVec[2][2]);
    dMinCell = fmin(d1, fmin(d2, d3)) * dLatConst;
    nCell = (int)(dMax / dMinCell + 1);
    nSCatom = nAtom * (int)(pow((nCell * 2 + 1), 3));
    if (nSCatom > 100000) {
        cout << "Supercell needs over 100000 atoms. Change the code to suit your calculation.\n";
        exit(1);
    }
    double dArrSC[100000][3];          // declare array containing supercell coordinates
    double dArrSCr[100000];            // declare array to hold temp radii from one center atom
    string strArrSC[100000];           // declare array containing elements in supercell
    count = 0;
    for (int i = -nCell; i <= nCell; i++) {
        for (int j = -nCell; j <= nCell; j++) {
            for (int k = -nCell; k <= nCell; k++) {
                for (int m = 0; m < nAtom; m++) {
                    dArrSC[count][0] = (dArrAtom[m][0] + i) * dArrVec[0][0]
                                     + (dArrAtom[m][1] + j) * dArrVec[1][0]
                                     + (dArrAtom[m][2] + k) * dArrVec[2][0];
                    dArrSC[count][1] = (dArrAtom[m][0] + i) * dArrVec[0][1]
                                     + (dArrAtom[m][1] + j) * dArrVec[1][1]
                                     + (dArrAtom[m][2] + k) * dArrVec[2][1];
                    dArrSC[count][2] = (dArrAtom[m][0] + i) * dArrVec[0][2]
                                     + (dArrAtom[m][1] + j) * dArrVec[1][2]
                                     + (dArrAtom[m][2] + k) * dArrVec[2][2];
                    strArrSC[count] = strArrAtom[m];
                    count++;
                }
            }
        }
    }
    if (count != nSCatom) {
        cout << "Error when generating supercell.\n";
        exit(1);
    }
    
    // Convert atomic coordinates to Cartesian
    for (int i = 0; i < nAtom; i++) {
        dArrAtom2[i][0] = dArrAtom[i][0] * dArrVec[0][0]
                        + dArrAtom[i][1] * dArrVec[1][0]
                        + dArrAtom[i][2] * dArrVec[2][0];
        dArrAtom2[i][1] = dArrAtom[i][0] * dArrVec[0][1]
                        + dArrAtom[i][1] * dArrVec[1][1]
                        + dArrAtom[i][2] * dArrVec[2][1];
        dArrAtom2[i][2] = dArrAtom[i][0] * dArrVec[0][2]
                        + dArrAtom[i][1] * dArrVec[1][2]
                        + dArrAtom[i][2] * dArrVec[2][2];
    }
    
    // Calculate pair distribution function
    nPairs = (nElem + 1) * nElem / 2;
    if (nPairs > 20) {
        cout << "Pair of elements exceeds 20. Change the code.\n";
        exit(1);
    }
    double dArrPDF[10000][20]; // declare the array containing PDF data
    for (int i = 0; i < nData; i++) {
        dArrPDF[i][0] = dMin + i * dStep;   // first column is r
        for (int j = 1; j <= nPairs + 1; j++) {
            dArrPDF[i][j] = 0.0;        // initialize PDF data
        }
    }
    int count2 = 0; // used to count pairs
    for (int i = 0; i < nElem; i++) {   // This and next loop scan all pairs
        for (int j = 0; j <= i; j++) {
            count2++;
            for (int atom = 0; atom < nAtom; atom++) {  // This loop scan all atoms as center atom
                if (strArrAtom[atom] == strArrElem[i]) {
                    count = 0;          // used to count number of atoms of supercell in the range
                    for (int SCatom = 0; SCatom < nSCatom; SCatom++) {// This loop scan all supercell atoms
                        if (strArrSC[SCatom] == strArrElem[j]) { // if this is the pair of elements we want
                            dTemp = pow((dArrSC[SCatom][0]-dArrAtom2[atom][0]), 2)
                                  + pow((dArrSC[SCatom][1]-dArrAtom2[atom][1]), 2)
                                  + pow((dArrSC[SCatom][2]-dArrAtom2[atom][2]), 2);
                            if (dTemp <= dMax2 && dTemp >= dMin2) {
                                dArrSCr[count] = sqrt(dTemp) * dLatConst;
                                count++;
                            }
                        }
                    }
                    bubbleSort(dArrSCr, count);
                    int num = 0;
                    for (int m = 0; m < count; m++) {   // check each recorded radius
                        for (int n = num; n < nData; n++) {
                            if (dArrSCr[m] >= dArrPDF[n][0] && dArrSCr[m] < dArrPDF[n][0] + dStep) {
                                dArrPDF[n][count2] += 1.0;
                                num = n;
                                break;
                            }
                        }
                    }
                }
            }
            for (int k = 0; k < nData; k++) {   // divide counted integers by volume and number of center atom
                dArrPDF[k][count2] /= (4*3.14159*pow(dArrPDF[k][0],2)*dStep*nArrElem[i]);
            }
        }
    }
    
    // Calculate total PDF
    for (int i = 0; i < nAtom; i++) {
        count = 0;
        for (int j = 0; j < nSCatom; j++) {
            dTemp = pow((dArrSC[j][0]-dArrAtom2[i][0]), 2)
                  + pow((dArrSC[j][1]-dArrAtom2[i][1]), 2)
                  + pow((dArrSC[j][2]-dArrAtom2[i][2]), 2);
            if (dTemp <= dMax2 && dTemp >= dMin2) {
                dArrSCr[count] = sqrt(dTemp) * dLatConst;
                count++;
            }
        }
        bubbleSort(dArrSCr, count);
        int num = 0;
        for (int m = 0; m < count; m++) {
            for (int n = num; n < nData; n++) {
                if (dArrSCr[m] >= dArrPDF[n][0] && dArrSCr[m] < dArrPDF[n][0] + dStep) {
                    dArrPDF[n][nPairs+1] += 1.0;
                    num = n;
                    break;
                }
            }
        }
    }
    for (int k = 0; k < nData; k++) {
        dArrPDF[k][nPairs+1] /= (4*3.14159*pow(dArrPDF[k][0],2)*dStep*nAtom);
    }
    
    // Output
    ofstream os;
    char output[] = "PDF_";
    strcat(output, argv[1]);
    strcat(output, ".txt");
    os.open(output);
    os << "  r  ";
    for (int i = 0; i < nElem; i++) {
        for (int j = 0; j <= i; j++) {
            os << setw(5) << strArrElem[i];
            os << "-";
            os << setw(2) << strArrElem[j];
        }
    }
    os << "  Total "<< "\n";
    for (int i = 0; i < nData; i++) {
        os << setw(5) << fixed << setprecision(3) << dArrPDF[i][0];
        for (int j = 1; j <= nPairs+1; j++) {
            os << setw(8) << fixed << setprecision(4) << dArrPDF[i][j];
        }
        os << "\n";
    }
    os.close();
    return 0;
}
