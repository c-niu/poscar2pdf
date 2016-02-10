# poscar2pdf
This code calculates the pair distribution function (PDF) from a POSCAR.
## Installation and Usage
1. Compile the c++ code: `c++ poscar2pdf.cpp -o poscar2pdf`
2. Run the code: `poscar2pdf POSCAR 0.1 0.05 5.0` (where 0.1 and 5.0 set the range of radius, and 0.05 is the size of steps).
3. Plot PDF using whatever visualization software you like with the data file `PDF_POSCAR.txt`. The first row has the titles of each column, which from left to right are the radius, PDFs of pair of elements, and total PDF.

_A small tip_

I find using a gnuplot script to plot the PDF is convenient, especially when I need to tune the size of steps to get a better looking PDF figure. In this case, all I need to do is change one variable in a bash script and re-run this script to get the .EPS file. An example is provided in the subdirectory `gnuplot`.
## Maths In The Code
### How to calculate the PDF
To calculate the PDF, or g(r), the following is done by the code:

1. Read the input for *r*: *min*, *dr*, *max*.
2. Expand the POSCAR so that the sphere defined by *max* is entirely covered by the expanded supercell.
3. There are PDF for any two elements and PDF of the whole structure, which are calculated in different ways.
	1. For the total PDF, consider every atom in the POSCAR, one at a time. Find all atoms in the supercell that are at a distance between *min* and *max* from the considered atom. Count the number of atoms that are between *r* and *r+dr*, where *r = min + n dr*. Repeat for every atom in the POSCAR. Then for each step *n* divide this number by the volume of that thin layer of shell, *4 pi r^2 dr*, and by the number of atoms in the POSCAR. This is g(r) for the entire POSCAR.
	2. For the PDF of two elements, consider every atom of one element in the POSCAR, one at a time. Find all atoms of the other element in the supercell that are at a distance between *min* and *max* from the considered atom. Count the number of atoms that are between *r* and *r+dr*, where *r = min + n dr*. Repeat for every atom of the first element in the POSCAR. Then for each step *n* divide this number by the volume of that thin layer of shell, *4 pi r^2 dr*, and by the number of atoms of the first element in the POSCAR. This is g(r) for the two elements. Repeat for other pair of elements.

### How to determine the size of supercell

POSCAR only contains a single periodic unit that represents the entire bulk material. We must expand this periodic unit cell to a supercell which is big enough to calculate the PDF within the desired range. So here comes the question: how big is the expanded supercell?

A 3-3-3 supercell has the original periodic unit cell in the center and eight cells surrounding it. The maximum range of the PDF is the radius of the largest sphere inside the supercell whose spherical center is anywhere in the periodic unit cell in the center. Therefore, this radius equals the smallest distance between any two opposite sides of the periodic unit cell. This is only true for the 3-3-3 supercell. If the desired maximum range of PDF is larger than this radius, a supercell larger than 3-3-3 is needed. Let’s put it this way:

<img src="Image-1.png" width="200" align="middle">

where *r\_max* is the maximum range of the desired PDF, *n* is the thickness of periodic unit cell surrounding the one in the center, and *d\_min* is the smallest distance between any two opposite sides of the periodic unit cell.
Now, we first need to calculate *d\_min*. The periodic unit cell is defined by three vectors, or three points plus the origin. The distance between two opposite sides equals the distance from one point to the plane defined by the other two points and the origin. This is a well-defined math problem. Let’s skip to the answers:

![](Image-2.png)

where *n\_1* is the normal vector to the plane defined by vector 2 and 3, *d\_1* is the distance from point 1 to the plane defined by vector 2 and 3, etc.
