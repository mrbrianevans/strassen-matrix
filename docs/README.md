 # The Strassen Method
This is a method of matrix multiplication which is more computationally efficient than the regular method. <br>
It can be used to multiply two square matrices of size nxn where n is a power of two.

### In the source code:
The [Matrix.java](https://github.com/mrbrianevans/strassen-matrix/blob/master/Matrix.java) file is my Java implementation.<br>
It contains:
- A matrix object (attributes and constructors)
- Matrix addition method
- Matrix subtraction method
- .equals() method
- .toString() formatting method for matrices
- Standard matrix dot product
- Strassen method of matrix multiplication method
### The Strassen method
The Strassen method splits the two matrices up into four quadrants recursivley, until it reaches a 2x2 matrix (the base case).
