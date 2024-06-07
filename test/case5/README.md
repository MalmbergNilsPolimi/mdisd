# Test Case 5
For more details -> <a href="https://github.com/MalmbergNilsPolimi/mdisd/blob/main/doc/mdisdReport.pdf" target="_blank">mdisdReport.pdf</a>


To use the test case, first modify the path to Eigen library in the Makefile. 

Then use :
```
make && make run
```

To clean the repository, use:
```
make clean
````

The user can choose if the interpolation is done on the 2D or 3D Franke function by changing lines 260 and 261. Moreover, the number of known points to test can be modify on lines 263 and 264. In the line 268, the user can choose how many points to interpolate for each interpolation.