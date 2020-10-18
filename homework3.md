## Question 1

if we have a folder directory structure of:  

~/dir1/dir2/dir3/dir4/dir5  

if we are in home, and use the tree command, the structure looks like so:  

```
.
+--- dir1
|   +--- dir2
|   |   +--- dir3
|   |   |   +--- dir4
|   |   |   |   +--- dir5
|   |   |   +--- folder3

```
 
If we are currently sitting in dir5, show the commands needed to move into dir3, delete folder3, then move into dir2 and create a folder called folder2

Answer:  

```
cd ../..  
rmdir folder3  
cd ..  
mkdir folder2  
```

## Question 2

Suppose you have a matrix, A, and a dataframe, B, in R.  

matrix A:  
```
     [,1] [,2] [,3]
[1,]   11   22   33
[2,]   44   55   66
[3,]   77   88   99
```
dataframe B:
```
  id    name      major
1  1   sally    biology
2  2 freddie    physics
3  3    jack statistics
```

What command would you use to get the following output:  
```
[1] 11 44 77
```
What two different commands could you use to get the following output:  
```
[1] 1 2 3
```
What command would you use to get the following output:  
```
  id
1  1
2  2
3  3
```
Answers:  

1) A[,1]  
2) B[,1] or B[[1]]  
3) B[1]  

## Question 3

Create a simple script called `scripty.sh` that assigns a value of "name" to the variable namey, then uses echo to print the content of the variable namey to the console

Answer:

```
#!/bin/bash

namey="name"
echo $namey
```

Now answer the following question:

Which of the following is a way, **without** using octal permissions, to provide `scripty.sh` executable access to everyone, but only writeable access to the owner (you):

* A) `chmod u+rwx,go+rx scripty.sh`  
* B) `chmod u-r,g+rwx,o-rx scripty.sh`  
* C) `chmod 777 scripty.sh`  
* D) `chmod ugo+rwx scripty.sh`  

Answer: A
