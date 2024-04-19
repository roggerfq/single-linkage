# single-linkage++

## Description
This algorithm constructs the hierarchical cluster of an n-dimensional dataset using the single-linkage algorithm described in [1].

## Running an example
Compile the file:
``` shell
g++ -o single_linkage single_linkage.cpp
```
Then, execute the compiled program:
``` shell
./single_linkage
```
This will display the following output:
```plaintext
Data:
A: 1, 1 
B: 1.5, 1.5 
C: 2, 2 
D: 1, 1 
E: 9, 9 
F: 1, 1 
G: 7, 7 
H: 9, 11 

Hierarchical cluster:
A,D,F,B,C,E,H,G: 3.9375, 4.1875 
               A,D,F,B,C: 1.3, 1.3 
                        A,D,F,B: 1.125, 1.125 
                               A,D,F: 1, 1 
                                    A,D: 1, 1 
                                       A: 1, 1 
                                       D: 1, 1 
                                    F: 1, 1 
                               B: 1.5, 1.5 
                        C: 2, 2 
               E,H,G: 8.33333, 9 
                    E,H: 9, 10 
                       E: 9, 9 
                       H: 9, 11 
                    G: 7, 7 
```
The result shown in the terminal could be interpreted as the following dendrogram:
![alt text](https://raw.githubusercontent.com/roggerfq/single-linkage/main/example.png)

## Using code in your own projects
The data shown above corresponds to the example of the function **check_example**. In this function, you can see how to add data using the class **Cluster** and how to construct the hierarchical cluster through the function **clustering**. The following is an example of how to add a single data point named "label_name" with a centroid of (1.3, 2.4):

```cpp
const int dim = 2;
std::set < Cluster * > clusters;  

// Add new data
{
 Cluster * new_cluster = new Cluster();
 // Configure new cluster
 new_cluster -> set_label("C");
 float point[dim] = {1.3, 2.4};
 new_cluster -> set_centroid(point, dim);
 clusters.insert(new_cluster); // Adding new data
}

```

After adding the data, we can call the function **clustering** to construct the hierarchical cluster:
```cpp
Cluster * big_cluster = clustering(clusters);
big_cluster -> show();// Show hierarchical cluster in terminal
```

The signature of the function **clustering** is:
```cpp
Cluster * clustering(const std::set < Cluster * > & clusters, norm_function norm = norml2)
```
Notice that by default, the norm function used to compute the distance between clusters is the l2-norm (refer to the source code file **norml2**). However, this behavior can be changed by implementing another function of type **norm_function**:
```cpp
typedef
const float( * norm_function)(Cluster * a, Cluster * b);

// l2-norm
const float norml2(Cluster * a, Cluster * b) {
  const float * centroid_a = a -> get_centroid();
  const float * centroid_b = b -> get_centroid();
  int dim_centroid = a -> get_dim_centroid();
  float dot_product = 0;
  
  for (int i = 0; i < dim_centroid; i++) {
    float diff = centroid_a[i] - centroid_b[i];
    dot_product += diff * diff;
  }
  
  return sqrt(dot_product);
}
```
## Author
Roger Figueroa Quintero - [Contact LinkedIn](https://www.linkedin.com/in/roger-figueroa-quintero/)

## License
This project is licensed under the [MIT License](LICENSE.md), which permits unrestricted use, modification, and distribution subject to the terms and conditions of the license.

## References
[1] Anil K. Jain and Richard C. Dubes. 1988. Algorithms for clustering data. Prentice-Hall, Inc., USA
