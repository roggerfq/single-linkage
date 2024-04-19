// Copyright (c) 2020
// Written by Roger Figueroa Quintero
  
// This single-linkage algorithm uses the method described in [1]. Here, the clustering is treated as a graph whose nodes correspond to individual clusters. Then, we calculate the distance matrix and sort it in decreasing order. Finally, we connect the nodes beginning with the minimum distance G(0) until we obtain a graph where each node has one link, denoted as G(k). Note that if there are multiple nodes with equal distances, they are grouped in pairs; that is, at each iteration k, we only merge two clusters. This is different from the example shown in Wikipedia [2], but both MATLAB [3] and other online resources [4] implement it this way.  

// [1] Anil K. Jain and Richard C. Dubes. 1988. Algorithms for clustering data. Prentice-Hall, Inc., USA.
// [2] https://en.wikipedia.org/wiki/Single-linkage_clustering
// [3] https://www.mathworks.com/help/stats/linkage.html 
// [4] https://people.revoledu.com/kardi/tutorial/Clustering/Online-Hierarchical-Clustering.html


#include <iostream>
#include <stdio.h>
#include <string.h>
#include <algorithm>    // std::sort
#include <vector>       // std::vector
#include <math.h>
#include <set>


/*
 
You have the flexibility to implement the norm function of your choice by utilizing the following function pointer type:

     typedef const float (*norm_function)(Cluster *a, Cluster *b);

and pass it as an argument before calling:

    Cluster *clustering(const std::set<Cluster *> &clusters, norm_function norm = norml2)

*/

class Cluster {

  Cluster * parent;
  Cluster * big_cluster;
  std::vector < Cluster * > children;
  std::vector < Cluster * > terminal_clusters;
  std::vector < Cluster * > all_partitions;
  std::string label;
  float * centroid;
  int dim_centroid;
  bool flag_terminal;

  int pos_cal_centroid;

  Cluster * get_new_big_cluster() {

    if (big_cluster -> parent == NULL)
      return big_cluster;
    else
      big_cluster = big_cluster -> parent -> get_new_big_cluster();

    return big_cluster;

  }

  float * get_new_centroid() {

    if (isTerminal() || (dim_centroid == 0) || (centroid == NULL))
      return centroid;

    int N = get_num_children();

    if (pos_cal_centroid == N)
      return centroid;

    float nt = float(get_size_terminal_clusters());
 

    for (int i = pos_cal_centroid; i < N; i++) {

      Cluster * tmp_child = children[i];
      float * tmp_centroid = tmp_child -> get_new_centroid();
      int ntc = tmp_child -> get_size_terminal_clusters();

      if (tmp_child -> isTerminal())
        ntc = 1;

      float mc = (i == pos_cal_centroid) ? float(pos_cal_centroid) / nt : 1.0;
      for (int j = 0; j < dim_centroid; j++) {

        centroid[j] = centroid[j] * mc + (tmp_centroid[j] * float(ntc) / nt);

      }

    }

    pos_cal_centroid = N;

    return centroid;

  }

  public:

    Cluster(): parent(NULL), label(""), flag_terminal(true), pos_cal_centroid(0) {

      dim_centroid = 0;
      centroid = NULL;
      big_cluster = this;
      terminal_clusters.push_back(this);

    }

  Cluster(const float * var_centroid,
      const int var_dim_centroid, std::string label = ""): parent(NULL), label(label), flag_terminal(true), pos_cal_centroid(0) {

      // Here we construct a terminal cluster
      dim_centroid = 0;
      centroid = NULL;
      set_centroid(var_centroid, var_dim_centroid);

      big_cluster = this;
      terminal_clusters.push_back(this);

    }

    ~Cluster() {


      if (centroid != NULL)
        delete[] centroid;

      for (int i = 0; i < children.size(); i++) {
        delete children[i];
        children[i] = NULL;
      }

    }

  void set_centroid(const float * var_centroid,
    const int var_dim_centroid) {

    if (flag_terminal) {

      dim_centroid = var_dim_centroid;

      if (centroid != NULL)
        delete[] centroid;

      centroid = new float[dim_centroid];
      memcpy(centroid, var_centroid, dim_centroid * sizeof(float));

    }

  }

  void add_cluster(Cluster * cluster) {

    // If we call this method then the cluster becomes in terminal	  

    if (flag_terminal) {

      dim_centroid = cluster -> get_dim_centroid();
      if (centroid != NULL)
        delete[] centroid;

      centroid = new float[dim_centroid];
      std::fill(centroid, centroid + dim_centroid, 0.0);

      flag_terminal = false;
      terminal_clusters.clear();

    }

    cluster -> parent = this;

    const std::vector < Cluster * > & new_terminal_clusters = cluster -> get_terminal_clusters();

    for (int i = 0; i < cluster -> get_size_terminal_clusters(); i++)
      terminal_clusters.push_back(new_terminal_clusters[i]);

    children.push_back(cluster);

  }

  const Cluster * get_parent() const {

    return parent;
  }

  Cluster * get_big_cluster() {

    return get_new_big_cluster();

  }

  const float * get_centroid() {

    return get_new_centroid();

  }

  const int get_dim_centroid() const {

    return dim_centroid;

  }

  const std::vector < Cluster * > get_children() const {

    return children;

  }

  const int get_num_children() const {

    return children.size();

  }

  const std::vector < Cluster * > & get_terminal_clusters() const {

    return terminal_clusters;

  }

  const std::vector < Cluster * > & get_all_partitions(bool flag_compute = true) {

    if (flag_compute == false)
      return all_partitions;

    all_partitions.clear();
    all_partitions.push_back(this);

    for (int i = 0; i < children.size(); i++) {

      const std::vector < Cluster * > & child_partitions = children[i] -> get_all_partitions(flag_compute);

      for (int j = 0; j < child_partitions.size(); j++)
        all_partitions.push_back(child_partitions[j]);

    }

    return all_partitions;

  }

  const int get_size_terminal_clusters() const {

    return terminal_clusters.size();

  }

  const std::string get_label() const {

    return label;

  }

  void set_label(std::string label) {

    this -> label = label;

  }

  const bool isTerminal() const {

    return flag_terminal;

  }

  void show_terminal_clusters() const {

    for (int i = 0; i < get_size_terminal_clusters(); i++)
      std::cout << terminal_clusters[i] -> get_label() << " "<<terminal_clusters[i] ->get_parent()->get_label()<<" \n";
    std::cout << "\n";

  }

  void show_centroid() {

    const float * tmp_centroid = get_centroid();

    
    std::cout << tmp_centroid[0];
    for (int i = 1; i < dim_centroid; i++)
      std::cout << ", "<< tmp_centroid[i] << " ";
    std::cout << "\n";

  }

  void show(int pos = 0) {

    std::cout << std::string(pos, ' ') << label << ": ";
    show_centroid();

    pos += label.length();

    for (int i = 0; i < children.size(); i++)
      children[i] -> show(pos);

  }

};

void fill_random(float * vector, int dim) {

  for (int i = 0; i < dim; i++)
    vector[i] = static_cast < float > (rand()) / static_cast < float > (RAND_MAX);

}

typedef
const float( * norm_function)(Cluster * a, Cluster * b);

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

struct PairDist {

  float dist;
  Cluster * a;
  Cluster * b;

};

bool comp_dist(const PairDist * a,
  const PairDist * b) {

  return ((a -> dist) < (b -> dist));

}

Cluster * clustering(const std::set < Cluster * > & clusters, norm_function norm = norml2) {

  // We construct the distance matrix

  int n = clusters.size();
  int numDist = n * (n - 1) / 2;

  PairDist * set_dist = new PairDist[numDist];

  int cnt = 0;
  std::set < Cluster * > ::iterator it_i = clusters.begin();

  for (int i = 0; i < n - 1; i++) {

    Cluster * cluster_i = * it_i;

    std::set < Cluster * > ::iterator it_j = clusters.begin();
    std::advance(it_j, i + 1);

    for (int j = i + 1; j < n; j++) {

      Cluster * cluster_j = * it_j;

      (set_dist + cnt) -> a = cluster_i;
      (set_dist + cnt) -> b = cluster_j;
      (set_dist + cnt) -> dist = norm(cluster_i, cluster_j);

      cnt++;
      it_j++;

    }

    it_i++;

  }

  PairDist ** ptr_dist = new PairDist * [numDist];

  for (int i = 0; i < numDist; i++)
    ptr_dist[i] = (set_dist + i);

  std::vector < PairDist * > vec_dist(ptr_dist, ptr_dist + numDist);

  // We order distances in decreasing order
  std::sort(vec_dist.begin(), vec_dist.end(), comp_dist);


  for (int i = 0; i < vec_dist.size(); i++) {

    Cluster * big_cluster_a = (vec_dist[i]) -> a -> get_big_cluster();
    Cluster * big_cluster_b = (vec_dist[i]) -> b -> get_big_cluster();

    if (big_cluster_a != big_cluster_b) {

      Cluster * merge_cluster = new Cluster();

      merge_cluster -> set_label(big_cluster_a -> get_label() + "," + big_cluster_b -> get_label());
      merge_cluster -> add_cluster(big_cluster_a);
      merge_cluster -> add_cluster(big_cluster_b);

    } else {

      if (big_cluster_a -> get_size_terminal_clusters() == n)
        break;

    }

  }

  // Delete heap
  vec_dist.clear();
  delete[] set_dist;
  delete[] ptr_dist;

  // At this point, any node has the same root node as its big_cluster.
  return ( * (clusters.begin())) -> get_big_cluster();

}


// Random example
void check_random_example() {

  std::srand(time(NULL));
  const int n = 10;
  const int dim = 5;

  float point[dim];

  std::set < Cluster * > clusters;

  for (int i = 0; i < n; i++) {

    Cluster * new_cluster = new Cluster();

    // Configure new data
    new_cluster -> set_label(std::to_string(i));
    fill_random(point, dim);
    new_cluster -> set_centroid(point, dim);

    clusters.insert(new_cluster);

  }

  std::cout << "Data:\n";

  std::set < Cluster * > ::iterator it = clusters.begin();
  for (int i = 0; i < clusters.size(); i++) {
    ( * it) -> show();
    std::advance(it, 1);
  }

  Cluster * big_cluster = clustering(clusters);

  std::cout << "\nHierarchical cluster:\n";
  big_cluster -> show();

  // Deleteting big_cluster call all children destructors in a recursive fashion
  delete big_cluster;

}

// This example was contrasted whit the algoritm of web site https://people.revoledu.com/kardi/tutorial/Clustering/Online-Hierarchical-Clustering.html
void check_example() {

  const int dim = 2;

  std::set < Cluster * > clusters;

  // Add new data
  {

    Cluster * new_cluster = new Cluster();

    // Configure new cluster
    new_cluster -> set_label("A");
    float point[dim] = {
      1,
      1
    };
    new_cluster -> set_centroid(point, dim);

    clusters.insert(new_cluster);

  }

  // Add new data
  {

    Cluster * new_cluster = new Cluster();

    // Configure new cluster
    new_cluster -> set_label("B");
    float point[dim] = {
      1.5,
      1.5
    };
    new_cluster -> set_centroid(point, dim);

    clusters.insert(new_cluster);

  }

  // Add new data
  {
    Cluster * new_cluster = new Cluster();

    // Configure new cluster
    new_cluster -> set_label("C");
    float point[dim] = {
      2,
      2
    };
    new_cluster -> set_centroid(point, dim);

    clusters.insert(new_cluster);

  }

  // Add new data
  {

    Cluster * new_cluster = new Cluster();

    // Configure new cluster
    new_cluster -> set_label("D");
    float point[dim] = {
      1,
      1
    };
    new_cluster -> set_centroid(point, dim);

    clusters.insert(new_cluster);

  }

  // Add new data
  {

    Cluster * new_cluster = new Cluster();

    // Configure new cluster
    new_cluster -> set_label("E");
    float point[dim] = {
      9,
      9
    };
    new_cluster -> set_centroid(point, dim);

    clusters.insert(new_cluster);

  }

  // Add new data
  {

    Cluster * new_cluster = new Cluster();

    // Configure new cluster
    new_cluster -> set_label("F");
    float point[dim] = {
      1,
      1
    };
    new_cluster -> set_centroid(point, dim);

    clusters.insert(new_cluster);

  }

  // Add new data
  {

    Cluster * new_cluster = new Cluster();

    // Configure new cluster
    new_cluster -> set_label("G");
    float point[dim] = {
      7,
      7
    };
    new_cluster -> set_centroid(point, dim);

    clusters.insert(new_cluster);

  }

  // Add new data
  {

    Cluster * new_cluster = new Cluster();

    // Configure new cluster
    new_cluster -> set_label("H");
    float point[dim] = {
      9,
      11
    };
    new_cluster -> set_centroid(point, dim);

    clusters.insert(new_cluster);

  }

  std::cout << "Data:\n";

  std::set < Cluster * > ::iterator it = clusters.begin();
  for (int i = 0; i < clusters.size(); i++) {
    ( * it) -> show();
    std::advance(it, 1);
  }

  Cluster * big_cluster = clustering(clusters);

  std::cout << "\nHierarchical cluster:\n";
  big_cluster -> show();
  
  // Deleteting big_cluster call all children destructors in a recursive fashion
  delete big_cluster;

}

int main() {

  //check_random_example();
  check_example();

  return 0;

}
