// Copyright 2022 Democritus University of Thrace
// Integrated Circuits Lab
// 
// Example using Fast-Float4HLS to build a Graph Convolutional Network
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <fstream>
#include <istream>
#include <sstream>
#include <string>

#include <vector>
#include <cstdlib>
#include <vector>

#include "defs.h"
#include "mc_scverify.h"

// print a 1d matrix
template <typename Type, int N>
void print1d (Type a[N]) {
  for (int i=0; i < N; i++) {
    std::cout << a[i] << "\t";
  }
  std::cout << std::endl;
}

// print a 1d matrix of fast-float datatype
template <int N>
void print1d (btype a[N]) {
  for (int i=0; i < N; i++) {
    std::cout << a[i].to_float() << "\t";
  }
  std::cout << std::endl;
}

// print a 2d matrix
template <typename Type, int N, int L>
void print2d (Type a[N][L]) {
  for (int i=0; i < N; i++) {
    for (int j=0; j < L; j++) {
      std::cout << a[i][j] << "\t";
    }
    std::cout << std::endl;
  }
}

// print a 2d matrix of fast-float datatype
template <int N, int L>
void print2d (btype a[N][L]) {
  for (int i=0; i < N; i++) {
    for (int j=0; j < L; j++) {
      std::cout << a[i][j].to_float() << "\t";
    }
    std::cout << std::endl;
  }
}


//read a 2d matrix from a txt file
template<typename Type, int N, int L>
void read_data(std::vector< std::vector<Type> > &img, std::string filename) {
  std::ifstream myFile(filename);
  if(myFile.is_open()) { //throw std::runtime_error("Could not open file");  // Make sure the file is open
  
  std::vector<Type> row;
  std::string line;
  Type val;
  int rowIdx = 0;

  // Read data, line by line
  while(std::getline(myFile, line)) {   
    std::stringstream ss(line);  // Create a stringstream of the current line
    int colIdx = 0;  // Keep track of the current column index
    
    row.clear();    

    // Extract each integer
    while(ss >> val){
        row.push_back(val);  // Write current input value
        if(ss.peek() == ',') ss.ignore(); // If the next token is a comma, ignore it and move on
               
        colIdx++;  // Increment the Column index
        if (colIdx == L) break;
    }
    img.push_back(row);
    rowIdx++;  // Increment the Row index
    if (rowIdx == N) break;
  }
    myFile.close();  // Close file
  }
}


//read a matrix in csr format from a txt file
//each one of three sub-matrices in csr format assumed to be in a different line
template<typename Type, int N, int nZ>
void read_adj(std::vector<int> &a_row, std::vector<int> &a_col, std::vector<Type> &a_val, std::string filename) {
  std::ifstream myFile(filename);
  if(myFile.is_open()) { //throw std::runtime_error("Could not open file");  // Make sure the file is open

    std::string line;
    int int_val;
    Type def_val;
    int rowIdx = 0;

    // Read the first line -- row_ind for CSR
    std::getline(myFile, line);   
    std::stringstream ss1(line);  // Create a stringstream of the current line
    int colIdx = 0;  // Keep track of the current column index
    
    // Extract each integer
    while(ss1 >> int_val){        
      a_row.push_back(int_val);
      if(ss1.peek() == ',') ss1.ignore(); // If the next token is a comma, ignore it and move on
      colIdx++;  // Increment the Column index
    }
    ss1.clear();
    
    // Read the second line -- col_ind for CSR
    std::getline(myFile, line);   
    std::stringstream ss2(line);  // Create a stringstream of the current line
    colIdx = 0;  // Keep track of the current column index
    
    // Extract each integer
    while(ss2 >> int_val){        
      a_col.push_back(int_val);
      if(ss2.peek() == ',') ss2.ignore(); // If the next token is a comma, ignore it and move on
      colIdx++;  // Increment the Column index
    }
    ss2.clear();

    // Read the third line -- value for CSR
    std::getline(myFile, line);   
    std::stringstream ss3(line);  // Create a stringstream of the current line
    colIdx = 0;  // Keep track of the current column index
    
    // Extract each integer
    while(ss3 >> def_val){        
      a_val.push_back(def_val);
      if(ss3.peek() == ',') ss3.ignore(); // If the next token is a comma, ignore it and move on
      colIdx++;  // Increment the Column index
    }
    ss3.clear();
    
    myFile.close();  // Close file
  }
}

