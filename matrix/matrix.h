#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <tuple>
#include <stdexcept>
#include <cmath>

using std::vector;
using std::tuple;
using std::ostream;

class Matrix {
    private:
    unsigned int rowsize;
    unsigned int colsize;
    vector<vector<long double>> matrix;

    public:
    // Constructor/Destructor
    Matrix(){;};
    Matrix(unsigned, unsigned, double);
    Matrix(const Matrix &);
    Matrix(int);

    ~Matrix(){};

    // Create matrix types

    // Matrix operations
    Matrix operator+(Matrix &);
    Matrix operator-(Matrix &);
    Matrix operator*(Matrix &);
    Matrix& operator*=(Matrix &);
    Matrix& operator+=(Matrix &);

    Matrix transpose();

    // Appends column matrices with same amount of rows and returns the lvalue for a new matrix.

    // Scalar operations
    Matrix operator+(long double);
    Matrix operator-(long double);
    Matrix operator*(long double);
    Matrix operator/(long double);

    // Aesthetic Methods
    long double& operator()(const unsigned &, const unsigned &);
    Matrix slice(std::tuple<int, int>, std::tuple<int, int>);
    Matrix append(Matrix &, int);
    double norm();


    unsigned getRows() const;
    unsigned getCols() const;
    friend std::ostream& operator<<(ostream &, Matrix &);

    // Inverter
    Matrix inv_4();
    Matrix inv_iter();
    Matrix pinv_4();
    Matrix pinv_iter();

    double det_3(unsigned, unsigned);
};

#endif