#include "matrix.h"

Matrix::Matrix(unsigned m_rowsize, unsigned m_colsize, double initial)
{
    rowsize = m_rowsize;
    colsize = m_colsize;
    matrix.resize(rowsize);
    for (unsigned i = 0; i < matrix.size(); i++)
    {
        matrix[i].resize(colsize, initial);
    }
}


Matrix::Matrix(const Matrix & B)
{
    colsize = B.getCols();
    rowsize = B.getRows();
    matrix.resize(rowsize);
    for (unsigned i=0; i<matrix.size(); i++)
    {
        matrix[i].resize(colsize, 0);
    }
    for (unsigned i=0; i<rowsize;i++)
    {
        for (unsigned j=0; j<colsize; j++)
        {
            (*this)(i,j) = B.matrix[i][j];
        }
    }
}

Matrix::Matrix(int shape)
{
    rowsize = shape;
    colsize = shape;
    matrix.resize(rowsize);
    for (unsigned i = 0; i < matrix.size(); i++)
    {
        matrix[i].resize(shape, 0);
        (*this)(i,i) = 1;
    }
}


// Addition of Two Matrices
Matrix Matrix::operator+(Matrix &B)
{
    Matrix sum(colsize, rowsize, 0.0);
    unsigned i,j;
    for (i = 0; i < rowsize; i++)
    {
        for (j = 0; j < colsize; j++)
        {
            sum(i,j) = matrix[i][j] + B(i,j);
        }
    }
    return sum;
}



// Subtraction of Two Matrices
Matrix Matrix::operator-(Matrix & B){
    Matrix diff(colsize, rowsize, 0.0);
    unsigned i,j;
    for (i = 0; i < rowsize; i++)
    {
        for (j = 0; j < colsize; j++)
        {
            diff(i,j) = matrix[i][j] - B(i,j);
        }
    }
    return diff;
}


// Multiplication of Two Matrices
Matrix Matrix::operator*(Matrix & B)
{
    Matrix multip(rowsize,B.getCols(),0.0);
    if(colsize == B.getRows())
    {
        unsigned i,j,k;
        double temp = 0.0;
        for (i = 0; i < rowsize; i++)
        {
            for (j = 0; j < B.getCols(); j++)
            {
                temp = 0.0;
                for (k = 0; k < colsize; k++)
                {
                    temp += matrix[i][k] * B(k,j);
                }
                multip(i,j) = temp;
            }
        }
        return multip;
    }
    else
    {
        throw std::invalid_argument( "received incompatible matrices" );
    }
}


Matrix& Matrix::operator+=(Matrix &B)
{
    *this = (*this) + B;
    return *this;
}



Matrix& Matrix::operator*=(Matrix &B)
{
    *this = (*this) * B;
    return *this;
}

// Scalar Addition
Matrix Matrix::operator+(long double scalar){
    Matrix result(rowsize,colsize,0.0);
    unsigned i,j;
    for (i = 0; i < rowsize; i++)
    {
        for (j = 0; j < colsize; j++)
        {
            result(i,j) = matrix[i][j] + scalar;
        }
    }
    return result;
}

// Scalar Subraction
Matrix Matrix::operator-(long double scalar){
    Matrix result(rowsize,colsize,0.0);
    unsigned i,j;
    for (i = 0; i < rowsize; i++)
    {
        for (j = 0; j < colsize; j++)
        {
            result(i,j) = this->matrix[i][j] - scalar;
        }
    }
    return result;
}

// Scalar Multiplication
Matrix Matrix::operator*(long double scalar){
    Matrix result(rowsize,colsize,0.0);
    unsigned i,j;
    for (i = 0; i < rowsize; i++)
    {
        for (j = 0; j < colsize; j++)
        {
            result(i,j) = matrix[i][j] * scalar;
        }
    }
    return result;
}

// Scalar Division
Matrix Matrix::operator/(long double scalar){
    Matrix result(rowsize,colsize,0.0);
    unsigned i,j;
    for (i = 0; i < rowsize; i++)
    {
        for (j = 0; j < colsize; j++)
        {
            result(i,j) = this->matrix[i][j] / scalar;
        }
    }
    return result;
}

// Returns value of given location when asked in the form A(x,y)
long double& Matrix::operator()(const unsigned &rowNo, const unsigned & colNo)
{
    if ((0 > rowNo) + (rowNo > this->getRows()) + (colNo > this->getCols()) +  (0 > colNo))
    {        
        throw std::invalid_argument( "Invalid column / row number");
    }
    return this->matrix[rowNo][colNo];
}


// No brainer - returns row #
unsigned Matrix::getRows() const
{
    return rowsize;
}

// returns col #
unsigned Matrix::getCols() const
{
    return colsize;
}

// Take any given matrices transpose and returns another matrix
Matrix Matrix::transpose()
{
    Matrix Transpose(colsize, rowsize, 0.0);
    for (unsigned i = 0; i < colsize; i++)
    {
        for (unsigned j = 0; j < rowsize; j++) {
            Transpose(i,j) = matrix[j][i];
        }
    }
    return Transpose;
}

// Prints the matrix beautifully

std::ostream & operator<<(std::ostream & os, Matrix & m)
{
    os << std::endl;
    for (unsigned i = 0; i < m.getRows(); i++) {
        for (unsigned j = 0; j < m.getCols(); j++) {
            double el = m(i,j);
            os << " " << el << " ";
        }
        os << std::endl;
    }
    return os;
}

// implementation GNU library
double Matrix::det_3(unsigned r, unsigned c)
{
    double det;
    double a[3][3];

    for (int i=0; i<4; i++)
    {
        for (int j=0; j<4; j++)
        {
            if ((i<r)*(j<c))
            {
                a[i][j] = (*this)(i,j);
            }
            if ((i<r)*(j>c))
            {
                a[i][j-1] = (*this)(i,j);
            }
            if ((i>r)*(j<c))
            {
                a[i-1][j] = (*this)(i,j);
            }
            if ((i>r)*(j>c))
            {
                a[i-1][j-1] = (*this)(i,j);
            }        
        }
    }
    for (int i=0; i<3;i++)
    {
            det = det + (a[0][i] * (a[1][(i+1)%3] * a[2][(i+2)%3] - a[1][(i+2)%3] * a[2][(i+1)%3]));
    }
    return det;

}

Matrix Matrix::inv_4()
{
    long double m[16];
    for (int i=0; i < 4; i++)
    {
        m[4*i] = (*this)(i,0);
        m[1+i*4] = (*this)(i,1);
        m[2+i*4] = (*this)(i,2);
        m[3+i*4] = (*this)(i,3);
    }

    long double invOut[16];
    long double inv[16], det;

    inv[0] = m[5]  * m[10] * m[15] - 
             m[5]  * m[11] * m[14] - 
             m[9]  * m[6]  * m[15] + 
             m[9]  * m[7]  * m[14] +
             m[13] * m[6]  * m[11] - 
             m[13] * m[7]  * m[10];

    inv[4] = -m[4]  * m[10] * m[15] + 
              m[4]  * m[11] * m[14] + 
              m[8]  * m[6]  * m[15] - 
              m[8]  * m[7]  * m[14] - 
              m[12] * m[6]  * m[11] + 
              m[12] * m[7]  * m[10];

    inv[8] = m[4]  * m[9] * m[15] - 
             m[4]  * m[11] * m[13] - 
             m[8]  * m[5] * m[15] + 
             m[8]  * m[7] * m[13] + 
             m[12] * m[5] * m[11] - 
             m[12] * m[7] * m[9];

    inv[12] = -m[4]  * m[9] * m[14] + 
               m[4]  * m[10] * m[13] +
               m[8]  * m[5] * m[14] - 
               m[8]  * m[6] * m[13] - 
               m[12] * m[5] * m[10] + 
               m[12] * m[6] * m[9];

    inv[1] = -m[1]  * m[10] * m[15] + 
              m[1]  * m[11] * m[14] + 
              m[9]  * m[2] * m[15] - 
              m[9]  * m[3] * m[14] - 
              m[13] * m[2] * m[11] + 
              m[13] * m[3] * m[10];

    inv[5] = m[0]  * m[10] * m[15] - 
             m[0]  * m[11] * m[14] - 
             m[8]  * m[2] * m[15] + 
             m[8]  * m[3] * m[14] + 
             m[12] * m[2] * m[11] - 
             m[12] * m[3] * m[10];

    inv[9] = -m[0]  * m[9] * m[15] + 
              m[0]  * m[11] * m[13] + 
              m[8]  * m[1] * m[15] - 
              m[8]  * m[3] * m[13] - 
              m[12] * m[1] * m[11] + 
              m[12] * m[3] * m[9];

    inv[13] = m[0]  * m[9] * m[14] - 
              m[0]  * m[10] * m[13] - 
              m[8]  * m[1] * m[14] + 
              m[8]  * m[2] * m[13] + 
              m[12] * m[1] * m[10] - 
              m[12] * m[2] * m[9];

    inv[2] = m[1]  * m[6] * m[15] - 
             m[1]  * m[7] * m[14] - 
             m[5]  * m[2] * m[15] + 
             m[5]  * m[3] * m[14] + 
             m[13] * m[2] * m[7] - 
             m[13] * m[3] * m[6];

    inv[6] = -m[0]  * m[6] * m[15] + 
              m[0]  * m[7] * m[14] + 
              m[4]  * m[2] * m[15] - 
              m[4]  * m[3] * m[14] - 
              m[12] * m[2] * m[7] + 
              m[12] * m[3] * m[6];

    inv[10] = m[0]  * m[5] * m[15] - 
              m[0]  * m[7] * m[13] - 
              m[4]  * m[1] * m[15] + 
              m[4]  * m[3] * m[13] + 
              m[12] * m[1] * m[7] - 
              m[12] * m[3] * m[5];

    inv[14] = -m[0]  * m[5] * m[14] + 
               m[0]  * m[6] * m[13] + 
               m[4]  * m[1] * m[14] - 
               m[4]  * m[2] * m[13] - 
               m[12] * m[1] * m[6] + 
               m[12] * m[2] * m[5];

    inv[3] = -m[1] * m[6] * m[11] + 
              m[1] * m[7] * m[10] + 
              m[5] * m[2] * m[11] - 
              m[5] * m[3] * m[10] - 
              m[9] * m[2] * m[7] + 
              m[9] * m[3] * m[6];

    inv[7] = m[0] * m[6] * m[11] - 
             m[0] * m[7] * m[10] - 
             m[4] * m[2] * m[11] + 
             m[4] * m[3] * m[10] + 
             m[8] * m[2] * m[7] - 
             m[8] * m[3] * m[6];

    inv[11] = -m[0] * m[5] * m[11] + 
               m[0] * m[7] * m[9] + 
               m[4] * m[1] * m[11] - 
               m[4] * m[3] * m[9] - 
               m[8] * m[1] * m[7] + 
               m[8] * m[3] * m[5];

    inv[15] = m[0] * m[5] * m[10] - 
              m[0] * m[6] * m[9] - 
              m[4] * m[1] * m[10] + 
              m[4] * m[2] * m[9] + 
              m[8] * m[1] * m[6] - 
              m[8] * m[2] * m[5];

    det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

    if (det == 0)
        return false;

    det = 1.0 / det;

    for (int i = 0; i < 16; i++)
        invOut[i] = inv[i] * det;

    Matrix inverted_matrix = Matrix(4,4,0);
    for (int i=0; i < 4; i++)
    {
        inverted_matrix(i,0) = invOut[4*i];
        inverted_matrix(i,1) = invOut[4*i+1];
        inverted_matrix(i,2) = invOut[4*i+2];
        inverted_matrix(i,3) = invOut[4*i+3];
    }

    return inverted_matrix;
}




Matrix Matrix::pinv_4()
{
    Matrix m_transposed = this->transpose();

    Matrix multiple = (*this) * m_transposed;
    Matrix multiple_inv = multiple.inv_4();
    Matrix pinv = m_transposed * multiple_inv;
    return pinv;
}


Matrix Matrix::slice(std::tuple<int, int> r_start_end, std::tuple<int, int> c_start_end)
{
    int rows = std::get<1>(r_start_end) - std::get<0>(r_start_end);
    int cols = std::get<1>(c_start_end) - std::get<0>(c_start_end);
    Matrix m_sliced = Matrix(rows, cols, 0);

    for (int row = 0; row < rows; row ++)
    {
        for (int col = 0; col < cols; col ++)
        {
            int row_o = row + std::get<0>(r_start_end);
            int col_o = col + std::get<0>(c_start_end);
            m_sliced(row,col) = (*this)(row_o, col_o);
        }
    }
    return m_sliced;
}


Matrix Matrix::append(Matrix & m_2, int axis)
{
    int rows1 = this->getRows();
    int rows2 = m_2.getRows();
    int cols1 = this->getCols();
    int cols2 = m_2.getCols();

    if ( ((cols1 != cols2) & (axis == 0)) + ((rows1 != rows2) & (axis == 1)) )
    {
        throw std::invalid_argument( "These matrices don't have an equal amount of rows");   
    }
    if ((axis != 0) * (axis != 1))
    {
        throw std::invalid_argument( "Invalid Axis");   
    }

    if (axis == 1)
    {
        Matrix new_matrix = Matrix(rows1, (cols1+cols2), 0);
        for(int i=0; i<rows1; i++)
        {
            for (int j=0; j<(cols1+cols2);j++)
            {
                if (j < cols1)
                {

                    new_matrix(i,j) = (*this)(i,j);
                }
                else 
                {
                    new_matrix(i,j) = m_2(i,j-cols1);
                }
            }
        }
        return new_matrix;
    }
    else if (axis == 0)
    {
        Matrix new_matrix = Matrix((rows1+rows2), cols1, 0);
        for(int i=0; i<(rows1+rows2); i++)
        {
            for (int j=0; j<cols1;j++)
            {
                if (j < rows1)
                {

                    new_matrix(i,j) = (*this)(i,j);
                }
                else 
                {
                    new_matrix(i,j) = m_2(i-rows2,j);
                }
            }
        }
        return new_matrix;
    }
    throw std::invalid_argument("Append function doesn't work");   
}

double Matrix::norm()
{
    double norm = 0;
    for (int i=0; i < (*this).getRows(); i++)
    {
        for (int j=0; j < (*this).getCols(); j++)
        {
            norm += ((*this)(i,j)*(*this)(i,j));
        }
    }
    return norm;
}