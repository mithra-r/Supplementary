// Simple linear algebra class that supports:
// dotVV: dot product between vectors
// dot: dot product between vectors or matrices
// normL2: L2 norm
// T: transpose
// addition/subtraction with other vector
// element wise multiply/divide with scalar
//
// Part of supplementary material of "Overhang control in topology optimization: 
// a comparison of continuous front propagation-based and discrete layer-by-layer
// overhang control", E. van de Ven, R. Maas, C. Ayas, M. Langelaar, F. van Keulen,
// 2019
//
// Code by Emiel van de Ven, 2019
// emiel@emielvandeven.nl

#include <vector>
#include <math.h>

class MatrixXX {
    public:
        std::vector<double> dat;
        int rows, cols;
        MatrixXX(int _rows = 0, int _cols = 0, bool eye = false): rows(_rows), cols(_cols) {
            dat.resize(rows*cols,0.0);
            if (eye) for (int i=0; i<dat.size(); i+=cols+1) dat[i] = 1.0;
        }
        MatrixXX(const MatrixXX &mat): rows(mat.rows), cols(mat.cols), dat(mat.dat) {}
        MatrixXX(int _rows, int _cols, std::initializer_list<double> l): rows(_rows), cols(_cols), dat(l) {}
        inline double dotVV(const MatrixXX &mat2) {
            double res = 0.0;
            for (int i=0; i<dat.size(); i++) res += dat[i]*mat2.dat[i];
            return res;
        }
        inline double normL2() {
            return sqrt(dotVV(*this));
        }
        inline MatrixXX dot(const MatrixXX &mat2) {
            MatrixXX mat3(rows, mat2.cols);
            for (int i=0; i<rows; i++) {
                for (int j=0; j<mat2.cols; j++) {
                    for (int k=0; k<cols; k++) {
                        mat3.dat[i*mat3.cols+j] += dat[i*cols+k]*mat2.dat[k*mat2.cols+j];
                    }
                }
            }
            return mat3;
        }
        inline MatrixXX T() {
            MatrixXX mat(cols, rows);
            for (int i=0; i<rows; i++) {
                for (int j=0; j<cols; j++) {
                    mat.dat[j*mat.cols + i] = dat[i*cols+j];
                }
            }
            return mat;
        }
        inline MatrixXX & operator += (const MatrixXX &mat2) {
            for (int i=0; i<dat.size(); i++) dat[i] += mat2.dat[i]; return *this; }
        inline MatrixXX & operator -= (const MatrixXX &mat2) { 
            for (int i=0; i<dat.size(); i++) dat[i] -= mat2.dat[i]; return *this; }
        inline MatrixXX & operator *= (double val) { 
            for (int i=0; i<dat.size(); i++) dat[i] *= val; return *this; }
        inline MatrixXX & operator /= (double val) { 
            for (int i=0; i<dat.size(); i++) dat[i] /= val; return *this; }
        MatrixXX operator + (const MatrixXX &mat2) {
            return MatrixXX(*this) += mat2; }
        MatrixXX operator - (const MatrixXX &mat2) {
            return MatrixXX(*this) -= mat2; }
        MatrixXX operator * (double val) {
            return MatrixXX(*this) *= val; }
        MatrixXX operator / (double val) {
            return MatrixXX(*this) /= val; }
};