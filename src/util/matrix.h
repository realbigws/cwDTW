#ifndef MATRIX_H__
#define MATRIX_H__

#include <iostream>
#include <cmath>
#include <cstring>
#include <cassert>
#include "exception.h"

namespace mx{

#define _FMX(X,args...)		(float[]){X, ##args}
#define _DMX(X,args...)		(double[]){X, ##args}
#define _IMX(X,args...)		(int[]){X, ##args}  
    
template<int D1, int D2, typename T=float>
class Matrix{
private:
    union{
	int cols;
	int width;
    }; 
    union{
	int rows;
	int height;
    }; 
private:
    T data[D1][D2];
public:
    Matrix(){height = D1; width = D2;}
    Matrix(const Matrix<D1, D2, T>& m){height = D1; width = D2;memcpy(data, m.data, sizeof(T)*rows*cols);}
    Matrix(const T* pm){height = D1; width = D2;memcpy(data, pm, sizeof(T)*rows*cols);}
public:
    inline int Width() const{return width;}
    inline int Cols() const{return cols;}
    inline int Height() const{return height;}
    inline int Rows() const{return rows;}
    inline T& V(int x, int y){return data[x][y];}
    inline T V(int x, int y) const{return data[x][y];}
    //if is column matrix. also be seen as [a,...,d]
    inline T& V(int x){assert(cols == 1); return data[x][0];}
    inline T V(int x) const{assert(cols == 1); return data[x][0];}
    //access the data as array
    inline T& D(int x){assert(x>=0 && x<rows*cols); return((T*)data)[x];}//return data[x/cols][x%cols];}
    inline T D(int x) const{assert(x>=0 && x<rows*cols); return((T*)data)[x];}//return data[x/cols][x%cols];}
    inline const T* D() const{return (T*)data;}
    inline T* D() {return (T*)data;}
    
    bool InvertGaussJordan()
    {
	assert(D1 ==  D2);
	//allocate the temp memory
	int* pnRow  =  new int[D1];
	int* pnCol  =  new int[D1];
	// elimination
	double d = 0, p = 0;
	for(int k = 0; k <= D1-1; k++){
	    d = 0.0;
	    for(int i = k; i <= D1-1; i++){
		for(int j = k; j <= D1-1; j++){
		    p = fabs(data[i][j]);
		    if(p>d){
			d = p;
			pnRow[k] = i;
			pnCol[k] = j;
		    }
		}
	    }
	    // failed
	    if(d == 0.0){
		delete [] pnRow; delete [] pnCol;
		return false;
	    }

	    if(pnRow[k] != k){
		for(int j = 0; j <= D1-1; j++){
		    p = data[k][j];
		    data[k][j] = data[pnRow[k]][j];
		    data[pnRow[k]][j] = p;
		}
	    }

	    if(pnCol[k] !=  k){
		for(int i = 0; i <= D1-1; i++){
		    p = data[i][k];
		    data[i][k] = data[i][pnCol[k]];
		    data[i][pnCol[k]] = p;
		}
	    }

	    data[k][k] = 1/data[k][k];
	    for(int j = 0; j <= D1-1; j++){
		if(j != k){
		    data[k][j] *= data[k][k];
		}
	    }

	    for(int i = 0; i <= D1-1; i++){
		if(i != k){
		    for(int j = 0; j <= D1-1; j++){
			if(j != k){
			    data[i][j] -= data[i][k]*data[k][j];
			}
		    }
		}
	    }

	    for(int i = 0; i <= D1-1; i++){
		if(i != k){
		    data[i][k] *= -data[k][k];
		}
	    }
	}
	//swap rows and columns
	for(int k = D1-1; k>= 0; k--){
	    if(pnCol[k] != k){
		for(int j = 0; j <= D1-1; j++){
		    p = data[k][j];
		    data[k][j] = data[pnCol[k]][j];
		    data[pnCol[k]][j] = p;
		}
	    }

	    if(pnRow[k] != k){
		for(int i = 0; i <= D1-1; i++){
		    p = data[i][k];
		    data[i][k] = data[i][pnRow[k]];
		    data[i][pnRow[k]] = p;
		}
	    }
	}

	delete[] pnRow;
	delete[] pnCol;
	
	return true;
    }
    
    void Print(std::ostream& s) const
    {
	s<<rows<<"x"<<cols<<"\n";
	for(int x = 0; x < rows; x++){
	    for(int y = 0; y < cols; y++){
		s<<data[x][y]<<" ";
	    }
	    s<<"\n";
	}
    }
    
    const Matrix<D1,D2,T>& operator = (const Matrix<D1,D2,T>& m)
    {
	assert(rows == m.rows && cols == m.cols);
        // copy the pointer
        memcpy(data, m.data, sizeof(T)*rows*cols);
        
        return *this;
    }
    
    const Matrix<D1,D2,T>& operator = (const T* pm)
    {
        memcpy(data, pm, sizeof(T)*rows*cols);
        return *this;
    }
    
    const Matrix<D1,D2,T>& operator *= (T v)
    {
        for(int i = 0; i < rows; i++){
            for (int j = 0; j < cols; j++){
		data[i][j] *= v;
	    }
        }

        return *this;
    }
    
public:
    template<typename K> 
    friend void Product333(const Matrix<3,3,K>& src1, const Matrix<3,3,K>& src2, Matrix<3,3,K>* destin);
    template <typename K>
    friend void Product444(const Matrix<4,4,K>& src1, const Matrix<4,4,K>& src2, Matrix<4,4,K>* destin);
    template <typename K>
    friend void Product121(const Matrix<1,2,K>& src1, const Matrix<2,1,K>& src2, Matrix<1,1,K>* destin);
    template <typename K>
    friend void Product131(const Matrix<1,3,K>& src1, const Matrix<3,1,K>& src2, Matrix<1,1,K>* destin);
    template <typename K>
    friend void Product331(const Matrix<3,3,K>& src1, const Matrix<3,1,K>& src2, Matrix<3,1,K>* destin);
    template <typename K>
    friend void Product341(const Matrix<3,4,K>& src1, const Matrix<4,1,K>& src2, Matrix<3,1,K>* destin);
    template <typename K>
    friend void Product441(const Matrix<4,4,K>& src1, const Matrix<4,1,K>& src2, Matrix<4,1,K>* destin);

};

template <typename T>
void Product333(const Matrix<3,3,T>& src1, const Matrix<3,3,T>& src2, Matrix<3,3,T>* destin)
{
    destin->data[0][0] = src1.data[0][0]*src2.data[0][0]+src1.data[0][1]*src2.data[1][0]+src1.data[0][2]*src2.data[2][0];
    destin->data[0][1] = src1.data[0][0]*src2.data[0][1]+src1.data[0][1]*src2.data[1][1]+src1.data[0][2]*src2.data[2][1];
    destin->data[0][2] = src1.data[0][0]*src2.data[0][2]+src1.data[0][1]*src2.data[1][2]+src1.data[0][2]*src2.data[2][2];
    
    destin->data[1][0] = src1.data[1][0]*src2.data[0][0]+src1.data[1][1]*src2.data[1][0]+src1.data[1][2]*src2.data[2][0];
    destin->data[1][1] = src1.data[1][0]*src2.data[0][1]+src1.data[1][1]*src2.data[1][1]+src1.data[1][2]*src2.data[2][1];
    destin->data[1][2] = src1.data[1][0]*src2.data[0][2]+src1.data[1][1]*src2.data[1][2]+src1.data[1][2]*src2.data[2][2];

    destin->data[2][0] = src1.data[2][0]*src2.data[0][0]+src1.data[2][1]*src2.data[1][0]+src1.data[2][2]*src2.data[2][0];
    destin->data[2][1] = src1.data[2][0]*src2.data[0][1]+src1.data[2][1]*src2.data[1][1]+src1.data[2][2]*src2.data[2][1];
    destin->data[2][2] = src1.data[2][0]*src2.data[0][2]+src1.data[2][1]*src2.data[1][2]+src1.data[2][2]*src2.data[2][2];
}

template <typename T>
void Product444(const Matrix<4,4,T>& src1, const Matrix<4,4,T>& src2, Matrix<4,4,T>* destin)
{
    destin->data[0][0] = src1.data[0][0]*src2.data[0][0]+src1.data[0][1]*src2.data[1][0]+src1.data[0][2]*src2.data[2][0]+src1.data[0][3]*src2.data[3][0];
    destin->data[0][1] = src1.data[0][0]*src2.data[0][1]+src1.data[0][1]*src2.data[1][1]+src1.data[0][2]*src2.data[2][1]+src1.data[0][3]*src2.data[3][1];
    destin->data[0][2] = src1.data[0][0]*src2.data[0][2]+src1.data[0][1]*src2.data[1][2]+src1.data[0][2]*src2.data[2][2]+src1.data[0][3]*src2.data[3][2];
    destin->data[0][3] = src1.data[0][0]*src2.data[0][3]+src1.data[0][1]*src2.data[1][3]+src1.data[0][2]*src2.data[2][3]+src1.data[0][3]*src2.data[3][3];
    
    destin->data[1][0] = src1.data[1][0]*src2.data[0][0]+src1.data[1][1]*src2.data[1][0]+src1.data[1][2]*src2.data[2][0]+src1.data[1][3]*src2.data[3][0];
    destin->data[1][1] = src1.data[1][0]*src2.data[0][1]+src1.data[1][1]*src2.data[1][1]+src1.data[1][2]*src2.data[2][1]+src1.data[1][3]*src2.data[3][1];
    destin->data[1][2] = src1.data[1][0]*src2.data[0][2]+src1.data[1][1]*src2.data[1][2]+src1.data[1][2]*src2.data[2][2]+src1.data[1][3]*src2.data[3][2];
    destin->data[1][3] = src1.data[1][0]*src2.data[0][3]+src1.data[1][1]*src2.data[1][3]+src1.data[1][2]*src2.data[2][3]+src1.data[1][3]*src2.data[3][3];

    destin->data[2][0] = src1.data[2][0]*src2.data[0][0]+src1.data[2][1]*src2.data[1][0]+src1.data[2][2]*src2.data[2][0]+src1.data[2][3]*src2.data[3][0];
    destin->data[2][1] = src1.data[2][0]*src2.data[0][1]+src1.data[2][1]*src2.data[1][1]+src1.data[2][2]*src2.data[2][1]+src1.data[2][3]*src2.data[3][1];
    destin->data[2][2] = src1.data[2][0]*src2.data[0][2]+src1.data[2][1]*src2.data[1][2]+src1.data[2][2]*src2.data[2][2]+src1.data[2][3]*src2.data[3][2];
    destin->data[2][3] = src1.data[2][0]*src2.data[0][3]+src1.data[2][1]*src2.data[1][3]+src1.data[2][2]*src2.data[2][3]+src1.data[2][3]*src2.data[3][3];

    destin->data[3][0] = src1.data[3][0]*src2.data[0][0]+src1.data[3][1]*src2.data[1][0]+src1.data[3][2]*src2.data[2][0]+src1.data[3][3]*src2.data[3][0];
    destin->data[3][1] = src1.data[3][0]*src2.data[0][1]+src1.data[3][1]*src2.data[1][1]+src1.data[3][2]*src2.data[2][1]+src1.data[3][3]*src2.data[3][1];
    destin->data[3][2] = src1.data[3][0]*src2.data[0][2]+src1.data[3][1]*src2.data[1][2]+src1.data[3][2]*src2.data[2][2]+src1.data[3][3]*src2.data[3][2];
    destin->data[3][3] = src1.data[3][0]*src2.data[0][3]+src1.data[3][1]*src2.data[1][3]+src1.data[3][2]*src2.data[2][3]+src1.data[3][3]*src2.data[3][3];
}

template <typename T>
void Product121(const Matrix<1,2,T>& src1, const Matrix<2,1,T>& src2, Matrix<1,1,T>* destin)
{
    destin->data[0][0] = src1.data[0][0]*src2.data[0][0]+src1.data[0][1]*src2.data[1][0];
}

template <typename T>
void Product131(const Matrix<1,3,T>& src1, const Matrix<3,1,T>& src2, Matrix<1,1,T>* destin)
{
    destin->data[0][0] = src1.data[0][0]*src2.data[0][0]+src1.data[0][1]*src2.data[1][0]+src1.data[0][2]*src2.data[2][0];
}

template <typename T>
void Product331(const Matrix<3,3,T>& src1, const Matrix<3,1,T>& src2, Matrix<3,1,T>* destin)
{
    destin->data[0][0] = src1.data[0][0]*src2.data[0][0]+src1.data[0][1]*src2.data[1][0]+src1.data[0][2]*src2.data[2][0];
    destin->data[1][0] = src1.data[1][0]*src2.data[0][0]+src1.data[1][1]*src2.data[1][0]+src1.data[1][2]*src2.data[2][0];
    destin->data[2][0] = src1.data[2][0]*src2.data[0][0]+src1.data[2][1]*src2.data[1][0]+src1.data[2][2]*src2.data[2][0];
}

template <typename T>
void Product341(const Matrix<3,4,T>& src1, const Matrix<4,1,T>& src2, Matrix<3,1,T>* destin)
{
    destin->data[0][0] = src1.data[0][0]*src2.data[0][0]+src1.data[0][1]*src2.data[1][0]+src1.data[0][2]*src2.data[2][0]+src1.data[0][3]*src2.data[3][0];
    destin->data[1][0] = src1.data[1][0]*src2.data[0][0]+src1.data[1][1]*src2.data[1][0]+src1.data[1][2]*src2.data[2][0]+src1.data[1][3]*src2.data[3][0];
    destin->data[2][0] = src1.data[2][0]*src2.data[0][0]+src1.data[2][1]*src2.data[1][0]+src1.data[2][2]*src2.data[2][0]+src1.data[2][3]*src2.data[3][0];
}

template <typename T>
void Product441(const Matrix<4,4,T>& src1, const Matrix<4,1,T>& src2, Matrix<4,1,T>* destin)
{
    destin->data[0][0] = src1.data[0][0]*src2.data[0][0]+src1.data[0][1]*src2.data[1][0]+src1.data[0][2]*src2.data[2][0]+src1.data[0][3]*src2.data[3][0];
    destin->data[1][0] = src1.data[1][0]*src2.data[0][0]+src1.data[1][1]*src2.data[1][0]+src1.data[1][2]*src2.data[2][0]+src1.data[1][3]*src2.data[3][0];
    destin->data[2][0] = src1.data[2][0]*src2.data[0][0]+src1.data[2][1]*src2.data[1][0]+src1.data[2][2]*src2.data[2][0]+src1.data[2][3]*src2.data[3][0];
    destin->data[3][0] = src1.data[3][0]*src2.data[0][0]+src1.data[3][1]*src2.data[1][0]+src1.data[3][2]*src2.data[2][0]+src1.data[3][3]*src2.data[3][0];
}

}

#endif