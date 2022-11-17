/************************************************************
 *                                                          *
 *  Permission is hereby granted  to  any  individual   or  *
 *  institution   for  use,  copying, or redistribution of  *
 *  this code and associated documentation,  provided       *
 *  that   such  code  and documentation are not sold  for  *
 *  profit and the  following copyright notice is retained  *
 *  in the code and documentation:                          *
 *     Copyright (c) held by Dianne Cook                    *
 *  All Rights Reserved.                                    *
 *                                                          *
 *  Questions and comments are welcome, and I request       *
 *  that you share any modifications with me.               *
 *                                                          *
 *                Dianne Cook                               *
 *             dicook@iastate.edu                           *
 *                                                          *
 ************************************************************/

#define PRECISION1 32768
#define PRECISION2 16384
/*#define PI 3.1415926535897932*/
#define MIN(x,y) ( (x) < (y) ? (x) : (y) )
#define MAX(x,y) ((x)>(y)?(x):(y))
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define MAXINT 2147483647
#define ASCII_TEXT_BORDER_WIDTH 4
#define MAXHIST 100
#define STEP0 0.01
#define FORWARD 1
#define BACKWARD -1
#define PROJ_DIM 5
#define True 1
#define False 0

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
//#include <cmat.h>
//#include <math.h>
//#include <cmath>
#include <cmath> 
#include <vector>
#include <random>
#include <ctime>
#include <iostream>
#include <omp.h>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <pthread.h>
#include <queue>
#include <limits>
#include "numeric"
#include <limits>
#include <time.h>
#include <bits/stdc++.h>




typedef struct {
	float x, y, z;
} fcoords;

typedef struct {
	long x, y, z;
} lcoords;

typedef struct {
	int x, y, z;
} icoords;

typedef struct {
	float min, max;
} lims;


/* grand tour history */
typedef struct hist_rec {
  struct hist_rec *prev, *next;
  float *basis[3];
  int pos;
} hist_rec;




//#pragma once



static double PYTHAG(double a, double b);
int dsvd(float **a, int m, int n, float *w, float **v);


/* 
 * svdcomp - SVD decomposition routine. 
 * Takes an mxn matrix a and decomposes it into udv, where u,v are
 * left and right orthogonal transformation matrices, and d is a 
 * diagonal matrix of singular values.
 *
 * This routine is adapted from svdecomp.c in XLISP-STAT 2.1 which is 
 * code from Numerical Recipes adapted by Luke Tierney and David Betz.
 *
 * Input to dsvd is as follows:
 *   a = mxn matrix to be decomposed, gets overwritten with u
 *   m = row dimension of a
 *   n = column dimension of a
 *   w = returns the vector of singular values of a
 *   v = returns the right orthogonal transformation matrix
*/



 
static double PYTHAG(double a, double b)
{
    double at = fabs(a), bt = fabs(b), ct, result;

    if (at > bt)       { ct = bt / at; result = at * sqrt(1.0 + ct * ct); }
    else if (bt > 0.0) { ct = at / bt; result = bt * sqrt(1.0 + ct * ct); }
    else result = 0.0;
    return(result);
}


int dsvd(float **a, int m, int n, float *w, float **v)
{
  /*
    Matrix a: [1..m][1..n]
      replaced by U.
    vector w: singular values
    matrix v: matrix V [1..n][1..n]
  */
    int flag, i, its, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;
    double *rv1;
  
    if (m < n) 
    {
        fprintf(stderr, "#rows must be > #cols \n");
        return(0);
    }
  
    rv1 = (double *)malloc((unsigned int) n*sizeof(double));

/* Householder reduction to bidiagonal form */
    for (i = 0; i < n; i++) 
    {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m) 
        {
            for (k = i; k < m; k++) 
                scale += fabs((double)a[k][i]);
            if (scale) 
            {
                for (k = i; k < m; k++) 
                {
                    a[k][i] = (float)((double)a[k][i]/scale);
                    s += ((double)a[k][i] * (double)a[k][i]);
                }
                f = (double)a[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][i] = (float)(f - g);
                if (i != n - 1) 
                {
                    for (j = l; j < n; j++) 
                    {
                        for (s = 0.0, k = i; k < m; k++) 
                            s += ((double)a[k][i] * (double)a[k][j]);
                        f = s / h;
                        for (k = i; k < m; k++) 
                            a[k][j] += (float)(f * (double)a[k][i]);
                    }
                }
                for (k = i; k < m; k++) 
                    a[k][i] = (float)((double)a[k][i]*scale);
            }
        }
        w[i] = (float)(scale * g);
    
        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != n - 1) 
        {
            for (k = l; k < n; k++) 
                scale += fabs((double)a[i][k]);
            if (scale) 
            {
                for (k = l; k < n; k++) 
                {
                    a[i][k] = (float)((double)a[i][k]/scale);
                    s += ((double)a[i][k] * (double)a[i][k]);
                }
                f = (double)a[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][l] = (float)(f - g);
                for (k = l; k < n; k++) 
                    rv1[k] = (double)a[i][k] / h;
                if (i != m - 1) 
                {
                    for (j = l; j < m; j++) 
                    {
                        for (s = 0.0, k = l; k < n; k++) 
                            s += ((double)a[j][k] * (double)a[i][k]);
                        for (k = l; k < n; k++) 
                            a[j][k] += (float)(s * rv1[k]);
                    }
                }
                for (k = l; k < n; k++) 
                    a[i][k] = (float)((double)a[i][k]*scale);
            }
        }
        anorm = MAX(anorm, (fabs((double)w[i]) + fabs(rv1[i])));
    }
  
    /* accumulate the right-hand transformation */
    for (i = n - 1; i >= 0; i--) 
    {
        if (i < n - 1) 
        {
            if (g) 
            {
                for (j = l; j < n; j++)
                    v[j][i] = (float)(((double)a[i][j] / (double)a[i][l]) / g);
                    /* double division to avoid underflow */
                for (j = l; j < n; j++) 
                {
                    for (s = 0.0, k = l; k < n; k++) 
                        s += ((double)a[i][k] * (double)v[k][j]);
                    for (k = l; k < n; k++) 
                        v[k][j] += (float)(s * (double)v[k][i]);
                }
            }
            for (j = l; j < n; j++) 
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }
  
    /* accumulate the left-hand transformation */
    for (i = n - 1; i >= 0; i--) 
    {
        l = i + 1;
        g = (double)w[i];
        if (i < n - 1) 
            for (j = l; j < n; j++) 
                a[i][j] = 0.0;
        if (g) 
        {
            g = 1.0 / g;
            if (i != n - 1) 
            {
                for (j = l; j < n; j++) 
                {
                    for (s = 0.0, k = l; k < m; k++) 
                        s += ((double)a[k][i] * (double)a[k][j]);
                    f = (s / (double)a[i][i]) * g;
                    for (k = i; k < m; k++) 
                        a[k][j] += (float)(f * (double)a[k][i]);
                }
            }
            for (j = i; j < m; j++) 
                a[j][i] = (float)((double)a[j][i]*g);
        }
        else 
        {
            for (j = i; j < m; j++) 
                a[j][i] = 0.0;
        }
        ++a[i][i];
    }

    /* diagonalize the bidiagonal form */
    for (k = n - 1; k >= 0; k--) 
    {                             /* loop over singular values */
        for (its = 0; its < 30; its++) 
        {                         /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--) 
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm) 
                {
                    flag = 0;
                    break;
                }
                if (fabs((double)w[nm]) + anorm == anorm) 
                    break;
            }
            if (flag) 
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++) 
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm) 
                    {
                        g = (double)w[i];
                        h = PYTHAG(f, g);
                        w[i] = (float)h; 
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++) 
                        {
                            y = (double)a[j][nm];
                            z = (double)a[j][i];
                            a[j][nm] = (float)(y * c + z * s);
                            a[j][i] = (float)(z * c - y * s);
                        }
                    }
                }
            }
            z = (double)w[k];
            if (l == k) 
            {                  /* convergence */
                if (z < 0.0) 
                {              /* make singular value nonnegative */
                    w[k] = (float)(-z);
                    for (j = 0; j < n; j++) 
                        v[j][k] = (-v[j][k]);
                }
                break;
            }
            if (its >= 30) {
                free((void*) rv1);
                fprintf(stderr, "No convergence after 30,000! iterations \n");
                return(0);
            }
    
            /* shift from bottom 2 x 2 minor */
            x = (double)w[l];
            nm = k - 1;
            y = (double)w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
          
            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++) 
            {
                i = j + 1;
                g = rv1[i];
                y = (double)w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < n; jj++) 
                {
                    x = (double)v[jj][j];
                    z = (double)v[jj][i];
                    v[jj][j] = (float)(x * c + z * s);
                    v[jj][i] = (float)(z * c - x * s);
                }
                z = PYTHAG(f, h);
                w[j] = (float)z;
                if (z) 
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++) 
                {
                    y = (double)a[jj][j];
                    z = (double)a[jj][i];
                    a[jj][j] = (float)(y * c + z * s);
                    a[jj][i] = (float)(z * c - y * s);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = (float)x;
        }
    }
    //delete((void*) rv1);
    return(1);
}

//#pragma once


using namespace std;
#define SORT_ON_X 0
#define SORT_ON_Y 1
#define SORT_ON_Z 2

namespace gs
{

  class KdTree;
  class Point;


  class Point
  {
  public:
    Point();
    Point(float x, float y, float z);
    inline Point(Point &rhs)
    {
      this->pos[0] = rhs.pos[0];
      this->pos[1] = rhs.pos[1];
      this->pos[2] = rhs.pos[2];
    }

    inline Point operator+(const Point& rhs)
    {
      Point r;
      r.pos[0] = pos[0] + rhs.pos[0];
      r.pos[1] = pos[1] + rhs.pos[1];
      r.pos[2] = pos[2] + rhs.pos[2];
      return r;
    }
    inline Point operator-(const Point& rhs)
    {
      Point r;
      r.pos[0] = pos[0] - rhs.pos[0];
      r.pos[1] = pos[1] - rhs.pos[1];
      r.pos[2] = pos[2] - rhs.pos[2];
      return r;
    }
    inline Point& operator=(const Point& rhs)
    {
      this->pos[0] = rhs.pos[0];
      this->pos[1] = rhs.pos[1];
      this->pos[2] = rhs.pos[2];

      return *this;
    }

    float pos[3];
  };

  class KdTree
  {
  public:
    KdTree(std::vector<Point*> &pointCloud);
    KdTree(std::vector<Point*> &pointCloud, int start, int end, int sortOn);
    virtual ~KdTree();
    void build(std::vector<Point*> &pointCloud, int start, int end, int sortOn);
    
    bool isLeaf();
    float split();
    KdTree* getChild(Point* searchPoint);
    float nodeX();
    float nodeY();
    float nodeZ();

    void search(Point* p, Point* result);
    void radiusSearch(Point* p, float* searchRadius, Point* result);

  private:
    gs::KdTree* __children[2];
    Point* __node;
    int __sortOn;
    int __start;
    int __end;

    void insertionSort(std::vector<Point*> &pointCloud, int start, int end, int sortOn);
    void mergeSort(std::vector<Point*> &pointCloud, int start, int end);
    void merge(std::vector<Point*> &pointCloud, int left, int mid, int right);

    int getNextSortOn(int sortOn);
    Point** tempArray;
  };
}



using namespace std;
gs::Point::Point()
{
  this->pos[0] = 0.0;
  this->pos[1] = 0.0;
  this->pos[2] = 0.0;
}
gs::Point::Point(float x, float y, float z)
{
  this->pos[0] = x;
  this->pos[1] = y;
  this->pos[2] = z;
}

gs::KdTree::KdTree(std::vector<Point*> &pointCloud)
{
  tempArray = new Point*[pointCloud.size()];
  build(pointCloud, 0, pointCloud.size(), 0);

  delete tempArray;
}
gs::KdTree::KdTree(std::vector<Point*> &pointCloud, int start, int end, int sortOn)
{
  build(pointCloud, start, end, sortOn);
}
gs::KdTree::~KdTree()
{
  delete __children[0];
  delete __children[1];
  delete __node;
}

void gs::KdTree::build(std::vector<Point*> &pointCloud, int start, int end, int sortOn)
{
  __children[0] = nullptr;
  __children[1] = nullptr;
  __node = nullptr;

  __sortOn = sortOn;
  __start = start;
  __end = end;

  if (end - start == 1)
  {
    __node = new Point(pointCloud[start]->pos[0], pointCloud[start]->pos[1], pointCloud[start]->pos[2]);
    return;
  }

  //insertionSort(pointCloud, start, end, sortOn);
  mergeSort(pointCloud, start, end);

  int numPoints = (end - start);
  int mid = (int)floor((float)numPoints*0.5) + start;

  __node = new Point(pointCloud[mid]->pos[0], pointCloud[mid]->pos[1], pointCloud[mid]->pos[2]);

  int numPointsHigh = (end - mid);
  int numPointsLow = (mid - start);

  if (numPointsHigh > 0)
  {
    __children[1] = new KdTree(pointCloud, mid, end, getNextSortOn(sortOn));
  }

  if (numPointsLow > 0)
  {
    __children[0] = new KdTree(pointCloud, start, mid, getNextSortOn(sortOn));
  }
}

void gs::KdTree::insertionSort(std::vector<Point*> &pointCloud, int start, int end, int sortOn)
{
  for (int i = start; i < end; i++)
  {
    for (int j = i + 1; j < end; j++)
    {
      if (pointCloud[j]->pos[sortOn] < pointCloud[i]->pos[sortOn])
      {
        Point* temp = pointCloud[i];
        pointCloud[i] = pointCloud[j];
        pointCloud[j] = temp;
      }
    }
  }
}

int gs::KdTree::getNextSortOn(int sortOn)
{
  switch (sortOn)
  {
  case SORT_ON_X:
    return SORT_ON_Y;
    break;
  case SORT_ON_Y:
    return SORT_ON_Z;
    break;
  case SORT_ON_Z:
    return SORT_ON_X;
    break;
  }
}

bool gs::KdTree::isLeaf()
{
  if (__children[0] == nullptr && __children[1] == nullptr)
    return true;
  
  return false;
}

float gs::KdTree::split()
{
  return __node->pos[__sortOn];
}

gs::KdTree* gs::KdTree::getChild(Point* searchPoint)
{
  if (searchPoint->pos[__sortOn] >= split())
  {
    return __children[1];
  }
  else
  {
    return __children[0];
  }
}

float gs::KdTree::nodeX()
{
  return __node->pos[0];
}
float gs::KdTree::nodeY()
{
  return __node->pos[1];
}
float gs::KdTree::nodeZ()
{
  return __node->pos[2];
}

void gs::KdTree::radiusSearch(Point* p, float* radius, Point* result)
{
  if (isLeaf())
  {
    float d = sqrt(pow(__node->pos[0] - p->pos[0], 2.0) + pow(__node->pos[1] - p->pos[1], 2.0) + pow(__node->pos[2] - p->pos[2], 2.0));
    if (d < *radius)
    {
      *radius = d;
      result->pos[0] = __node->pos[0];
      result->pos[1] = __node->pos[1];
      result->pos[2] = __node->pos[2];
      return;
    }
  }
  else
  {
    if (std::abs(__node->pos[__sortOn] - p->pos[__sortOn]) < *radius)
    {
      __children[0]->radiusSearch(p, radius, result);
      __children[1]->radiusSearch(p, radius, result);
    }
    else
    {
      if (p->pos[__sortOn] >= __node->pos[__sortOn])
      {
        __children[1]->radiusSearch(p, radius, result);
      }
      else
      {
        __children[0]->radiusSearch(p, radius, result);
      }
    }
  }
}

void gs::KdTree::search(Point* p, Point* result)
{
  // get closets node
  KdTree* tree = this;
  while (!tree->isLeaf())
  {
    tree = tree->getChild(p);
  }
  result->pos[0] = tree->nodeX();
  result->pos[1] = tree->nodeY();
  result->pos[2] = tree->nodeZ();

  float radius = sqrt(pow(p->pos[0] - result->pos[0], 2.0) + pow(p->pos[1] - result->pos[1], 2.0) + pow(p->pos[2] - result->pos[2], 2.0));

  radiusSearch(p, &radius, result);
}

void gs::KdTree::mergeSort(std::vector<Point*> &pointCloud, int start, int end)
{
  int mid;
  if (start > end)
  {
    mid = (int)floor((end + start)*0.5);
    mergeSort(pointCloud, start, mid);
    mergeSort(pointCloud, mid + 1, end);

    merge(pointCloud, start, mid + 1, end);
  }
}

void gs::KdTree::merge(std::vector<Point*> &pointCloud, int left, int mid, int right)
{
  int i, leftEnd, numElements, tempPos;

  leftEnd = mid;
  tempPos = 0;
  numElements = right - left;

  while (left < leftEnd && mid <= right)
  {
    if (pointCloud[left]->pos[__sortOn] <= pointCloud[mid]->pos[__sortOn])
    {
      tempArray[tempPos] = pointCloud[left];
      tempPos++;
      left++;
    }
    else
    {
      tempArray[tempPos] = pointCloud[mid];
      tempPos++;
      mid++;
    }
  }

  while (left < leftEnd)
  {
    tempArray[tempPos] = pointCloud[left];
    left++;
    tempPos++;
  }
  while (mid <= right)
  {
    tempArray[tempPos] = pointCloud[mid];
    mid++;
    tempPos++;
  }

  for (int i = tempPos - 1; i >= 0; i--)
  {
    pointCloud[right] = tempArray[i];
    right--;
  }
}


/*
Iterative Closest Point

Performs the iterative closest point algorithm. A point cloud 'dynamicPointCloud' is transformed such that
it best "matches" the point cloud 'staticPointCloud'.

this program uses the method outlined in [1].

1. for each point in the dynamic point cloud, find the closest point in the static point cloud.
2. solve for an affine transform which minimizes the distance between all dynamic points and their respective static points.
3. transform the dynamic points.
4. goto -> 1 until stopping conditions are met.

affine transforms are solved using SVD.

@author Greg Smith 2017

[1] Arun, K. Somani, Thomas S. Huang, and Steven D. Blostein. "Least-squares fitting of two 3-D point sets." 
  IEEE Transactions on pattern analysis and machine intelligence 5 (1987): 698-700.

*/

//#pragma once



namespace gs
{

  /*
  void icp
    transforms the point cloud 'dynamicPointCloud' such that it best matches 'staticPointCloud'

  @param std::vector<Point*> dynamicPointCloud : point cloud to be rotated and translated to match 'staticPointCloud'
  @param std::vector<Point*> staticPointCloud : reference point cloud.
  */
  void icp(std::vector<Point*> &dynamicPointCloud, std::vector<Point*> &staticPointCloud);

  inline void computeCloudMean(std::vector<Point*> &cloud, gs::Point* mean)
  {
    mean->pos[0] = 0.0;
    mean->pos[1] = 0.0;
    mean->pos[2] = 0.0;
    for (int i = 0; i < cloud.size(); i++)
    {
      for (int j = 0; j < 3; j++)
      {
        mean->pos[j] += cloud[i]->pos[j];
      }
    }
    mean->pos[0] = mean->pos[0] / (float)cloud.size();
    mean->pos[1] = mean->pos[1] / (float)cloud.size();
    mean->pos[2] = mean->pos[2] / (float)cloud.size();
  }

  inline void clearTranslation(float* translation)
  {
    translation[0] = 0.0;
    translation[1] = 0.0;
    translation[2] = 0.0;
  }

  inline void clearRotation(float* rotation)
  {
    rotation[0] = 0.0;
    rotation[1] = 0.0;
    rotation[2] = 0.0;

    rotation[3] = 0.0;
    rotation[4] = 0.0;
    rotation[5] = 0.0;

    rotation[6] = 0.0;
    rotation[7] = 0.0;
    rotation[8] = 0.0;
  }

  inline void clearMatrix(float* mat)
  {
    mat[0] = 0.0;
    mat[1] = 0.0;
    mat[2] = 0.0;

    mat[3] = 0.0;
    mat[4] = 0.0;
    mat[5] = 0.0;

    mat[6] = 0.0;
    mat[7] = 0.0;
    mat[8] = 0.0;
  }

  inline void rotate(gs::Point* p, float* rotationMatrix, gs::Point* result)
  {
    result->pos[0] = p->pos[0] * rotationMatrix[0] + p->pos[1] * rotationMatrix[1] + p->pos[2] * rotationMatrix[2];
    result->pos[1] = p->pos[0] * rotationMatrix[3] + p->pos[1] * rotationMatrix[4] + p->pos[2] * rotationMatrix[5];
    result->pos[2] = p->pos[0] * rotationMatrix[6] + p->pos[1] * rotationMatrix[7] + p->pos[2] * rotationMatrix[8];
  }

  inline void translate(gs::Point* p, float* translationVector, gs::Point* result)
  {
    result->pos[0] = p->pos[0] + translationVector[0];
    result->pos[1] = p->pos[1] + translationVector[1];
    result->pos[2] = p->pos[2] + translationVector[2];
  }

  inline void outerProduct(gs::Point* a, gs::Point* b, float* mat)
  {
    mat[0] = a->pos[0] * b->pos[0];
    mat[1] = a->pos[0] * b->pos[1];
    mat[2] = a->pos[0] * b->pos[2];

    mat[3] = a->pos[1] * b->pos[0];
    mat[4] = a->pos[1] * b->pos[1];
    mat[5] = a->pos[1] * b->pos[2];

    mat[6] = a->pos[2] * b->pos[0];
    mat[7] = a->pos[2] * b->pos[1];
    mat[8] = a->pos[2] * b->pos[2];
  }

  inline void matrixMult(float* a, float* b, float* result)
  {
    result[0] = a[0] * b[0] + a[1] * b[3] + a[2] * b[6];
    result[1] = a[0] * b[1] + a[1] * b[4] + a[2] * b[7];
    result[2] = a[0] * b[2] + a[1] * b[5] + a[2] * b[8];

    result[3] = a[3] * b[0] + a[4] * b[3] + a[5] * b[6];
    result[4] = a[3] * b[1] + a[4] * b[4] + a[5] * b[7];
    result[5] = a[3] * b[2] + a[4] * b[5] + a[5] * b[8];

    result[6] = a[6] * b[0] + a[7] * b[3] + a[8] * b[6];
    result[7] = a[6] * b[1] + a[7] * b[4] + a[8] * b[7];
    result[8] = a[6] * b[2] + a[7] * b[5] + a[8] * b[8];
  }

  inline void transpose(float* a)
  {
    float temp;

    temp = a[1];
    a[1] = a[3];
    a[3] = temp;

    temp = a[2];
    a[2] = a[6];
    a[6] = temp;

    temp = a[5];
    a[5] = a[7];
    a[7] = temp;
  }

  inline void addMatrix(float* a, float* b, float* result)
  {
    result[0] = a[0] + b[0];
    result[1] = a[1] + b[1];
    result[2] = a[2] + b[2];

    result[3] = a[3] + b[3];
    result[4] = a[4] + b[4];
    result[5] = a[5] + b[5];

    result[6] = a[6] + b[6];
    result[7] = a[7] + b[7];
    result[8] = a[8] + b[8];
  }

  inline float error(Point* ps, Point* pd, float* r, float* t)
  {
    Point res;
    rotate(pd, r, &res);
    float err = pow(ps->pos[0] - res.pos[0] - t[0],2.0);
    err += pow(ps->pos[1] - res.pos[1] - t[1], 2.0);
    err += pow(ps->pos[2] - res.pos[2] - t[2], 2.0);
    return err;
  }

  inline void copyMatToUV(float* mat, float** result)
  {
    result[0][0] = mat[0];
    result[0][1] = mat[1];
    result[0][2] = mat[2];

    result[1][0] = mat[3];
    result[1][1] = mat[4];
    result[1][2] = mat[5];

    result[2][0] = mat[6];
    result[2][1] = mat[7];
    result[2][2] = mat[8];
  }

  inline void copyUVtoMat(float** mat, float* result)
  {
    result[0] = mat[0][0];
    result[1] = mat[0][1];
    result[2] = mat[0][2];

    result[3] = mat[1][0];
    result[4] = mat[1][1];
    result[5] = mat[1][2];

    result[6] = mat[2][0];
    result[7] = mat[2][1];
    result[8] = mat[2][2];
  }
}


int index_to_ckeck = 1237;
#define point_dim 3
using namespace std;
using namespace gs;
int num_calcs;
//double basis[3] = {26.19372829,12.88486764,7.69639231};
//double basis[3] = {4.41840962,2.42806299,1.49673601};
//double basis[3] = {6.89697076,3.52467804,2.04844848};
//double basis[3] = {1.87220823,0.76542574,0.41296428};
//double basis[3] = {1.09387599, -0.96979575 , 0.94387252};
//double basis[3] = {1,1,-1};
//double basis[3] = {1.51664444, 0.66704231, 0.51779803};
//double basis[3] = {4.07571969, 1.99038848, 1.32392432};
//double basis[3] = {1.12409238, -0.90933214,  0.94841502};
// double basis[3] = {0.11346113, -0.71007125,  1.34066972};
//double basis[3] = {0.29276996, -0.6530365 ,  1.33276279};
//double basis[3] = {-0.12370598, -0.10563619, -1.31032003};
//double basis[3]  = {
 //double basis[3] = {0.11346113, -0.71007125,  1.34066972}
//double basis[3] = {1,1,1};
//double basis[3] = {-0.51475309, -0.29207345,  0.58999927};
//double basis[3]  = {0.46084765, 0.18950231, 0.68975714};
//double basis[3] = {-0.06889851,  0.00933936, -0.71927585};
//double basis[3] = {0.02665797,  0.02011677, -0.11730251};
//double basis[3] = {-0.47216308, -0.06279136,  1.00655682};
//double basis[3] = {-0.10831018, -0.65576258,  2.53081975};
//double basis[3] = {-0.63842603,  0.02794823,  2.83228815};
//double basis[3] = {0.05957336, -0.35304188,  0.96837945};
//double basis[3] = {0,-1,-1};
//double basis[3] = {0.39195903, 0.57994389, 0.72309869};
//double basis[3] = {1,1,1};
//double basis[3] = {0.01105764, 0.50933852, 0.74604819};
double basis[3] = {1,1,1};
double basis_size = sqrt(basis[0] * basis[0] + basis[1] * basis[1] + basis[2] * basis[2]);
bool iter_save_points = true;
bool main_save_points = false;
bool random_init_points = false;
int num_initial_points =200;
int svd_size = 150;


vector<vector <double>> init_queries;
vector<vector <double>> init_references(num_initial_points);
vector <double> temp_ref(3);
bool use_proposed;
enum dist_metric
{
    Modified_Manhattan,
    Euclidean,
    Manhattan
};


void shuffle_array(int arr[], int n)
{
 
    // To obtain a time-based seed
    unsigned seed = 0;
 
    // Shuffling our array
    shuffle(arr, arr + n,
            default_random_engine(seed));
 
    // Printing our array
    //for (int i = 0; i < n; ++i)
        //cout << arr[i] << " ";
    //cout << endl;
}



double calc_distance (double *v1, vector<double>v2, dist_metric type)
{
    //cout<<"distance of ";
    //print_vector_double(v1);
    //print_vector_double(v2);
    if (type == Modified_Manhattan)
    {
        double sum1 = 0;
        double sum2 = 0;
        for(int i = 0; i<point_dim;i++)
            sum1+=v1[i];
        for(int i = 0; i<point_dim;i++)
            sum2+=v2[i];
        return (abs(sum2 - sum1));
    }
    else
    {
        double sum = 0;
        for(int i = 0; i<point_dim;i++)
        {
            if (type==Euclidean)
            {
              sum+= pow(abs(v1[i] - v2[i]), 2);
            }
            if (type==Manhattan)
            sum+= abs(v1[i] - v2[i]);
        }
        double result = sum;
        if (type == Euclidean)
            result = sqrt(result);
        //cout<<result;
        
        return(result);
        }
}

class Frame{
    //later: change to private
public:
    vector<vector<double>> row_data;
    int num_points; 
    double ** data;
    Frame(string file_adress, int max_points=0)
    {
      ifstream fin(file_adress, ios::binary);
      fin.seekg(0, ios::end);
      size_t num_elements = fin.tellg() / sizeof(double);
      cout<<file_adress<<file_adress<< num_elements<<endl;
      if (max_points!=0) num_elements = (max_points*point_dim);
      num_points = num_elements/point_dim;
      cout<<num_points<<endl;
      fin.seekg(0, ios::beg);
      row_data = vector<vector<double>> (num_points , vector<double> (point_dim, 0));
      for (int c = 0 ; c<num_points; c++)
      {
          if (c%200 == 0) 
              {cout<<c<<endl;}
          fin.read(reinterpret_cast<char*>(&row_data[c][0]), point_dim*sizeof(double));
      }
      allocate_data();
      }
      void allocate_data()
      {
          //allocating 
          double * temp = new double [num_points*(point_dim+1)];
          data = new double*[num_points];
          for (int i = 0 ; i < num_points;i++)
          {
              data[i] = (temp+i*(point_dim+1));
          }
    }
    Frame(vector<Point*>* point_data)
    {
      num_points = point_data->size();
      row_data = vector<vector<double>> (point_data->size() , vector<double> (point_dim, 0));
      for (int i = 0 ; i <point_data->size();i++)
      {
        for (int j = 0 ; j < point_dim; j++)
        {
          row_data[i][j] = (*point_data)[i]->pos[j];
        }
      }
      allocate_data();
    }
    void create_reference_data()
    {
    vector<double> reference_projected(num_points);
    for (int i =0 ; i<num_points;i++)
    {
        reference_projected[i] = 0;
        for (int j = 0; j<point_dim;j++)
        {
        reference_projected[i] += ((row_data[i][j]) * basis[j]);
        }
    }
    vector<int> reference_order(num_points);
    iota(reference_order.begin(),reference_order.end(),0); //Initializing
    sort( reference_order.begin(),reference_order.end(), [&](int i,int j){return reference_projected[i]<reference_projected[j];} );
    double sum;
    for (int i = 0; i<num_points;i++)
    {
        sum = 0;
        for (int j = 0; j<point_dim;j++)
        {
            data[i][j] = row_data[reference_order[i]][j];
            sum += ((data[i][j])*basis[j]);
        }
        data[i][point_dim] = sum;
    }
    }

    void create_query_data(vector<double>* query_projected)
    {
        for (int i =0 ; i<num_points;i++)
    {
        (*query_projected)[i] = 0;
        for (int j = 0; j<point_dim;j++)
        {
        (*query_projected)[i] += row_data[i][j];
        }
    }
    }
};


int binary_search (double ** reference, double query, int begin, int end)
{
  //cout<<"query: "<<query<<" begin: "<<begin<<" end: "<<end<<endl;
    int length = end - begin+1;
    int end_orig = end;
    int middle_index = (begin + end) / 2;
    double middle = reference[(int)((begin + end) / 2)][point_dim];
    while (end >= begin)
    {
        middle_index = (begin + end) / 2;
        middle = reference[(int)((begin + end) / 2)][point_dim];

        if (query == middle) 
        {
          //*(nearest_projected) = reference[middle_index][point_dim];
            return (middle_index);
        }
        else if (query > middle) 
        {
            begin = middle_index+1;
        }
        else if(query < middle) 
            {

                end = middle_index-1;
            }
        }
        double diff1 = abs(query - middle);
        double diff2;
        double diff3;
        if (middle_index < end_orig)
        {
            diff2 = abs(query - reference[(middle_index+1)][point_dim]);
        }
        else {
            diff2 =numeric_limits<double>::max() ;
        }
        if (middle_index > 0)
        {
            diff3 = abs(query - reference[middle_index-1][point_dim]);
        }
        else
        {
            diff3 = numeric_limits<double>::max();
        }
        if ((diff1 <= diff2) && (diff1 <= diff3))  {
        //*(nearest_projected) = reference[middle_index][point_dim];
        return(middle_index);
        }
        else if ((diff2 <= diff1) && (diff2 <= diff3))
        {
          //*(nearest_projected) = reference[middle_index+1][point_dim];
            return(middle_index+1);
        }
        else if((diff3 <= diff2) && (diff3 <= diff1)) 
        {
        //*(nearest_projected) = reference[middle_index-1][point_dim];
        return(middle_index-1);
        }
}

int exact_knn_projected(const Frame* reference,vector<double>query, double query_projected, int nearest_index, int K, int row, int num_ref_points, double * NN_points, int result_index, vector<int>&search_range, int bucket_size = 0)
{
  if (iter_save_points) init_queries.push_back(query);
  search_range[0] = nearest_index;
  search_range[1] = nearest_index;
  //cout<<endl<<"init_references: "<<init_references.size();
  //cout<<endl<<"a new call for "<<nearest_index<<endl;
    double NN_dist = calc_distance(reference->data[nearest_index], query, Euclidean);
    
    /*if (dist_to_prev_NN != 0)
    {
      if (dist_to_prev_NN < NN_dist)
        NN_dist = dist_to_prev_NN;
    }*/
    int NN_index = nearest_index;    
    double dist;
    bool right_progress = true;
    bool left_progress = true;
    bool original_threshold_visited = false;
    bool bidirectional_cont = true;
    int right_candidate = nearest_index;
    int left_candidate = nearest_index;
    int next;
    bool search_cont = true;
    if (left_candidate == -1)
    {
        left_candidate++;
        bidirectional_cont = false;
        left_progress = false;
        search_range[0] = left_candidate;
    }
    if (right_candidate == num_ref_points)
    {
        right_candidate--;
        bidirectional_cont = false;
        right_progress = false;
        search_range[1] = right_candidate;
    }
    if (abs(reference->data[right_candidate][point_dim] - query_projected) <abs(reference->data[left_candidate][point_dim] - query_projected))
            next = right_candidate;
    else
        next = left_candidate;
    while(search_cont)
    {
      //cout<<endl<<"in processing "<<nearest_index<<" "<< left_progress<<" , "<<right_progress<<" ----------: "<<next;
        dist = calc_distance(reference->data[next], query, Euclidean);
        num_calcs++;
        if ((bucket_size) && (num_calcs>bucket_size))
        {
          search_cont = false;
        }
        if (dist < NN_dist)
        {
          NN_index = next;
          NN_dist = dist;
        }

        if  ((abs( reference->data[next][point_dim] - query_projected ) > (basis_size*NN_dist)))
          {search_cont = false;}

        if (left_progress && right_progress)
        {
                        if (left_candidate == -1)
            {
                left_candidate++;
                bidirectional_cont = false;
                left_progress = false;
                search_range[0] = left_candidate;
            }
            if (right_candidate == num_ref_points)
            {
                right_candidate--;
                bidirectional_cont = false;
                right_progress = false;
                search_range[1] = right_candidate;
            }
             if (abs(reference->data[right_candidate][point_dim] - query_projected) <abs(reference->data[left_candidate][point_dim] - query_projected))
            {next = right_candidate;
                right_candidate++;
                search_range[1] =right_candidate;
            }
            else
            {
            next = left_candidate;
            left_candidate--;
            search_range[0] = left_candidate;
        }
        }
        else
        {
        
        if (!left_progress)
        {
            right_candidate++;
            next = right_candidate;
            search_range[1] = right_candidate;
            if (right_candidate == num_ref_points)
            {
                right_candidate--;
                search_range[1] = right_candidate;
                bidirectional_cont = false;
                right_progress = false;
            }

        }
        if (!right_progress)
        {
            left_candidate--;
            search_range[0] = left_candidate;
            next = left_candidate;
            if (left_candidate == -1)
            {
                left_candidate++;
                bidirectional_cont = false;
                left_progress = false;
                search_range[0] = left_candidate;
            }
        }
        if ((!left_progress) && (!right_progress))          
        {
          search_cont = 0;
          break;
        }

    }
    }
    for (int i = 0 ; i < point_dim; i++)
    {
      NN_points[result_index *point_dim+i] = reference->data[NN_index][i];
      temp_ref[i]= reference->data[NN_index][i];
    }
    if (iter_save_points)
    {init_references[result_index] = temp_ref;}
    //cout<<endl<<num_calcs<<endl;
    return NN_index;
}





/*int exact_knn_projected_(const Frame* reference,vector<double>query, double query_projected, int nearest_index, int K, int row, int num_ref_points)
{
    
    double NN_dist = calc_distance(reference->data[nearest_index], query, Euclidean);
    int NN_index = nearest_index;
    double max_dist = dist;
    
    
    if (right_arrow<num_ref_points)
    {
    while( abs( reference->data[right_arrow][point_dim] - query_projected ) <= (sqrt(3)*max_dist)    )
    {
        dist = calc_distance(reference->data[right_arrow], query, Euclidean);
        if (dist < max_dist)
        {
          NN_index = right_arrow;
          NN_dist = dist;
          max_dist = dist;
        }
        right_arrow++;
        if (right_arrow == num_ref_points)
            break;
      }
  }
if (left_arrow>0)
{
        while(abs(reference->data[left_arrow][point_dim] - query_projected) <= (sqrt(3)*max_dist))
    {
        dist = calc_distance(reference->data[left_arrow], query, Euclidean);
        if (dist < max_dist)
        {
          NN_index = left_arrow;
          NN_dist = dist;
          max_dist = dist;
        }
        left_arrow--;
        if (left_arrow<0) break;
    }
}
  int c = 0;

}
*/




void gs::icp(std::vector<Point*> &dynamicPointCloud, std::vector<Point*> &staticPointCloud)
{
  

  //fout.close();
  
  int test_sum_calcs = 0;
  float rotationMatrix[9];
  float translation[3];
  int number_of_in_ranges = 0;

  std::vector<Point*> staticPointCloudCopy;
  
  gs::Point dynamicMid(0.0,0.0,0.0);
  gs::Point staticMid(0.0,0.0,0.0);

  // copy the static point cloud
  for (int i = 0; i < staticPointCloud.size(); i++)
  {
    Point* pCopy = new Point(staticPointCloud[i]->pos[0], staticPointCloud[i]->pos[1], staticPointCloud[i]->pos[2]);
    staticPointCloudCopy.push_back(pCopy);
    //cout<<endl<<"Point "<<i<<" copied";
  }
      /*cout<<endl<<"staticPointCloudCopy "<< staticPointCloudCopy.size();
  cout<<endl<<"dynamicPointCloud "<< dynamicPointCloud.size();
  for (int i =0 ; i < 10; i++)
  {
    cout<<endl<<staticPointCloudCopy[i]->pos[0]<<" , "<<staticPointCloudCopy[i]->pos[1]<<" , " <<staticPointCloudCopy[i]->pos[2]<<endl;    
    cout<<endl<<dynamicPointCloud[i]->pos[0]<<" , "<<dynamicPointCloud[i]->pos[1]<<" , " <<dynamicPointCloud[i]->pos[2]<<endl;    

  }

  exit(0);*/

  // create the kd tree
  KdTree* tree = new KdTree(staticPointCloudCopy);
  Frame reference(&staticPointCloudCopy);
  double create_reference_data_time = -omp_get_wtime();
  reference.create_reference_data();
  create_reference_data_time += omp_get_wtime();
  cout<<endl<<create_reference_data_time<<endl;
  

  for (int i = 0 ; i < 5; i++)
  {
    for (int j =0 ; j < point_dim+1; j++)
    {
      //cout<<endl<<"i: "<<i<<" j: "<<j<<" "<<reference.data[i][j]<<endl;
    }
  }
  size_t numDynamicPoints = dynamicPointCloud.size();

  computeCloudMean(staticPointCloud, &staticMid);
  computeCloudMean(dynamicPointCloud, &dynamicMid);

  // initialize the translation vector
  clearTranslation(translation);

  // initialize the rotation matrix
  clearRotation(rotationMatrix);

  const int maxIterations = 30;
  //const int numRandomSamples = dynamicPointCloud.size();
  const int numRandomSamples = svd_size;
  const float eps = 1e-8;
  gs::Point p;
  gs::Point x;
  int j1 = 0;
  int point_cloud_size = dynamicPointCloud.size();

  //vector<int> random_indices1;
  int randSample;
  int * random_indices2 = new int[point_cloud_size];
  for (int i = 0 ; i < point_cloud_size; i++)
  {
    random_indices2[i] = i;
  }
  //shuffle_array(random_indices2, point_cloud_size);
  /*
  for (int i = 0; i < svd_size; i++)
    {
      randSample = std::rand() % dynamicPointCloud.size();
      //random_indices1.push_back(randSample);
      random_indices1.push_back(i);
      //for (int j2 = 0 ; j2<1000;j2++)
      //j1++;
      //cout<<"index: " <<random_indices1[random_indices1.size()-1]<<endl;
    }*/

    if (random_init_points)
    {


    //vector<vector <double>> init_first_pair(num_initial_points);
    //vector<vector <double>> init_second_pair(num_initial_points);
    //vector <double> temp_pair;
    vector<vector<double>> init_first_pair(num_initial_points , vector<double> (point_dim, 0));
    vector<vector<double>> init_second_pair(num_initial_points , vector<double> (point_dim, 0));
    vector<double> init_distances(num_initial_points);
    int randSample1;
    int randSample2;
    //vector<int> init_indices(num_initial_points*2) [2];
    for (int i = 0; i < num_initial_points; i++)
    {
      randSample1 = std::rand() % dynamicPointCloud.size();
      randSample2 = std::rand() % dynamicPointCloud.size();
      for (int j = 0 ; j < point_dim; j++)
      {
        init_first_pair[i][j] = reference.data[randSample1][j];
        init_second_pair[i][j] = reference.data[randSample2][j];
      }
      double sum = 0;
      for (int j = 0 ; j < point_dim; j++)
      {
        sum+= pow((reference.data[randSample1][j]) - (reference.data[randSample2][j]),2);
      }
      init_distances[i] = sqrt(sum);
    }
    ofstream random_points("/home/behnam/phd/Research/ICP/init_points/random_points");
    for(int i = 0; i < num_initial_points; i++)
      {
        cout<<endl<<"init first pair:  "<<i<<" 'th point: ["<<init_first_pair[i][0]<<" , "<<init_first_pair[i][1]<<" , "<<init_first_pair[i][2]<<"]"<<endl;
        cout<<endl<<"init second pair: "<<i<<" 'th point: ["<<init_second_pair[i][0]<<" , "<<init_second_pair[i][1]<<" , "<<init_second_pair[i][2]<<"]"<<endl;
        cout<<endl<<"distance        : "<<init_distances[i]<<endl;
      }


      for(int i = 0; i < num_initial_points; i++)
      {
        for (int j = 0 ; j <point_dim;j++)
        {
          random_points<<init_first_pair[i][j]<<endl;
          //cout<<endl<<"reference "<<i<<" "<<j<<" wrote"<<endl;
        }
        for (int j = 0 ; j <point_dim;j++)
        {
          random_points<<init_second_pair[i][j]<<endl;
        }
      }


      
exit(0);
    }



  float cost = 1.0;
  std::srand(std::time(0));

  gs::Point qd;
  gs::Point qs;

  float U[9];
  float w[9];
  float sigma[3];
  float V[9];

  float** uSvd = new float*[3];
  float** vSvd = new float*[3];
  uSvd[0] = new float[3];
  uSvd[1] = new float[3];
  uSvd[2] = new float[3];

  vSvd[0] = new float[3];
  vSvd[1] = new float[3];
  vSvd[2] = new float[3];
  int num_iterations = 0;
  int sum_calcs = 0;
  double iter_time = 0;
  double kd_tree_search_time;
  double sum_kd_search_time = 0;
  double sum_proposed_time = 0;
  //vector<double> NN_projected  (dynamicPointCloud.size());
  //to do : free this space
  double * NN_projected =  new double[dynamicPointCloud.size()];
  //double NN_projected[39000];
  double * NN_points = new double [dynamicPointCloud.size() * point_dim];
  //double NN_points[39000*3];
  //double * displacement = new double[dynamicPointCloud.size()];
  
  
  double x_prev,y_prev,z_prev;
  double itrative_threshold;
  double prev_NN_projected;
  double p_projected;
  int switch_iteration = -1;
  int num_rights = 0;
  int num_lefts = 0;
  use_proposed = false;
  vector<int> search_range;
  search_range.push_back(0);
  search_range.push_back(1);
  int * past_nearest_indices = new int[point_cloud_size];
  float sum_ratio = 0;
  for (int iter = 0; iter < maxIterations && abs(cost) > eps; iter++)
  {
    sum_ratio =0;
    number_of_in_ranges = 0;
    num_lefts = 0;
    num_rights = 0;
    if (iter > switch_iteration)
      use_proposed = true;
    iter_time = -omp_get_wtime();
    sum_calcs = 0;
    num_iterations++;
    cost = 0.0;
    //clearRotation(rotationMatrix);
    clearMatrix(U);
    clearMatrix(V);
    clearMatrix(w);
    computeCloudMean(dynamicPointCloud, &dynamicMid);
    int nearest_index = 0;
    if (!use_proposed)
    {
      for (int i = 0; i < numRandomSamples; i++)
      {
        int randSample = std::rand() % dynamicPointCloud.size();
        // sample the dynamic point cloud
        p = *dynamicPointCloud[randSample];

        // get the closest point in the static point cloud
        kd_tree_search_time = -omp_get_wtime();
        tree->search(&p, &x);
        kd_tree_search_time += omp_get_wtime();
        sum_kd_search_time += kd_tree_search_time;
        qd = p - dynamicMid;
        qs = x - staticMid;
        outerProduct(&qs, &qd, w);
        addMatrix(w, U, U);
        cost += error(&x, &p, rotationMatrix, translation);
      }
    }
    else
    {
      if (iter == maxIterations-1)
        svd_size = 500;
      for (int i = 0; i < svd_size; i++)
      {
        if (iter == maxIterations-1)
        //cout<<"i:"<<i<<endl;
        //if (iter < 2) cout<<endl<<"processing point: "<<i;
        int randSample = std::rand() % dynamicPointCloud.size();
        // sample the dynamic point cloud
        //p = *dynamicPointCloud[randSample];
        p = *dynamicPointCloud[random_indices2[i]];
        p_projected = (p.pos[0]*basis[0]) + (p.pos[1]*basis[1]) + (p.pos[2]*basis[2]);
        
        vector <double> p_vec(3);
        for (int d = 0 ; d < point_dim;d++)
        {
          p_vec[d] = p.pos[d];
        }
        double binary_search_time = -omp_get_wtime();
        int nearest_index = binary_search (reference.data,p_projected, 0, reference.num_points-1);
        binary_search_time += omp_get_wtime();
        /*if (i ==random_indices1[5])
        {
          cout<<endl<<"************************"<<endl;
          cout<<"in iteration: "<<iter<<" nearest_index is: "<<nearest_index;
          cout<<endl<<"************************"<<endl;
        }*/

        num_calcs = 0;
        if (i<num_initial_points) {iter_save_points = true;}
        else {iter_save_points = false;}
        double exact_search_time = -omp_get_wtime();
        int bucket_size = 0;
            exact_knn_projected(&reference,p_vec,p_projected, nearest_index, 1, 0,reference.num_points, NN_points, i, search_range, bucket_size);
            //cout<<endl<<num_calcs<<endl;
            //cout<<endl<<"iteration: "<<iter<<" in processing point: "<<i<<" "<<search_range[0]<<" to "<<search_range[1]<<" : " <<search_range[1] - search_range[0]<< " "<<num_calcs;
            if (iter>0)
            {
              //cout<<" "<<past_nearest_indices[i]<<" ";
              if ((past_nearest_indices[i] >= search_range[0]) && (past_nearest_indices[i] <= search_range[1]))
              {
                //cout<<"Yes"<<endl;
                number_of_in_ranges++;
              }
              else
              {
                //cout<<"No"<<endl;
              }

            }
            else
            {
              //cout<<endl;
            }
            past_nearest_indices[i] = nearest_index;
            exact_search_time += omp_get_wtime();
            //cout<<endl<<"point: "<<i<<" ratio: "<<binary_search_time / exact_search_time;
            sum_ratio += (binary_search_time / exact_search_time);
            //cout<<endl<<"in iteration: "<<iter<<" for point " <<i<<" : "<<nearest_index<<" , " <<output_index[i]<<" direction: "<< output_index[i] - nearest_index<<endl;
            sum_calcs += num_calcs;
            kd_tree_search_time = -omp_get_wtime();
        x.pos[0] = NN_points[i * point_dim];
        x.pos[1] = NN_points[i * point_dim + 1];
        x.pos[2] = NN_points[i * point_dim + 2];
      
        qd = p - dynamicMid;
        qs = x - staticMid;
        outerProduct(&qs, &qd, w);
        addMatrix(w, U, U);
        cost += error(&x, &p, rotationMatrix, translation);
        /*if (i<num_initial_points)
        {
          cout<<"point"<<i<<"added "<<init_queries.size()<<" "<<init_references.size()<<endl;;
          init_queries[i] = p_vec;
          temp_ref[0] = x.pos[0];
          temp_ref[1] = x.pos[1];
          temp_ref[2] = x.pos[2];
          init_references[i] = temp_ref;
        }*/

        //***********************************
        //saving the points when processing num_initial_points's point
        //***********************************

        if ((iter == 0) &&main_save_points && (i == num_initial_points-1))
        {
        ofstream fout("/home/behnam/phd/Research/ICP/init_points/points", ios::app);
        for (int i3 = 0 ;i3 < num_initial_points; i3++ )
        {
          cout<<endl<<i3<<":"<<endl;
          cout<<endl<<"inja: reference point: ";
          for (int i4 = 0 ; i4 < point_dim; i4++)
          {
            cout<<init_references[i3][i4]<<" , ";
            fout<<init_references[i3][i4]<<endl;
            cout<<"neveshte shod";

          }
          cout<<endl<<"inja: query point: ";
          for (int i4 = 0 ; i4 < point_dim; i4++)
          {
            cout<<init_queries[i3][i4]<<" , ";
            fout<<init_queries[i3][i4]<<endl;
          }
          cout<<endl;
        }
        fout.close();
        exit(0);
        }


    }
        //shuffle_array(random_indices2, point_cloud_size);
    cout<<endl<<"iteration level ratio: "<<sum_ratio/svd_size<<endl;
  }
    copyMatToUV(U, uSvd);
    dsvd(uSvd, 3, 3, sigma, vSvd);
    copyUVtoMat(uSvd, U);
    copyUVtoMat(vSvd, V);

    transpose(V);
    matrixMult(U, V, rotationMatrix);

    gs::Point t(0.0, 0.0, 0.0);
    rotate(&dynamicMid, rotationMatrix, &t);
    translation[0] = staticMid.pos[0] - t.pos[0];
    translation[1] = staticMid.pos[1] - t.pos[1];
    translation[2] = staticMid.pos[2] - t.pos[2];
    //update the point cloud
    for (int i = 0; i < point_cloud_size; i++)
    {
      x_prev = dynamicPointCloud[i]->pos[0];
      y_prev = dynamicPointCloud[i]->pos[1];
      z_prev = dynamicPointCloud[i]->pos[2];
      rotate(dynamicPointCloud[i], rotationMatrix, &p);
      translate(&p, translation, dynamicPointCloud[i]);     
      //dist_to_prev_NN[i] = sqrt(pow(NN_points[(i*point_dim)] - dynamicPointCloud[i]->pos[0],2) + pow(NN_points[(i*point_dim+1)]-dynamicPointCloud[i]->pos[1],2) +pow(NN_points[i*point_dim+2] - dynamicPointCloud[i]->pos[2],2));
      //displacement[i] = sqrt(pow(x_prev - dynamicPointCloud[i]->pos[0],2) + pow(y_prev-dynamicPointCloud[i]->pos[1],2) +pow(z_prev - dynamicPointCloud[i]->pos[2],2));
    }
    cout<<endl<<"in iteration "<<iter<<" average number of calcs: "<<((float)sum_calcs / svd_size)<<endl;
    cout<<endl<<"in iteration "<<iter<<" cost is: "<<cost<<endl;
    cout<<endl<<"in iteration "<<iter<<" number_of_in_ranges is: "<<(float)number_of_in_ranges/svd_size<<endl;
    test_sum_calcs += ((float)sum_calcs / svd_size);

    //cout<<endl<<"num_lefts: "<<num_lefts<<" num_rights: "<<num_rights<<endl;
    //cout<<endl<<"in iteration "<<iter<<" time: "<<iter_time+omp_get_wtime()<<endl;
    //cout<<endl<<"in iteration "<<iter<<" cost: "<<cost<<endl;
    //cout<<endl<<endl<<endl<<"first iteration is ok:"<<endl;
    //cout<<endl<<"and the results are:"<<endl;
    cout<<endl<<init_references.size()<<" "<<init_queries.size()<<endl;

    
      
    
    /*if (save_points)
    {
    for (int i3 = 0 ;i3 < num_initial_points; i3++ )
    {
      cout<<endl<<i3<<":"<<endl;
      cout<<endl<<"reference point: ";
      for (int i4 = 0 ; i4 < point_dim; i4++)
      {
        cout<<init_references[i3][i4]<<" , ";
        fout<<init_references[i3][i4]<<endl;

      }
      cout<<endl<<"query point: ";
      for (int i4 = 0 ; i4 < point_dim; i4++)
      {
        cout<<init_queries[i3][i4]<<" , ";
        fout<<init_queries[i3][i4]<<endl;
      }
      cout<<endl;

    }
    
    //exit(0);
    }*/
    //char temp;
    //cin>>temp;
  }

  //cout<<endl<<"sum_kd_search_time"<<sum_kd_search_time<<endl;

  cout<<endl<<num_iterations<<" maxIterations ";
  

  staticPointCloudCopy.clear();
  delete tree;
  
  delete[] uSvd[0]; 
  delete[] uSvd[1];
  delete[] uSvd[2];
  delete[] uSvd;

  delete[] vSvd[0];
  delete[] vSvd[1];
  delete[] vSvd[2];
  delete[] vSvd;
  cout<<endl<<"test_sum_calcs: "<<test_sum_calcs<<endl;
  cout<<endl<<endl<<"----deleting-----"<<endl<<endl;
  delete [] NN_projected;
  delete [] NN_points; 
  //delete [] displacement;
  //delete [] dist_to_prev_NN ;



}



/*

whole kd-tree search:

total_error0.00118645
whole_time: 0.277023
sh: 1: pause: not found



after 15 iterations:

Alignment Error: 0.00000 
total_error5.14011e-07
whole_time: 0.903092
sh: 1: pause: not found



after 17 iterations:
Alignment Error: 0.00002 
total_error0.000408422
whole_time: 0.690345
sh: 1: pause: not found


double basis[3] = {1,1,-1};
in iteration 0 average number of calcs: 15760.5

double basis[3] = {0.11346113, -0.71007125,  1.34066972};
in iteration 0 average number of calcs: 22121.8

double basis[3] = {0.29276996, -0.6530365 ,  1.33276279};
22216.2


double basis[3] = {-0.12370598, -0.10563619, -1.31032003};
10046.3

double basis[3] = {1,1,1};
12795.5

double basis[3] = {-0.51475309, -0.29207345,  0.58999927};
16699.8

double basis[3] = {-0.06889851,  0.00933936, -0.71927585};
14020.4

double basis[3] = {0.02665797,  0.02011677, -0.11730251};
21420.4

four times of expand:
in iteration 99 average number of calcs: 3014.09

in iteration 99 cost is: 206.711

200 15000

100 maxIterations 
test_sum_calcs: 379976


----deleting-----

Alignment Error: 103796.72656 

each_test_time: 550.876
total_error103797

whole_time: 550.895


when shufflig in each iteration: eval size= 50000
513412

519036

5639.71875


150000  not shuffling every time
in iteration 49 average number of calcs: 6029.03
in iteration 49 cost is: 1.57389e+06


150000, shuffling every time
in iteration 49 average number of calcs: 6084.35
in iteration 49 cost is: 1.54776e+06




kollan beddone shuffle (svd = 150 and 1500):
155906
yek bar shuffle (svd = 150 and 1500):
426753



49988
i:149989
i:149990
i:149991
i:149992
i:149993
i:149994
i:149995
i:149996
i:149997
i:149998
i:149999

in iteration 49 average number of calcs: -4393.93

in iteration 49 cost is: 6.71558e+08

in iteration 49 number_of_in_ranges is: 0.00348

200 10000

50 maxIterations 
test_sum_calcs: 259609


----deleting-----

Alignment Error: 6305.30127 

each_test_time: 7533.38
total_error6305.3

whole_time: 7533.38
sh: 1: pause: not found
*/



#include <iostream>
#include <fstream>
#include <omp.h>
using namespace gs;
using namespace std;
#define point_dim 3


float total_error = 0;
/*
Create a test box of points.
*/

void print_vector_2D (vector<vector<float>>input){
    for (int i = 0; i< input.size();i++)
    {
        for(int j = 0; j<input[0].size();j++)
        {
            cout<<input[i][j]<<" ";
        }
        cout<<endl;
    }

}
class Frame1{
    //later: change to private
public:
    vector<vector<float>> data;
    Frame1(string file_adress)
    {
    ifstream fin(file_adress);
    
    bool finished = false;
    string a1 = "start";
    int counter = 0;
    while(!finished)
    {
        //cout<<a1<<" has been read" ;
        getline(fin, a1, ',');      
        if (a1[0]!='e')
        {
            //cout<<stof(a1)<<" has been read"<<endl ;
            counter++;
            if ((counter%100) == 0)
            {
                //cout<<counter<<endl;
            }
            //cout<<stof(a1)<<" has been read"<<endl ;
            data.push_back(vector<float>(point_dim));       
            data[data.size()-1][0] = stof(a1);
            for (int c = 1 ;c<point_dim;c++)
            {
                getline(fin, a1, ',');      
                data[data.size()-1][c] = stof(a1);
            }
        }
        else finished = true;
        
    }
    }
};






void createPoints_backup(std::vector<Point*>& points, Frame1 * reference)
{
  for (int i = 0 ; i < reference->data.size(); i++)
  {
    points.push_back(new Point(reference->data[i][0], reference->data[i][1], reference->data[i][2]));   
  }
}


void createPoints1(std::vector<Point*>& points)
{
  points.push_back(new Point(3.0f, -2.0f, -7.0f));
  points.push_back(new Point(1.0f, -1.0f, -1.0f));
  points.push_back(new Point(0.0f, 1.0f, -5.0f));
  points.push_back(new Point(2.6f, 1.7f, 9.8f));
}

/*
Apply an afine transofrm to a point cloud
*/
void applyAffineTransform(std::vector<Point*>& points, float* rotationMatrix, float* translation)
{
  Point pRot;
  for (int i = 0; i < points.size(); i++)
  {
    pRot.pos[0] = rotationMatrix[0] * points[i]->pos[0] + rotationMatrix[1] * points[i]->pos[1] + rotationMatrix[2] * points[i]->pos[2] + translation[0];
    pRot.pos[1] = rotationMatrix[3] * points[i]->pos[0] + rotationMatrix[4] * points[i]->pos[1] + rotationMatrix[5] * points[i]->pos[2] + translation[1];
    pRot.pos[2] = rotationMatrix[6] * points[i]->pos[0] + rotationMatrix[7] * points[i]->pos[1] + rotationMatrix[8] * points[i]->pos[2] + translation[2];

    *points[i] = pRot;
  }
}

/*
ICP is used to minimize the distance between 2 point clouds.

This example computes a point cloud of a unit box (Static point cloud)
and a second unit box with a linear transform applied (dynamic point cloud).
ICP is used to transform the dynamic point cloud to best match the static point cloud.
*/

using namespace std;

void expand_point_cloud (vector <Point*>* input)
{
  //double minx, maxx, miny, maxy, minx, maxz;
  double min_cordinates[point_dim];
  double max_cordinates[point_dim];
  for (int dim =0 ; dim < point_dim; dim++)
  {
    min_cordinates[dim] = (*input)[0]->pos[dim];
    max_cordinates[dim] = (*input)[0]->pos[dim];
  }

  for (int i = 1 ; i<input->size(); i++)
  {
    for(int j = 0 ; j < point_dim; j++)
    {
      if ((*input)[i]->pos[j] < min_cordinates[j])
      {
        min_cordinates[j] = (*input)[i]->pos[j];
      }
      if ((*input)[i]->pos[j] > max_cordinates[j])
      {
        max_cordinates[j] = (*input)[i]->pos[j];
      }
    }
  }
  cout<<endl<<min_cordinates[0]<<" to "<<max_cordinates[0]<<endl;
  cout<<endl<<min_cordinates[1]<<" to "<<max_cordinates[1]<<endl;
  cout<<endl<<min_cordinates[2]<<" to "<<max_cordinates[2]<<endl;
  int index =0;
  /*while(1)
  {
    cout<<"index: "; cin>>index;
    cout<<endl<<"["<<(*input)[index]->pos[0]<<" , "<<(*input)[index]->pos[1]<<" , "<<(*input)[index]->pos[2]<<"]"<<endl;
  }
  exit(0);*/
  int input_size = input->size();
  for (int i = 0 ; i < input_size; i++)
  {
    input->push_back(new Point
      ((*input)[i]->pos[0]+max_cordinates[0],(*input)[i]->pos[1],(*input)[i]->pos[2]));
  }
  

}


int read_frame_size(const char * file_adress)
{
    cout<<endl<<"reding size of file "<<file_adress;
    FILE *f1 = fopen(file_adress, "rb");
    fseek(f1, 0L, SEEK_END); // seek to end of file
    int frame_size = ftell(f1);    //get current file pointer

    frame_size /= (point_dim * sizeof(float));
    //vector<vector<float>> result (frame_size , vector<float> (point_dim, 0));
    fseek(f1, 0, SEEK_SET);
    return frame_size;
    
}

void createPoints(vector<Point*>& PointCloud ,const char * file_adress, int frame_size)
{  
    FILE *f1 = fopen(file_adress, "rb");
    fseek(f1, 0, SEEK_SET);
    float * row_data = new float[frame_size* point_dim];
    cout<<endl<<"---------- reading file ------------"<<endl;
    fread((row_data), sizeof(float),frame_size*point_dim, f1);
    cout<<"file has been read";
    cout<<endl<<row_data[(frame_size - 1)*3 + 2];

    for (int i = 0 ; i < frame_size; i++)
    {
        
        
        if (i % 20000 == 0)
        {
            cout<<endl<<i;
        //cout<<endl<<i;
        //cout<<endl<<"point "<<i<<"from "<<frame_size<<" has been read";
        
        }
        //cout<<endl<<row_data[i*3]<<" , "<<row_data[i*3 + 1]<<" , "<< row_data[i*3+2];
        PointCloud[i] = new Point;
        
        //cout<<endl<<"point "<<i<<" dim[0] set to "<<row_data[i*3];
        PointCloud[i]->pos[0] = row_data[i*3];
        //cout<<endl<<"point "<<i<<" dim[1] set to "<<row_data[i*3+1];
        PointCloud[i]->pos[1] = row_data[i*3+1];
        //cout<<endl<<"point "<<i<<" dim[2] set to "<<row_data[i*3+2];
        PointCloud[i]->pos[2] = row_data[i*3+2];
        //cout<<endl<<"point "<<i<<"from "<<frame_size<<" has been read";
        
        //cout<<endl<<"point "<<i<<"from "<<frame_size<<" has been read";
        //points.push_back(new Point(row_data[i*3], row_data[i*3 + 1], row_data[i*3+2]));   
    }
    fclose(f1);
    //for (int i = 0 ; i < reference->data.size(); i++)
    //{
        //
    //}

/*
    for (int i = 0 ; i < 30; i++)
      {
    cout<<endl<<i<<" : "<<reference.data[i/3][i%3];
    }
*/
/*
    
    for (int i = 0 ; i < size; i++)
    {
        cout<<endl<<i<<" : "<<row_data[i];
    }


  for (int i = 0 ; i < reference->data.size(); i++)
  {
    points.push_back(new Point(reference->data[i][0], reference->data[i][1], reference->data[i][2]));   
  }
  */
}


void icpExample(const char * ref_file, const char * query_file)
{
  
  //create a static box point cloud used as a reference.
    int reference_size = read_frame_size(ref_file);
    int query_size = read_frame_size(query_file);
    vector<Point*> dynamicPointCloud(query_size);
    vector<Point*> staticPointCloud(reference_size);
    createPoints(staticPointCloud, ref_file, reference_size);
    cout<<endl<<"first frame has been read";
    createPoints(dynamicPointCloud,  query_file, query_size);
  
  //expand_point_cloud(&staticPointCloud);
  //expand_point_cloud(&staticPointCloud);
  //expand_point_cloud(&staticPointCloud);
  //expand_point_cloud(&staticPointCloud);
  //expand_point_cloud(&staticPointCloud);

  //create a dynamic box point cloud.
  //this point cloud is transformed to match the static point cloud.
  
  //expand_point_cloud(&dynamicPointCloud);
  //expand_point_cloud(&dynamicPointCloud);
  //expand_point_cloud(&dynamicPointCloud);
  //expand_point_cloud(&dynamicPointCloud);
  //expand_point_cloud(&dynamicPointCloud);
  //cout<<endl<<"expansion done"<<endl;

  //apply an artitrary rotation and translation to the dynamic point cloud to misalign the point cloud.
  //float rotation[] = { 1.0f, 0.0f, 0.0f,  0.0f, 0.70710678f, -0.70710678f,  0.0f, 0.70710678f, 0.70710678f };
  //float translation[] = { -0.75f, 0.5f, -0.5f };
  //applyAffineTransform(dynamicPointCloud, rotation, translation);

  /*printf("Static point Cloud: \n");
  for (int i = 0; i < staticPointCloud.size(); i++)
  {
    printf("%0.2f, %0.2f, %0.2f \n", staticPointCloud[i]->pos[0], staticPointCloud[i]->pos[1], staticPointCloud[i]->pos[2]);
  }
  printf("\n");*/

  /*printf("Dynamic point Cloud: \n");
  for (int i = 0; i < dynamicPointCloud.size(); i++)
  {
    printf("%0.2f, %0.2f, %0.2f \n", dynamicPointCloud[i]->pos[0], dynamicPointCloud[i]->pos[1], dynamicPointCloud[i]->pos[2]);
  }
  printf("\n");*/

  //use iterative closest point to transform the dynamic point cloud to best align the static point cloud.

  icp(dynamicPointCloud, staticPointCloud);

  /*printf("Dynamic point Cloud Transformed: \n");
  for (int i = 0; i < dynamicPointCloud.size(); i++)
  {
    printf("%0.2f, %0.2f, %0.2f \n", dynamicPointCloud[i]->pos[0], dynamicPointCloud[i]->pos[1], dynamicPointCloud[i]->pos[2]);
  }
  printf("\n");*/

  float alignmentError = 0.0f;
  for (int i = 0; i < min(dynamicPointCloud.size(), staticPointCloud.size()); i++)
  {
    //cout<<endl<<"for point "<<i<<" it is ok "<<endl;
    alignmentError += pow(dynamicPointCloud[i]->pos[0] - staticPointCloud[i]->pos[0], 2.0f);
    alignmentError += pow(dynamicPointCloud[i]->pos[1] - staticPointCloud[i]->pos[1], 2.0f);
    alignmentError += pow(dynamicPointCloud[i]->pos[2] - staticPointCloud[i]->pos[2], 2.0f);
  }

  alignmentError /= (float)dynamicPointCloud.size();
  total_error += alignmentError;

  printf("Alignment Error: %0.5f \n", alignmentError);
}

void save_frame1(Frame1 & input_frame, const char * file_adress)
{
    FILE *f = fopen(file_adress, "w");
    for (int i=0 ; i < input_frame.data.size(); i++)
    {
        for (int j =0 ;j <point_dim; j++)
        {
            fwrite(&(input_frame.data[i][j]), sizeof(float), 1, f);
            //fwrite(test_data,              sizeof(float), data_size, f);
        }
    }
    fclose(f);
}

void temp_show_file(const char * file_adress, int size)
{
    float * row_data = new float[size];
    FILE *f1 = fopen(file_adress, "rb");
    fseek(f1, 0, SEEK_SET);
    fread((row_data), sizeof(float), size, f1);
    fclose(f1);
    cout<<endl<<"-------------"<<endl;
    for (int i = 0 ; i < size; i++)
    {
        cout<<endl<<i<<" : "<<row_data[i];
    }
}

int main()
{
  //Frame reference("reformed_dataset/PC1.txt");
  //Frame1 reference("reformed_dataset/0_gr.txt");
  //Frame1 query("reformed_dataset/1_gr.txt");
  //cout<<endl<<reference.data.size();
  
  //--------------- run once to create a frame--------------:
  //const char file_adress1[] = "reformed_dataset/0_gr_test.data";
  /*
  Frame1 reference("reformed_dataset/rad_and_black_0.txt");
  Frame1 query("reformed_dataset/rad_and_black_1.txt");
  const char file_adress1[] = "reformed_dataset/rad_and_black_0.data";
  const char file_adress2[] = "reformed_dataset/rad_and_black_1.data";
  save_frame1(reference, file_adress1);
  save_frame1(query, file_adress2);
  exit(0);*/
  
  //temp_show_file(file_adress, 30);
  /*
  for (int i = 0 ; i < 30; i++)
  {
    cout<<endl<<i<<" : "<<reference.data[i/3][i%3];
  }
  */

  //exit(0);

  //Frame1 reference("reformed_dataset/rad_and_black_0.txt");
  //Frame reference("reformed_dataset/rad_and_black_0.txt");
  //Frame1 query("reformed_dataset/rad_and_black_1.txt");
  //query.create_query_data();
  //Frame1 query("reformed_dataset/1_gr.txt");
  const char ref_file[]   = "reformed_dataset/rad_and_black_0.data";
  const char query_file[] = "reformed_dataset/rad_and_black_1.data";
  //cout<<reference.data.size()<<endl;
  //cout<<query.data.size();
  double whole_time = -omp_get_wtime();
  double each_test_time = 0;
  int num_tests = 1;
  for (int i = 0 ;i< num_tests; i++)
  {
    cout<<"test number "<<i<<endl;
    each_test_time = -omp_get_wtime();
    icpExample(ref_file, query_file);
    each_test_time += omp_get_wtime();
    cout<<endl<<"each_test_time: "<<each_test_time<<endl;
  }
  total_error /= num_tests;
  cout<<"total_error"<<total_error<<endl;
  whole_time += omp_get_wtime();
  cout<<endl<<"whole_time: "<<whole_time<<endl;
  system("pause");
  return 0;
}