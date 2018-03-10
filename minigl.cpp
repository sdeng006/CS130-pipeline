/**
 * minigl.cpp
 * -------------------------------
 * Implement miniGL here.
 *
 * You may include minigl.h and any of the standard C++ libraries.
 * No other includes are permitted.  Other preprocessing directives
 * are also not permitted.  These requirements are strictly
 * enforced.  Be sure to run a test grading to make sure your file
 * passes the sanity tests.
 *
 * The behavior of the routines your are implenting is documented here:
 * https://www.opengl.org/sdk/docs/man2/
 * Note that you will only be implementing a subset of this.  In particular,
 * you only need to implement enough to pass the tests in the suite.
 */

#include "minigl.h"
#include "vec.h"
#include "mat.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstring>
#include <iostream>

using namespace std;

/**
 * Useful data types
 */
typedef mat<MGLfloat,4,4> mat4; //data structure storing a 4x4 matrix, see mat.h
typedef mat<MGLfloat,3,3> mat3; //data structure storing a 3x3 matrix, see mat.h
typedef vec<MGLfloat,4> vec4;   //data structure storing a 4 dimensional vector, see vec.h
typedef vec<MGLfloat,3> vec3;   //data structure storing a 3 dimensional vector, see vec.h
typedef vec<MGLfloat,2> vec2;   //data structure storing a 2 dimensional vector, see vec.h

// A structure that stores a 4d position and a color
struct vertex {
    vec4 pos;
    vec3 color;
};

// A triangle structure that stores the 3 vertices
struct triangle {
    vertex a;
    vertex b;
    vertex c;
};

// list of vertex and triangles
vector<vertex> vertices;
vector<triangle> triangles;

MGLpoly_mode poly_Mode;
MGLmatrix_mode matrix_Mode;
vec3 color; 


void Rasterize_Triangle(const triangle& tri, int width,
                        int height, MGLpixel* data);

mat4 I = {{1.0f,0.0f,0.0f,0.0f,
            0.0f,1.0f,0.0f,0.0f,
            0.0f,0.0f,1.0f,0.0f,
            0.0f,0.0f,0.0f,1.0f}};
vector<mat4> Vec_projection = {I};
vector<mat4> Vec_modelview = {I};
vector<vector<float>> z_buffer;

/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}

/**
 * Read pixel data starting with the pixel at coordinates
 * (0, 0), up to (width,  height), into the array
 * pointed to by data.  The boundaries are lower-inclusive,         CHECK HERE
 * that is, a call with width = height = 1 would just read
 * the pixel at (0, 0).
 *
 * Rasterization and z-buffering should be performed when
 * this function is called, so that the data array is filled
 * with the actual pixel values that should be displayed on
 * the two-dimensional screen.
 */
void mglReadPixels(MGLsize width,
                   MGLsize height,
                   MGLpixel *data)
{
    // Allocate size, set data black, and make z_buffer all 2
    z_buffer.resize(width);
    for(unsigned int i = 0; i < width; i++) {
        z_buffer[i].resize(height);
        for(unsigned int j = 0; j < height; j++) {
            data[i + j*width] = Make_Pixel(0, 0, 0);
            z_buffer[i][j] = 2.0f;
        }
    }
    
    for(unsigned int i = 0; i < triangles.size(); i++){
        Rasterize_Triangle(triangles[i], width, height, data);
    }
    triangles.clear();
}

/**
 * helper function to get area of 3 points and return their area    DONE
 */
MGLfloat area(vec2 a, vec2 b, vec2 c){
    // Look at LAB 6 EQ. 8
    // return (0.5 * (a[0]*b[1] - a[1]*b[0] + b[0]*c[1] - b[1]*c[0] + c[0]*a[1] - c[1]*a[0]));
    return a[0] * (b[1] - c[1]) + a[1] * (c[0] - b[0]) + (b[0] * c[1] - b[1] * c[0]);
}

/**
 * my helper function=====================================================/
 */
 
// Rasterize function described in assignment file
void Rasterize_Triangle(const triangle& tria,
                        int width,
                        int height,
                        MGLpixel* data)
{
    // leave it here from mglVertex3
    triangle tri = tria;
    for(unsigned int i = 0; i < 3; i++) {
        tri.a.pos[i] = tri.a.pos[i] / tri.a.pos[3];
        tri.b.pos[i] = tri.b.pos[i] / tri.b.pos[3];
        tri.c.pos[i] = tri.c.pos[i] / tri.c.pos[3];
    }
    
    // a[0] means Ai, a[0] means Aj
    // x,y are the triangle's three vertex position 0 and 1
    vec2 A, B, C;
    A[0] = (tri.a.pos[0] + 1) * width * 0.5 - 0.5;
    A[1] = (tri.a.pos[1] + 1) * height * 0.5 - 0.5;
    B[0] = (tri.b.pos[0] + 1) * width * 0.5 - 0.5;
    B[1] = (tri.b.pos[1] + 1) * height * 0.5 - 0.5;
    C[0] = (tri.c.pos[0] + 1) * width * 0.5 - 0.5;
    C[1] = (tri.c.pos[1] + 1) * height * 0.5 - 0.5;
    
    MGLfloat xmin, xmax;
    xmin = max( min( min(A[0], B[0]), C[0]) , 0.0f);
    xmax = min( max( max(A[0], B[0]), C[0]) , (float)width);
    
    // min size: find the any of x point of the triangle or 0.0f(border)
    // max size: locate x max or width  
    for(int i = xmin; i < xmax; i++){
        for(int j = 0; j < height; j++){
            vec2 P = vec2(i, j);
            
            MGLfloat a = area(P,B,C) / area(A,B,C);
            MGLfloat b = area(A,P,C) / area(A,B,C);
            MGLfloat c = area(A,B,P) / area(A,B,C);
            
            if(a >= 0 && b >= 0 && c >=0) {
                
                MGLfloat curr_z = a * tri.a.pos[2] + b * tri.b.pos[2] + c * tri.c.pos[2];
                
                if(curr_z < z_buffer[i][j] && 
                   (-1 <= curr_z) && (curr_z <= 1)) {
                    
                    MGLfloat div = (a / tri.a.pos[3]) + (b / tri.b.pos[3]) + (c / tri.c.pos[3]);
                    MGLfloat alpha = a / tri.a.pos[3] / div;
                    MGLfloat beta = b / tri.b.pos[3] / div;
                    MGLfloat gamma = c / tri.c.pos[3] / div;
                    
                    MGLfloat Z = alpha * tri.a.pos[2] + beta * tri.b.pos[2] + gamma * tri.c.pos[2];
                    if(Z > 1 || Z < -1) continue;
                    
                    data[i + j*width] = Make_Pixel(255.0f * (tri.a.color[0] * alpha + tri.b.color[0] * beta + tri.c.color[0] * gamma), 
                                                   255.0f * (tri.a.color[1] * alpha + tri.b.color[1] * beta + tri.c.color[1] * gamma),
                                                   255.0f * (tri.a.color[2] * alpha + tri.b.color[2] * beta + tri.c.color[2] * gamma));
                    z_buffer[i][j] = curr_z;
                }
            }
        }
    }
}

// current matrix on the top, times with input
void multiple_current_matrix(mat4 a) {
    
    mat4 temp;
    
    if(matrix_Mode == MGL_PROJECTION){
        temp = Vec_projection.back() * a;
        Vec_projection.pop_back();
        Vec_projection.push_back(temp);
    }
    else if(matrix_Mode == MGL_MODELVIEW){
        temp = Vec_modelview.back() * a;
        Vec_modelview.pop_back();
        Vec_modelview.push_back(temp);
    }
    else {
        cout << "ERROR IN MULTIPLE_CURRENT_MATRIX" << endl;
        exit(0);
    }
}

/**=====================================================================**/

/**
 * Start specifying the vertices for a group of primitives,     **DONE
 * whose type is specified by the given mode.                   **DONE
 */
void mglBegin(MGLpoly_mode mode)
{
    // if(mode != MGL_TRIANGLES && mode != MGL_QUADS) {
    //     cout << "mglBegin: Error - mode not vaild" << endl;
    // }
    poly_Mode = mode;
}

/**
 * Stop specifying the vertices for a group of primitives.      **DONE
 */
void mglEnd()
{
    if(poly_Mode == MGL_TRIANGLES){
        for(unsigned int i = 0; i < vertices.size(); i = i+3){
            
            // double secure
            if(i > vertices.size()) {
                break;
            }
            
            triangle tri;
            tri.a = vertices[i];
            tri.b = vertices[i+1];
            tri.c = vertices[i+2];
    
            triangles.push_back(tri);
        }
    }
    else if(poly_Mode == MGL_QUADS){
        for(unsigned int i = 0; i < vertices.size(); i = i+4){
            
            // double secure
            if(i > vertices.size()) {
                break;
            }
            
            triangle tri1;
            triangle tri2;
            
            tri1.a = vertices[i];
            tri2.a = vertices[i];
            tri1.b = vertices[i+2];
            tri2.b = vertices[i+2];
            
            // cuzing bug? need to be checked
            tri1.c = vertices[i+1];
            tri2.c = vertices[i+3];
            
            triangles.push_back(tri1);
            triangles.push_back(tri2);
        }
    }
    vertices.clear();
}

/**
 * Specify a two-dimensional vertex; the x- and y-coordinates   DONE
 * are explicitly specified, while the z-coordinate is assumed  DONE
 * to be zero.  Must appear between calls to mglBegin() and     DONE
 * mglEnd().
 */
void mglVertex2(MGLfloat x,
                MGLfloat y)
{
    mglVertex3(x,y,0);
}

/**
 * Specify a three-dimensional vertex.  Must appear between     **DONE
 * calls to mglBegin() and mglEnd().                            
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{
    vec4 v = vec4(x, y, z, 1.0f);
    
    vec4 tempv = Vec_projection.back() * Vec_modelview.back() * v;

    vertex temp;
    temp.pos = tempv;
    // temp.pos = tempv/tempv[3];
    temp.color = color;
    
    vertices.push_back(temp);
}

/**
 * Set the current matrix mode (modelview or projection).  **DONE
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
    // if(mode != MGL_MODELVIEW && mode != MGL_PROJECTION) {
    //     cout << "mglMatrixMode: Error - mode not vaild" << endl;
    // }
    matrix_Mode = mode;
}

/**
 * Push a copy of the current matrix onto the stack for the     **DONE
 * current matrix mode.
 */
void mglPushMatrix()
{
    if(matrix_Mode == MGL_PROJECTION){
        Vec_projection.push_back(Vec_projection.back());
    }
    else if(matrix_Mode == MGL_MODELVIEW){
        Vec_modelview.push_back(Vec_modelview.back());
    }
}

/**
 * Pop the top matrix from the stack for the current matrix     **DONE
 * mode.
 */
void mglPopMatrix()
{
    if(matrix_Mode == MGL_PROJECTION) {
        Vec_projection.pop_back();
    }
    else if(matrix_Mode == MGL_MODELVIEW) {
        Vec_modelview.pop_back();
    }
}

/**
 * Replace the current matrix with the identity.                **DONE
 */
void mglLoadIdentity()
{
    if(matrix_Mode == MGL_PROJECTION){
        Vec_projection.pop_back();
        Vec_projection.push_back(I);
    }
    else if(matrix_Mode == MGL_MODELVIEW){
        Vec_modelview.pop_back();
        Vec_modelview.push_back(I);
    }
}

/**
 * Replace the current matrix with an arbitrary 4x4 matrix,  **DONE
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglLoadMatrix(const MGLfloat *matrix)
{
    mat4 temp;
    for(unsigned int i = 0; i < 16; i++) {
        temp.values[i] = matrix[i];
    }
    
    if(matrix_Mode == MGL_PROJECTION){
        Vec_projection.pop_back();
        Vec_projection.push_back(temp);
    }
    else if(matrix_Mode == MGL_MODELVIEW){
        Vec_modelview.pop_back();
        Vec_modelview.push_back(temp);
    }
}

/**
 * Multiply the current matrix by an arbitrary 4x4 matrix,  **DONE
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglMultMatrix(const MGLfloat *matrix)
{
    mat4 temp;
    for(unsigned int i = 0; i < 16; i++) {
        temp.values[i] = matrix[i];
    }
    
    multiple_current_matrix(temp);
}

/**
 * Multiply the current matrix by the translation matrix        **DONE
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
    mat4 translate = {{1.0f, 0.0f, 0.0f, 0.0f,
                   0.0f, 1.0f, 0.0f, 0.0f,
                   0.0f, 0.0f, 1.0f, 0.0f,
                   x, y, z, 1.0f}};
    
    multiple_current_matrix(translate);
}

/**
 * Multiply the current matrix by the rotation matrix           **DONE
 * for a rotation of (angle) degrees about the vector
 * from the origin to the point (x, y, z).
 */
void mglRotate(MGLfloat angle,
               MGLfloat x,
               MGLfloat y,
               MGLfloat z)
{
    float a = angle * M_PI / 180.0f;
    float s = sin(a);
    float c = cos(a);
    
    float n = sqrt(x*x + y*y + z*z);
    x = x / n;
    y = y / n;
    z = z / n;

    mat4 rotate = {{ x*x*(1.0f-c)+c,   y*x*(1-c)+z*s, x*z*(1-c)-y*s, 0.0f,
                     x*y*(1.0f-c)-z*s, y*y*(1-c)+c,   y*z*(1-c)+x*s, 0.0f,
                     x*z*(1.0f-c)+y*s, y*z*(1-c)-x*s, z*z*(1-c)+c,   0.0f,
                     0.0f, 0.0f, 0.0f, 1.0f}};
                     
    multiple_current_matrix(rotate);
}

/**
 * Multiply the current matrix by the scale matrix          **DONE
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
    mat4 scale = {{x, 0.0f, 0.0f ,0.0f,
                   0.0f, y, 0.0f, 0.0f,
                   0.0f, 0.0f, z, 0.0f,
                   0.0f, 0.0f, 0.0f, 1.0f}};
    
    multiple_current_matrix(scale);
}

/**
 * Multiply the current matrix by the perspective matrix    **DONE
 * with the given clipping plane coordinates.               **DONE
 */
void mglFrustum(MGLfloat left,
                MGLfloat right,
                MGLfloat bottom,
                MGLfloat top,
                MGLfloat near,
                MGLfloat far)
{
        /*
            https://msdn.microsoft.com/en-us/library
            /windows/desktop/dd373537(v=vs.85).aspx
        */
        
        MGLfloat A = (right + left) / (right - left);
        MGLfloat B = (top + bottom) / (top - bottom);
        MGLfloat C = -(far + near) / (far - near);
        MGLfloat D = -(2.0f * far * near) / (far - near);
        
        // CAREFUL! WARNING!
        // CAREFUL! COLUMN ORDER FIRST!
        mat4 frustum = {{2.0f*near/(right-left),0.0f,0.0f,0.0f,
                      0.0f,2.0f*near/(top-bottom),0.0f,0.0f,
                      A,B,C,-1.0f,
                      0.0f,0.0f,D,0.0f}};

        multiple_current_matrix(frustum);
}

/**
 * Multiply the current matrix by the orthographic matrix   **DONE
 * with the given clipping plane coordinates.               **DONE
 */
void mglOrtho(MGLfloat left,
              MGLfloat right,
              MGLfloat bottom,
              MGLfloat top,
              MGLfloat near,
              MGLfloat far)
{
    /* 
        https://msdn.microsoft.com/en-us/library/
        windows/desktop/dd373965(v=vs.85).aspx
    */
    
    MGLfloat tx = -(right + left) / (right - left);
    MGLfloat ty = -(top + bottom) / (top - bottom);
    MGLfloat tz = -(far + near) / (far - near);
    
    // CAREFUL! WARNING!
    // CAREFUL! COLUMN ORDER FIRST!
    mat4 ortho = {{2.0f/(right-left),0.0f,0.0f,0.0f,
                  0.0f,2.0f/(top-bottom),0.0f,0.0f,
                  0.0f,0.0f,-2.0f/(far-near),0.0f,
                  tx,ty,tz,1.0f}};
    
    multiple_current_matrix(ortho);
}

/**
 * Set the current color for drawn shapes.              **DONE
 */
void mglColor(MGLfloat red,
              MGLfloat green,
              MGLfloat blue)
{
    color[0] = red;
    color[1] = green;
    color[2] = blue;
}