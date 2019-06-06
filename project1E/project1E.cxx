#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <math.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <algorithm>
#include <vtkDataSetWriter.h>


using std::cerr;
using std::endl;

double ceil441(double f)
{
    return ceil(f-0.00001);
}

double floor441(double f)
{
    return floor(f+0.00001);
}


vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void
WriteImage(vtkImageData *img, const char *filename)
{
   std::string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
   writer->SetInputData(img);
   writer->SetFileName(full_filename.c_str());
   writer->Write();
   writer->Delete();
}

class Triangle
{
  public:
      double         X[3];
      double         Y[3];
      double         Z[3];
      double colors[3][3];
      double normals[3][3];
  // would some methods for transforming the triangle in place be helpful?
};

class Screen
{
  public:
      unsigned char   *buffer;
      int width, height;

  // would some methods for accessing and setting pixels be helpful?
};




std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1e_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    cerr << "Done reading" << endl;
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();

    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
    double *color_ptr = var->GetPointer(0);
    //vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    //float *color_ptr = var->GetPointer(0);
    vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
    float *normals = n->GetPointer(0);
    std::vector<Triangle> tris(numTris);
    vtkIdType npts;
    vtkIdType *ptIds;
    int idx;
    for (idx = 0, cells->InitTraversal() ; cells->GetNextCell(npts, ptIds) ; idx++)
    {
        if (npts != 3)
        {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        double *pt = NULL;
        pt = pts->GetPoint(ptIds[0]);
        tris[idx].X[0] = (pt[0]+10)*50.0;
        tris[idx].Y[0] = (pt[1]+10)*50.0;
        tris[idx].Z[0] = (pt[2]-10)*0.05;
        tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
        tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
        tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
        pt = pts->GetPoint(ptIds[1]);
        tris[idx].X[1] = (pt[0]+10)*50.0;
        tris[idx].Y[1] = (pt[1]+10)*50.0;
        tris[idx].Z[1] = (pt[2]-10)*0.05;
        tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
        tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
        tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
        pt = pts->GetPoint(ptIds[2]);
        tris[idx].X[2] = (pt[0]+10)*50.0;
        tris[idx].Y[2] = (pt[1]+10)*50.0;
        tris[idx].Z[2] = (pt[2]-10)*0.05;
        tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
        tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
        tris[idx].normals[2][2] = normals[3*ptIds[2]+2];
        
        // 1->2 interpolate between light blue, dark blue
        // 2->2.5 interpolate between dark blue, cyan
        // 2.5->3 interpolate between cyan, green
        // 3->3.5 interpolate between green, yellow
        // 3.5->4 interpolate between yellow, orange
        // 4->5 interpolate between orange, brick
        // 5->6 interpolate between brick, salmon
        double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
        double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
        unsigned char RGB[8][3] = { { 71, 71, 219 },
                                    { 0, 0, 91 },
                                    { 0, 255, 255 },
                                    { 0, 128, 0 },
                                    { 255, 255, 0 },
                                    { 255, 96, 0 },
                                    { 107, 0, 0 },
                                    { 224, 76, 76 }
                                  };
        for (int j = 0 ; j < 3 ; j++)
        {
            float val = color_ptr[ptIds[j]];
            int r;
            for (r = 0 ; r < 7 ; r++)
            {
                if (mins[r] <= val && val < maxs[r])
                    break;
            }
            if (r == 7)
            {
                cerr << "Could not interpolate color for " << val << endl;
                exit(EXIT_FAILURE);
            }
            double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
            tris[idx].colors[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
            tris[idx].colors[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
            tris[idx].colors[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
        }
    }

    return tris;
}
struct LightingParameters
{
    LightingParameters(void)
    {
         lightDir[0] = -0.6;
         lightDir[1] = 0;
         lightDir[2] = -0.8;
         Ka = 0.3;
         Kd = 0.7;
         Ks = 5.3;
         alpha = 7.5;
    };


    double lightDir[3]; // The direction of the light source
    double Ka;           // The coefficient for ambient lighting.
    double Kd;           // The coefficient for diffuse lighting.
    double Ks;           // The coefficient for specular lighting.
    double alpha;        // The exponent term for specular lighting.
};

LightingParameters lp;

double interp(double f_a, double f_b, int a, int b, double x){
	double t = (x-a)/(b-a);
        return f_a + (t * (f_b - f_a));
}

//L is from the lightparameters
//N is from the normal
double applyShading(double x, double y, double z){
        double r[3];
        double l =  lp.lightDir[0]*lp.lightDir[0] + lp.lightDir[1]*lp.lightDir[1] + lp.lightDir[2]*lp.lightDir[2];
        double n = x*x + y*y + z*z;
        double dif = abs(l*n);
        double dotProd = (lp.lightDir[0]*x + lp.lightDir[1]*y + lp.lightDir[2]*z);
        double _x = lp.lightDir[0] - 2*(dotProd * x);
        double _y = lp.lightDir[1] - 2*(dotProd * y);
        double _z = lp.lightDir[2] - 2*(dotProd * z);
        r[0] = _x;
        r[1] = _y;
        r[2] = _z;
        double spec = std::max(0.0,lp.Ks * pow((r[0] * x + r[1] * y + r[2] * z) , lp.alpha));
        return lp.Ka + (lp.Ks * spec) + (lp.Kd * dif);
        
}


int main()
{

   cout << "ka " << lp.Ka << endl;
   vtkImageData *image = NewImage(1000, 1000);
   unsigned char *buffer = 
     (unsigned char *) image->GetScalarPointer(0,0,0);
   int npixels = 1000*1000;
   for (int i = 0 ; i < npixels*3 ; i++)
       buffer[i] = 0;

   std::vector<Triangle> triangles = GetTriangles();
   int size = triangles.size();
   cout << "size of triangles vector " << size << endl;
   Screen screen;
   screen.buffer = buffer;
   screen.width = 1000;
   screen.height = 1000;
   double y1, m1, b1, m2, b2, midLine;
   double rowMin, rowMax, leftEnd, rightEnd;
   int triangleCase;
   double (*colors)[3];
   double bottomLeft[3], bottomRight[3], topVert[3];
   double BL[3], BR[3], TV[3];
   double BL_N[3], BR_N[3], TV_N[3];
   int br, bl, tv;
   for(int i =0 ; i <size ; i++){
       std::cout.setstate(std::ios_base::failbit); 
       cout << "TRIANGLE: " << i << endl;
       cout << "v0 (x,y,z) : ("<<std::setprecision(16) << triangles[i].X[0] << ", " << std::setprecision(16) << triangles[i].Y[0] << ", " << triangles[i].Z[0] << ")" <<endl;
       cout << "v1 (x,y,z) : (" << triangles[i].X[1] << ", "<< std::setprecision(16)  << triangles[i].Y[1] << ", " << triangles[i].Z[1] << ")" <<endl;
       cout << "v2 (x,y,z) : (" << triangles[i].X[2] << ", "  << triangles[i].Y[2] << ", " << triangles[i].Z[2] << ")" <<endl;
       if(triangles[i].Y[0] == triangles[i].Y[1]){
           cout << " CASE 3!!" << endl;
           if(triangles[i].X[0] > triangles[i].X[1] && triangles[i].X[2] > triangles[i].X[1]){
               tv = 0;
               br = 2;
               bl = 1;
           }
           else if(triangles[i].X[0] > triangles[i].X[1] && triangles[i].X[2] <= triangles[i].X[1]){
               tv = 1;
               br = 0;
               bl = 2;
           }
           else if(triangles[i].X[0] < triangles[i].X[1] && triangles[i].X[2] > triangles[i].X[0]){
               tv = 1;
               br = 2;
               bl = 0;
           }
           else if(triangles[i].X[0] < triangles[i].X[1] && triangles[i].X[2] <= triangles[i].X[0]){
               tv = 0;
               br = 1;
               bl = 2;
           }
           triangleCase = 3;
       }
       else if(triangles[i].Y[1] == triangles[i].Y[2]){
           cout << " CASE 4!!" << endl;
           if(triangles[i].X[1] < triangles[i].X[2] && triangles[i].X[0] > triangles[i].X[1]){
               tv = 2;
               br = 0;
               bl = 1;
           }
           else if(triangles[i].X[1] < triangles[i].X[2] && triangles[i].X[0] <= triangles[i].X[1]){
               tv = 1;
               br = 2;
               bl = 0;
           }
           if(triangles[i].X[1] > triangles[i].X[2] && triangles[i].X[0] > triangles[i].X[2]){
               tv = 1;
               br = 0;
               bl = 2;
           }
           else if(triangles[i].X[1] > triangles[i].X[2] && triangles[i].X[0] <= triangles[i].X[2]){
               tv = 2;
               br = 1;
               bl = 0;
           }
           triangleCase = 4;
       }
       else if(triangles[i].Y[2] == triangles[i].Y[0]){
           cout << " CASE 5!!" << endl;
           if(triangles[i].X[2] < triangles[i].X[0] && triangles[i].X[1] > triangles[i].X[2]){
               tv = 0;
               br = 1;
               bl = 2;
           }
           
           else if(triangles[i].X[2] < triangles[i].X[0] && triangles[i].X[1] <= triangles[i].X[2]){
               tv = 2;
               br = 0;
               bl = 1;
           }
           if(triangles[i].X[2] > triangles[i].X[0] && triangles[i].X[1] > triangles[i].X[0]){
               tv = 2;
               br = 1;
               bl = 0;
           }
           else if(triangles[i].X[2] > triangles[i].X[0] && triangles[i].X[1] <= triangles[i].X[0]){
               tv = 0;
               br = 2;
               bl = 1;
           }
           triangleCase = 5;
       }
       if(triangles[i].Y[0] > triangles[i].Y[1] && triangles[i].Y[0] > triangles[i].Y[2] ){
           tv = 0;
	   triangleCase = 0;
           cout << " CASE 0!!" << endl;
       }
       else if(triangles[i].Y[1] > triangles[i].Y[0] && triangles[i].Y[1] > triangles[i].Y[2] ){
           tv = 1;
	   triangleCase = 1;
           cout << " CASE 1!!" << endl;
       }
       else if(triangles[i].Y[2] > triangles[i].Y[0] && triangles[i].Y[2] > triangles[i].Y[1] ){
           tv = 2;
	   triangleCase = 2;
           cout << " CASE 2!!" << endl;
       }     
       //finds the left and right vertex
       if(triangleCase == 0){
           if(triangles[i].X[1] < triangles[i].X[2]){
               bl = 1;
               br = 2;
           }
           else{
               br = 1;
               bl = 2;
           }
       }
       if(triangleCase == 1){
           if(triangles[i].X[0] < triangles[i].X[2]){
               br = 2;
               bl = 0;
           }
           else{
               br = 0;
               bl = 2;
           }
       }
       if(triangleCase == 2){
           if(triangles[i].X[0] < triangles[i].X[1]){
               br = 1;
               bl = 0;
           }
           else{
               br = 0;
               bl = 1;
           }
       }
       bottomRight[0] = triangles[i].X[br];
       bottomRight[1] = triangles[i].Y[br];
       bottomRight[2] = triangles[i].Z[br];
       topVert[0] = triangles[i].X[tv];
       topVert[1] = triangles[i].Y[tv];
       topVert[2] = triangles[i].Z[tv];
       bottomLeft[0] = triangles[i].X[bl];
       bottomLeft[1] = triangles[i].Y[bl];
       bottomLeft[2] = triangles[i].Z[bl];
       TV[0] = triangles[i].colors[tv][0];
       TV[1] = triangles[i].colors[tv][1];
       TV[2] = triangles[i].colors[tv][2];
       BL[0] = triangles[i].colors[bl][0];
       BL[1] = triangles[i].colors[bl][1];
       BL[2] = triangles[i].colors[bl][2];
       BR[0] = triangles[i].colors[br][0];
       BR[1] = triangles[i].colors[br][1];
       BR[2] = triangles[i].colors[br][2];
       TV_N[0] = triangles[i].normals[tv][0];
       TV_N[1] = triangles[i].normals[tv][1];
       TV_N[2] = triangles[i].normals[tv][2];
       BL_N[0] = triangles[i].normals[bl][0];
       BL_N[1] = triangles[i].normals[bl][1];
       BL_N[2] = triangles[i].normals[bl][2];
       BR_N[0] = triangles[i].normals[br][0];
       BR_N[1] = triangles[i].normals[br][1];
       BR_N[2] = triangles[i].normals[br][2];
       cout << "topVert x,y : (" << topVert[0] << ", "  << topVert[1] << ")" <<endl;
       cout << "bottomLeft x,y : (" << bottomLeft[0] << ", "  << bottomLeft[1] << ")" <<endl;
       cout << "bottomRight x,y : (" << bottomRight[0] << ", "  << bottomRight[1] << ")" <<endl;

       cout << "TV: {" << TV[0]<< " , " << TV[1] << " , " <<  TV[2] << " )" << endl;
       cout << "BL: {" << BL[0]<< " , " << BL[1] << " , " <<  BL[2] << " )" << endl;
       cout << "BR: {" << BR[0]<< " , " << BR[1] << " , " <<  BR[2] << " )" << endl;
      

//find the middle vertex to split the triangle in half
       int side;
       if(bottomLeft[1] > bottomRight[1]){
           midLine = floor441(bottomLeft[1]);
           rowMin = ceil441(bottomRight[1]);
           //left side is higher
           side = 0;
       }else{
           rowMin = ceil441(bottomLeft[1]);
           midLine = floor441(bottomRight[1]);
           //right side is higher
           side = 1;
       }
       rowMax = floor441(topVert[1]);
       int c;
       //bottom half of triangle
       double left_shade[3], right_shade[3], _m1, _m2, _b1, _b2, _leftEnd, _rightEnd, _leftColor[3], _rightColor[3], _leftDepth, _rightDepth;
       cout << "rowMin : " << rowMin <<endl;
       cout << "midline: " << midLine << endl;
       cout << "rowMax : " << rowMax <<endl;
       for(double r1 = rowMin ; r1 <= midLine; r1++){
           cout << "first loop" << endl;
             if(side == 0){

               if(bottomLeft[0] != bottomRight[0]){
                   _m1 = (bottomLeft[1] - bottomRight[1]) / (bottomLeft[0] - bottomRight[0]);
                   _b1 = bottomRight[1] - (_m1 * bottomRight[0]);
                   _leftEnd = ceil441(((r1 - _b1) / _m1));
               }else{
                   _leftEnd = bottomLeft[0];
               }
               if(topVert[0] != bottomRight[0]){
                   _m2 = (topVert[1] - bottomRight[1]) / (topVert[0] - bottomRight[0]);
                   _b2 = bottomRight[1] - (_m2 * bottomRight[0]);
                   _rightEnd = floor441(((r1 - _b2) / _m2));
               }else{
                   _rightEnd = bottomRight[0];
               }
	       
               _leftColor[0] = interp(BR[0], BL[0], bottomRight[1], bottomLeft[1], r1);
               _leftColor[1] = interp(BR[1], BL[1], bottomRight[1], bottomLeft[1], r1);
               _leftColor[2] = interp(BR[2], BL[2], bottomRight[1], bottomLeft[1], r1);
               _rightColor[0] = interp(TV[0], BR[0], topVert[1], bottomRight[1], r1);
               _rightColor[1] = interp(TV[1], BR[1], topVert[1], bottomRight[1], r1);
               _rightColor[2] = interp(TV[2], BR[2], topVert[1], bottomRight[1], r1);
               _leftDepth = interp(bottomRight[2], bottomLeft[2], bottomRight[1], bottomLeft[1], r1);
               _rightDepth = interp(topVert[2], bottomRight[2], topVert[1], bottomRight[1], r1);
               left_shade[0] = interp(BR_N[0], BL_N[0], bottomRight[1], bottomLeft[1], r1);
               left_shade[1] = interp(BR_N[1], BL_N[1], bottomRight[1], bottomLeft[1], r1);
               left_shade[2] = interp(BR_N[2], BL_N[2], bottomRight[1], bottomLeft[1], r1);
               right_shade[0] = interp(TV_N[0], BR_N[0], topVert[1], bottomRight[1], r1);
               right_shade[1] = interp(TV_N[1], BR_N[1], topVert[1], bottomRight[1], r1);
               right_shade[2] = interp(TV_N[2], BR_N[2], topVert[1], bottomRight[1], r1);

           }else if(side ==1){
               if(topVert[0] != bottomLeft[0]){
	           _m1 = (topVert[1] - bottomLeft[1]) / (topVert[0] - bottomLeft[0]);
                   _b1 = bottomLeft[1] - (_m1 * bottomLeft[0]);
                   _leftEnd = ceil441(((r1 - _b1) / _m1));
               }else{
                   _leftEnd = topVert[0];
               }
	       if(bottomLeft[0] != bottomRight[0]){
                   _m2 = (bottomLeft[1] - bottomRight[1]) / (bottomLeft[0] - bottomRight[0]);
                   _b2 = bottomRight[1] - (_m2 * bottomRight[0]);
                   _rightEnd = floor441(((r1 - _b2) / _m2));
               }else{
                   _rightEnd = bottomRight[0];
               }
               _leftColor[0] = interp(TV[0], BL[0], topVert[1], bottomLeft[1], r1);               
               _leftColor[1] = interp(TV[1], BL[1], topVert[1], bottomLeft[1], r1);               
               _leftColor[2] = interp(TV[2], BL[2], topVert[1], bottomLeft[1], r1);               
               _rightColor[0] = interp(BR[0], BL[0], bottomRight[1], bottomLeft[1], r1);
               _rightColor[1] = interp(BR[1], BL[1], bottomRight[1], bottomLeft[1], r1);
               _rightColor[2] = interp(BR[2], BL[2], bottomRight[1], bottomLeft[1], r1);
               _leftDepth = interp(topVert[2], bottomLeft[2], topVert[1], bottomLeft[1], r1);
               _rightDepth = interp(bottomRight[2], bottomLeft[2], bottomRight[1], bottomLeft[1], r1);
               left_shade[0] = interp(TV_N[0], BL_N[0], topVert[1], bottomLeft[1], r1);              
               left_shade[1] = interp(TV_N[1], BL_N[1], topVert[1], bottomLeft[1], r1);              
               left_shade[2] = interp(TV_N[2], BL_N[2], topVert[1], bottomLeft[1], r1);               
               right_shade[0] = interp(BR_N[0], BL_N[0], bottomRight[1], bottomLeft[1], r1);
               right_shade[1] = interp(BR_N[1], BL_N[1], bottomRight[1], bottomLeft[1], r1);
               right_shade[2] = interp(BR_N[2], BL_N[2], bottomRight[1], bottomLeft[1], r1);
               
               }

           _rightEnd = floor441(_rightEnd);
           _leftEnd = floor441(_leftEnd);
           cout << "r1 : " << r1 << endl;
           cout << "case : "  << side << endl;
           cout << "left end : " << _leftEnd << endl;
           cout << "lright end : " << _rightEnd << endl;

	
           for(c = _leftEnd; c <= _rightEnd; c++){
               cout << "second loop" << endl;
               if(c >= 0 && r1 >= 0 && r1 < 1000 && c < 1000){
                   unsigned char *buf = (unsigned char *) image->GetScalarPointer(c,r1,0);   

                   double shading[3];
		   double shading_factor;
                   shading[0] = interp(left_shade[0], right_shade[0], _leftEnd, _rightEnd, c);
                   shading[1] = interp(left_shade[1], right_shade[1], _leftEnd, _rightEnd, c);
                   shading[2] = interp(left_shade[2], right_shade[2], _leftEnd, _rightEnd, c);
                   shading_factor = applyShading(shading[0], shading[1], shading[2]);
                   buf[0] = (unsigned char) ceil441(255.0 * std::min(1.0, shading_factor * interp(_leftColor[0], _rightColor[0], _leftEnd, _rightEnd, c)));
                   buf[1] = (unsigned char) ceil441(255.0 * std::min(1.0, shading_factor * interp(_leftColor[1], _rightColor[1], _leftEnd, _rightEnd, c)));
                   buf[2] = (unsigned char) ceil441(255.0 * std::min(1.0, shading_factor *  interp(_leftColor[2], _rightColor[2], _leftEnd, _rightEnd, c)));
                   double z = interp(_leftDepth, _rightDepth , _leftEnd, _rightEnd, c);
              //     Assign(r1,c,buf,z);
                 }
           } 
       }
       for(int r1 = midLine+1; r1 <= rowMax; r1++){
               cout << "third loop" << endl;
	       m1 = (topVert[1] - bottomLeft[1]) / (topVert[0] - bottomLeft[0]);
               b1 = bottomLeft[1] - (m1 * bottomLeft[0]);
               leftEnd = ceil441(((r1-b1) / m1));
               if(topVert[0] == bottomLeft[0]){
	            leftEnd = floor441(bottomLeft[0]);
               }
               _leftColor[0] = interp(BL[0],TV[0], bottomLeft[1], topVert[1], r1);
               _leftColor[1] = interp(BL[1],TV[1], bottomLeft[1], topVert[1], r1);
               _leftColor[2] = interp(BL[2],TV[2], bottomLeft[1], topVert[1], r1);

               m2 = (topVert[1] - bottomRight[1])/(topVert[0] - bottomRight[0]);
               b2 = bottomRight[1] - (m2 * bottomRight[0]);
	       rightEnd = floor441(((r1-b2) /m2));
               if(topVert[0] == bottomRight[0]){
	            rightEnd = floor441(bottomRight[0]);
               }
               _rightColor[0] = interp(BR[0],TV[0], bottomRight[1], topVert[1], r1);
               _rightColor[1] = interp(BR[1],TV[1], bottomRight[1], topVert[1], r1);
               _rightColor[2] = interp(BR[2],TV[2], bottomRight[1], topVert[1], r1);
               

               left_shade[0] = interp(TV_N[0], BL_N[0], topVert[1], bottomLeft[1], r1);              
               left_shade[1] = interp(TV_N[1], BL_N[1], topVert[1], bottomLeft[1], r1);              
               left_shade[2] = interp(TV_N[2], BL_N[2], topVert[1], bottomLeft[1], r1);              
               right_shade[0] = interp(TV_N[0], BR_N[0], topVert[1], bottomRight[1], r1);
               right_shade[1] = interp(TV_N[1], BR_N[1], topVert[1], bottomRight[1], r1);
               right_shade[2] = interp(TV_N[2], BR_N[2], topVert[1], bottomRight[1], r1);
               
                _leftDepth = interp(topVert[2], bottomLeft[2], topVert[1], bottomLeft[1], r1);
               _rightDepth = interp(topVert[2], bottomRight[2], topVert[1], bottomRight[1], r1);
           cout << "rightEnd: " <<  rightEnd << endl;
           cout << "rleftEnd: " << leftEnd << endl;
           for(c = leftEnd; c <= rightEnd; c++){
               if(c >= 0 && r1 >= 0 && r1 < 1000 && c < 1000){
                   unsigned char *buf = (unsigned char *) image->GetScalarPointer(c,r1,0);
                   double shading[3];
		   double shading_factor;
                   shading[0] = interp(left_shade[0], right_shade[0], leftEnd, rightEnd, c);
                   shading[1] = interp(left_shade[1], right_shade[1], leftEnd, rightEnd, c);
                   shading[2] = interp(left_shade[2], right_shade[2], leftEnd, rightEnd, c);
                   shading_factor = applyShading(shading[0], shading[1], shading[2]);
                   buf[0] = (unsigned char) ceil441(255.0 * std::min(1.0, shading_factor * interp( _leftColor[0] , _rightColor[0], leftEnd, rightEnd, c))); 
                   buf[1] = (unsigned char) ceil441(255.0 * std::min(1.0, shading_factor * interp( _leftColor[1] , _rightColor[1], leftEnd, rightEnd, c))); 
                   buf[2] = (unsigned char) ceil441(255.0 * std::min(1.0, shading_factor * interp( _leftColor[2] , _rightColor[2], leftEnd, rightEnd, c))); 
                   double z = interp(_leftDepth, _rightDepth , leftEnd, _rightEnd, c);
               //    Assign(r1,c,buf,z);

               }
           }
       }
}
        

   // YOUR CODE GOES HERE TO DEPOSIT TRIANGLES INTO PIXELS USING THE SCANLINE ALGORITHM

   WriteImage(image, "duck");
}
