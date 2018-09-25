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
      unsigned char colors[3][3];

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
    rdr->SetFileName("proj1d_geometry.vtk");
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
    vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    float *color_ptr = var->GetPointer(0);
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
        tris[idx].X[0] = pts->GetPoint(ptIds[0])[0];
        tris[idx].X[1] = pts->GetPoint(ptIds[1])[0];
        tris[idx].X[2] = pts->GetPoint(ptIds[2])[0];
        tris[idx].Y[0] = pts->GetPoint(ptIds[0])[1];
        tris[idx].Y[1] = pts->GetPoint(ptIds[1])[1];
        tris[idx].Y[2] = pts->GetPoint(ptIds[2])[1];
        tris[idx].Z[0] = pts->GetPoint(ptIds[0])[2];
        tris[idx].Z[1] = pts->GetPoint(ptIds[1])[2];
        tris[idx].Z[2] = pts->GetPoint(ptIds[2])[2];
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

double interp(double f_a, double f_b, int a, int b, double x){
	double t = (x-a)/(b-a);
        return f_a + (t * (f_b - f_a));
}

int main()
{
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
   unsigned char (*colors)[3];
   double bottomLeft[3], bottomRight[3], topVert[3];
   unsigned char *BL, *BR, *TV;
   for(int i =0 /*2098251*/; i < 2/* 2098252*/; i++){
       //std::cout.setstate(std::ios_base::failbit); 
       cout << "TRIANGLE: " << i << endl;
       cout << "v1 (x,y,z) : ("<<std::setprecision(16) << triangles[i].X[0] << ", " << std::setprecision(16) << triangles[i].Y[0] << ", " << triangles[i].Z[0] << ")" <<endl;
       cout << "v2 (x,y,z) : (" << triangles[i].X[1] << ", "<< std::setprecision(16)  << triangles[i].Y[1] << ", " << triangles[i].Z[1] << ")" <<endl;
       cout << "v3 (x,y,z) : (" << triangles[i].X[2] << ", "  << triangles[i].Y[2] << ", " << triangles[i].Z[2] << ")" <<endl;

       //cout << i << endl;
//finds the top vertext
       colors = triangles[i].colors;
       if(triangles[i].Y[0] > triangles[i].Y[1] && triangles[i].Y[0] > triangles[i].Y[2] ){
           topVert[0] = triangles[i].X[0];
	   topVert[1] = triangles[i].Y[0];
	   topVert[2] = triangles[i].Z[0];
           TV = colors[0];
	   triangleCase = 0;
       }
       else if(triangles[i].Y[1] > triangles[i].Y[0] && triangles[i].Y[1] > triangles[i].Y[2] ){
           topVert[0] = triangles[i].X[1];
	   topVert[1] = triangles[i].Y[1];
	   topVert[2] = triangles[i].Z[1];
           TV = colors[1];
	   triangleCase = 1;
       }
       else if(triangles[i].Y[2] > triangles[i].Y[0] && triangles[i].Y[2] > triangles[i].Y[1] ){
           topVert[0] = triangles[i].X[2];
	   topVert[1] = triangles[i].Y[2];
	   topVert[2] = triangles[i].Z[2];
           TV = colors[2];
	   triangleCase = 2;
       }
       else if(triangles[i].Y[0] == triangles[i].Y[1]){
           bottomLeft[0] = triangles[i].X[2];
           bottomLeft[1] = triangles[i].Y[2];
           BL = colors[2];
           if(triangles[i].X[0] < triangles[i].X[1]){
               bottomRight[0] = triangles[i].X[0];
               bottomRight[1] = triangles[i].Y[0];
               bottomRight[2] = triangles[i].Z[0];
               BR = colors[0];
               topVert[0] = triangles[i].X[1];
               topVert[1] = triangles[i].Y[1];
               topVert[2] = triangles[i].Z[1];
               TV = colors[1];
           }
           else{
               bottomRight[0] = triangles[i].X[1];
               bottomRight[1] = triangles[i].Y[1];
               bottomRight[2] = triangles[i].Z[1];
               topVert[0] = triangles[i].X[0];
               topVert[1] = triangles[i].Y[0];
               topVert[2] = triangles[i].Z[0];
               BR = colors[1];
               TV = colors[0];
           }
           triangleCase = 3;
       }
       else if(triangles[i].Y[1] == triangles[i].Y[2]){
           bottomLeft[0] = triangles[i].X[0];
           bottomLeft[1] = triangles[i].Y[0];
           bottomLeft[2] = triangles[i].Z[0];
           BL = colors[0];
           if(triangles[i].X[1] < triangles[i].X[2]){
               bottomRight[0] = triangles[i].X[2];
               bottomRight[1] = triangles[i].Y[2];
               bottomRight[2] = triangles[i].Z[2];
               topVert[0] = triangles[i].X[1];
               topVert[1] = triangles[i].Y[1];
               topVert[2] = triangles[i].Z[1];
               BR = colors[2];
               TV = colors[1];
           }
           else{
               bottomRight[0] = triangles[i].X[1];
               bottomRight[1] = triangles[i].Y[1];
               bottomRight[2] = triangles[i].Z[1];
               topVert[0] = triangles[i].X[2];
               topVert[1] = triangles[i].Y[2];
               topVert[2] = triangles[i].Z[2];
               BR = colors[1];
               TV = colors[2];
           }
           triangleCase = 4;
       }
       else if(triangles[i].Y[0] == triangles[i].Y[2]){
           bottomLeft[0] = triangles[i].X[1];
           bottomLeft[1] = triangles[i].Y[1];
           bottomLeft[2] = triangles[i].Z[1];
           BL = colors[1];
           if(triangles[i].X[0] < triangles[i].X[2]){
               bottomRight[0] = triangles[i].X[2];
               bottomRight[1] = triangles[i].Y[2];
               bottomRight[2] = triangles[i].Z[2];
               topVert[0] = triangles[i].X[0];
               topVert[1] = triangles[i].Y[0];
               topVert[2] = triangles[i].Z[0];
               BR = colors[2];
               TV = colors[0];
           }
           else{
               bottomRight[0] = triangles[i].X[0];
               bottomRight[1] = triangles[i].Y[0];
               bottomRight[2] = triangles[i].Z[0];
               topVert[0] = triangles[i].X[2];
               topVert[1] = triangles[i].Y[2];
               topVert[2] = triangles[i].Z[2];
               TV = colors[2];
               BR = colors[0];
           }
           triangleCase = 5;
        }


//finds the left and right vertex
       if(triangleCase == 0){
           if(triangles[i].X[1] < triangles[i].X[2]){
               bottomLeft[0] = triangles[i].X[1];
               bottomLeft[1] = triangles[i].Y[1];
               bottomLeft[2] = triangles[i].Z[1];
               bottomRight[0] = triangles[i].X[2];
               bottomRight[1] = triangles[i].Y[2];
               bottomRight[2] = triangles[i].Z[2];
               BL = colors[1];
               BR = colors[2];
           }
           else{
           
               bottomRight[0] = triangles[i].X[1];
               bottomRight[1] = triangles[i].Y[1];
               bottomRight[2] = triangles[i].Z[1];
               bottomLeft[0] = triangles[i].X[2];
               bottomLeft[1] = triangles[i].Y[2];
               bottomLeft[2] = triangles[i].Z[2];
               BL = colors[2];
               BR = colors[1];
           }
       }
       if(triangleCase == 1){
           if(triangles[i].X[0] < triangles[i].X[2]){
               bottomLeft[0] = triangles[i].X[0];
               bottomLeft[1] = triangles[i].Y[0];
               bottomLeft[2] = triangles[i].Z[0];
               bottomRight[0] = triangles[i].X[2];
               bottomRight[1] = triangles[i].Y[2];
               bottomRight[2] = triangles[i].Z[2];
               BL = colors[0];
               BR = colors[2];
           }
           else{
           
               bottomRight[0] = triangles[i].X[0];
               bottomRight[1] = triangles[i].Y[0];
               bottomRight[2] = triangles[i].Z[0];
               bottomLeft[0] = triangles[i].X[2];
               bottomLeft[1] = triangles[i].Y[2];
               bottomLeft[2] = triangles[i].Z[2];
               BL = colors[2];
               BR = colors[0];
           }
       }
       if(triangleCase == 2){
           if(triangles[i].X[0] < triangles[i].X[1]){
               bottomLeft[0] = triangles[i].X[0];
               bottomLeft[1] = triangles[i].Y[0];
               bottomLeft[2] = triangles[i].Z[0];
               bottomRight[0] = triangles[i].X[1];
               bottomRight[1] = triangles[i].Y[1];
               bottomRight[2] = triangles[i].Z[1];
               BL = colors[0];
               BR = colors[1];
           }
           else{
           
               bottomRight[0] = triangles[i].X[0];
               bottomRight[1] = triangles[i].Y[0];
               bottomRight[2] = triangles[i].Z[0];
               bottomLeft[0] = triangles[i].X[1];
               bottomLeft[1] = triangles[i].Y[1];
               bottomLeft[2] = triangles[i].Z[1];
               BL = colors[1];
               BR = colors[0];
           }
       }
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
       if(rowMax > 1785){ rowMax = 1785;}
       //bottom half of triangle
       cout << "rowMin : " << rowMin <<endl;
       cout << "midline: " << midLine << endl;
       cout << "rowMax : " << rowMax <<endl;
       double _m1, _m2, _b1, _b2, _leftEnd, _rightEnd;
       for(double r1 = rowMin ; r1 <= midLine; r1++){
           if(side == 0){
               _m1 = (bottomLeft[1] - bottomRight[1]) / (bottomLeft[0] - bottomRight[0]);
               _b1 = bottomRight[1] - (_m1 * bottomRight[0]);
               _leftEnd = ceil441(((r1 - _b1) / _m1));
           
               _m2 = (topVert[1] - bottomRight[1]) / (topVert[0] - bottomRight[0]);
               _b2 = bottomRight[1] - (_m2 * bottomRight[0]);
               _rightEnd = floor441(((r1 - _b2) / _m2));
           
           }else if(side ==1){
	       _m1 = (topVert[1] - bottomLeft[1]) / (topVert[0] - bottomLeft[0]);
               _b1 = bottomLeft[1] - (_m1 * bottomLeft[0]);
               _leftEnd = ceil441(((r1 - _b1) / _m1));
               
               _m2 = (bottomLeft[1] - bottomRight[1]) / (bottomLeft[0] - bottomRight[0]);
               _b2 = bottomRight[1] - (_m2 * bottomRight[0]);
               _rightEnd = floor441(((r1 - _b2) / _m2));
           }
	   else{
               cout << "not case one or 2" << endl;
           }
           cout << "case : "  << side << endl;
           cout << "left end : " << _leftEnd << endl;
           cout << "lright end : " << _rightEnd << endl;
           for(c = _leftEnd; c <= _rightEnd; c++){
                       cout << "inside this loop" <<endl;
               if(c >= 0 && r1 >= 0 && r1 < 1344 && c < 1786){
                   unsigned char *buf = (unsigned char *) image->GetScalarPointer(c,r1,0);   
           //        buf[0] = triangles[i].colors[0];
             //      buf[1] = triangles[i].colors[1];
               //    buf[2] = triangles[i].colors[2];
                   cout << "PIXEL DEPOSITED 1, c  :" << c << ", r : "  << r1 << endl;
               }
           } 
       }
       for(int r1 = midLine+1; r1 <= rowMax; r1++){
           if(side == 0 || side == 1){
	       m1 = (topVert[1] - bottomLeft[1]) / (topVert[0] - bottomLeft[0]);
               b1 = bottomLeft[1] - (m1 * bottomLeft[0]);
               leftEnd = ceil441(((r1-b1) / m1));
               
               m2 = (topVert[1] - bottomRight[1])/(topVert[0] - bottomRight[0]);
               b2 = bottomRight[1] - (m2 * bottomRight[0]);
	       rightEnd = floor441(((r1-b2) /m2));
           }

           cout << "left end : " << leftEnd << endl;
           cout << "lright end : " << rightEnd << endl;
           for(c = leftEnd; c <= rightEnd; c++){
               if(c >= 0 && r1 >= 0 && r1 < 1344 && c < 1786){
                   unsigned char *buf = (unsigned char *) image->GetScalarPointer(c,r1,0);   
            //       buf[0] = triangles[i].colors[0];
            //       buf[1] = triangles[i].colors[1];
            //       buf[2] = triangles[i].colors[2];

                   cout << "PIXEL DEPOSITED 2, :" << c << ", "  << r1 << endl;
               }
           }
       }
}
        

   // YOUR CODE GOES HERE TO DEPOSIT TRIANGLES INTO PIXELS USING THE SCANLINE ALGORITHM

   WriteImage(image, "duck");
}
