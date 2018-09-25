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
      unsigned char color[3];

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
    rdr->SetFileName("proj1c_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();
    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkFloatArray *colors = (vtkFloatArray *) pd->GetPointData()->GetArray("color_nodal");
    float *color_ptr = colors->GetPointer(0);
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
        tris[idx].color[0] = (unsigned char) color_ptr[4*ptIds[0]+0];
        tris[idx].color[1] = (unsigned char) color_ptr[4*ptIds[0]+1];
        tris[idx].color[2] = (unsigned char) color_ptr[4*ptIds[0]+2];
    }
    cerr << "Done reading" << endl;

    return tris;

   /*std::vector<Triangle> rv(100);

   unsigned char colors[6][3] = { {255,128,0}, {255, 0, 127}, {0,204,204}, 
                                  {76,153,0}, {255, 204, 204}, {204, 204, 0}};
   for (int i = 0 ; i < 100 ; i++)
   {
       int idxI = i%10;
       int posI = idxI*100;
       int idxJ = i/10;
       int posJ = idxJ*100;
       int firstPt = (i%3);
       rv[i].X[firstPt] = posI;
       if (i == 50)
           rv[i].X[firstPt] = -10;
       rv[i].Y[firstPt] = posJ;
       rv[i].X[(firstPt+1)%3] = posI+99;
       rv[i].Y[(firstPt+1)%3] = posJ;
       rv[i].X[(firstPt+2)%3] = posI+i;
       rv[i].Y[(firstPt+2)%3] = posJ+10*(idxJ+1);
       if (i == 95)
          rv[i].Y[(firstPt+2)%3] = 1050;
       rv[i].color[0] = colors[i%6][0];
       rv[i].color[1] = colors[i%6][1];
       rv[i].color[2] = colors[i%6][2];
   }

   return rv;*/
}

//set vertexs



int main()
{
   vtkImageData *image = NewImage(1786, 1344);
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
   double bottomLeft[2], bottomRight[2], topVert[2];
   for(int i =0 /*2098251*/; i <size/* 2098252*/; i++){
       std::cout.setstate(std::ios_base::failbit); 
       cout << "TRIANGLE: " << i << endl;
       cout << "v1 (x,y) : ("<<std::setprecision(16) << triangles[i].X[0] << ", " << std::setprecision(16) << triangles[i].Y[0] << ")" <<endl;
       cout << "v2 (x,y) : (" << triangles[i].X[1] << ", "<< std::setprecision(16)  << triangles[i].Y[1] << ")" <<endl;
       cout << "v3 (x,y) : (" << triangles[i].X[2] << ", "  << triangles[i].Y[2] << ")" <<endl;

       //cout << i << endl;
//finds the top vertext

       if(triangles[i].Y[0] > triangles[i].Y[1] && triangles[i].Y[0] > triangles[i].Y[2] ){
           topVert[0] = triangles[i].X[0];
	   topVert[1] = triangles[i].Y[0];
	   triangleCase = 0;
       }
       else if(triangles[i].Y[1] > triangles[i].Y[0] && triangles[i].Y[1] > triangles[i].Y[2] ){
           topVert[0] = triangles[i].X[1];
	   topVert[1] = triangles[i].Y[1];
	   triangleCase = 1;
       }
       else if(triangles[i].Y[2] > triangles[i].Y[0] && triangles[i].Y[2] > triangles[i].Y[1] ){
           topVert[0] = triangles[i].X[2];
	   topVert[1] = triangles[i].Y[2];
	   triangleCase = 2;
       }
       else if(triangles[i].Y[0] == triangles[i].Y[1]){
           bottomLeft[0] = triangles[i].X[2];
           bottomLeft[1] = triangles[i].Y[2];
           if(triangles[i].X[0] < triangles[i].X[1]){
               bottomRight[0] = triangles[i].X[0];
               bottomRight[1] = triangles[i].Y[0];
               topVert[0] = triangles[i].X[1];
               topVert[1] = triangles[i].Y[1];
           }
           else{
               bottomRight[0] = triangles[i].X[1];
               bottomRight[1] = triangles[i].Y[1];
               topVert[0] = triangles[i].X[0];
               topVert[1] = triangles[i].Y[0];
           }
           triangleCase = 3;
       }
       else if(triangles[i].Y[1] == triangles[i].Y[2]){
           bottomLeft[0] = triangles[i].X[0];
           bottomLeft[1] = triangles[i].Y[0];
           if(triangles[i].X[1] < triangles[i].X[2]){
               bottomRight[0] = triangles[i].X[2];
               bottomRight[1] = triangles[i].Y[2];
               topVert[0] = triangles[i].X[1];
               topVert[1] = triangles[i].Y[1];
           }
           else{
               bottomRight[0] = triangles[i].X[1];
               bottomRight[1] = triangles[i].Y[1];
               topVert[0] = triangles[i].X[2];
               topVert[1] = triangles[i].Y[2];
           }
           triangleCase = 4;
       }
       else if(triangles[i].Y[0] == triangles[i].Y[2]){
           bottomLeft[0] = triangles[i].X[1];
           bottomLeft[1] = triangles[i].Y[1];
           if(triangles[i].X[0] < triangles[i].X[2]){
               bottomRight[0] = triangles[i].X[2];
               bottomRight[1] = triangles[i].Y[2];
               topVert[0] = triangles[i].X[0];
               topVert[1] = triangles[i].Y[0];
           }
           else{
               bottomRight[0] = triangles[i].X[0];
               bottomRight[1] = triangles[i].Y[0];
               topVert[0] = triangles[i].X[2];
               topVert[1] = triangles[i].Y[2];
           }
           triangleCase = 5;
        }
//finds the left and right vertex
       if(triangleCase == 0){
           if(triangles[i].X[1] < triangles[i].X[2]){
               bottomLeft[0] = triangles[i].X[1];
               bottomLeft[1] = triangles[i].Y[1];
               bottomRight[0] = triangles[i].X[2];
               bottomRight[1] = triangles[i].Y[2];
           }
           else{
           
               bottomRight[0] = triangles[i].X[1];
               bottomRight[1] = triangles[i].Y[1];
               bottomLeft[0] = triangles[i].X[2];
               bottomLeft[1] = triangles[i].Y[2];
           }
       }
       if(triangleCase == 1){
           if(triangles[i].X[0] < triangles[i].X[2]){
               bottomLeft[0] = triangles[i].X[0];
               bottomLeft[1] = triangles[i].Y[0];
               bottomRight[0] = triangles[i].X[2];
               bottomRight[1] = triangles[i].Y[2];
           }
           else{
           
               bottomRight[0] = triangles[i].X[0];
               bottomRight[1] = triangles[i].Y[0];
               bottomLeft[0] = triangles[i].X[2];
               bottomLeft[1] = triangles[i].Y[2];
           }
       }
       if(triangleCase == 2){
           if(triangles[i].X[0] < triangles[i].X[1]){
               bottomLeft[0] = triangles[i].X[0];
               bottomLeft[1] = triangles[i].Y[0];
               bottomRight[0] = triangles[i].X[1];
               bottomRight[1] = triangles[i].Y[1];
           }
           else{
           
               bottomRight[0] = triangles[i].X[0];
               bottomRight[1] = triangles[i].Y[0];
               bottomLeft[0] = triangles[i].X[1];
               bottomLeft[1] = triangles[i].Y[1];
           }
       }
       cout << "topVert x,y : (" << topVert[0] << ", "  << topVert[1] << ")" <<endl;
       cout << "bottomLeft x,y : (" << bottomLeft[0] << ", "  << bottomLeft[1] << ")" <<endl;
       cout << "bottomRight x,y : (" << bottomRight[0] << ", "  << bottomRight[1] << ")" <<endl;
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
	       cout<< "m1: " << bottomLeft[1] << " - " << bottomRight[1] << "  /  " << bottomLeft[0] << " - " << bottomRight[0] << endl;
               _b1 = bottomRight[1] - (_m1 * bottomRight[0]);
               cout << "b1 : " << bottomRight[1] << " - " <<  _m1 << " * " <<  bottomRight[0] << endl;
               _leftEnd = ceil441(((r1 - _b1) / _m1));
               cout << "leftEnd: " <<  r1 << " - " <<  _b1 << " / " <<  _m1 << endl;
           
               _m2 = (topVert[1] - bottomRight[1]) / (topVert[0] - bottomRight[0]);
               _b2 = bottomRight[1] - (_m2 * bottomRight[0]);
               _rightEnd = floor441(((r1 - _b2) / _m2));
           
           }else if(side ==1){
	       _m1 = (topVert[1] - bottomLeft[1]) / (topVert[0] - bottomLeft[0]);
	       cout<< "m1: " << topVert[1] << " - " << bottomLeft[1] << "  /  " << topVert[0] << " - " << bottomLeft[0] << endl;
               _b1 = bottomLeft[1] - (_m1 * bottomLeft[0]);
               cout << "b1 : " << bottomLeft[1] << " - " <<  _m1 << " * " <<  bottomLeft[0] << endl;
               _leftEnd = ceil441(((r1 - _b1) / _m1));
               cout << "leftEnd: " <<  r1 << " - " <<  _b1 << " / " <<  _m1 << endl;
               
               _m2 = (bottomLeft[1] - bottomRight[1]) / (bottomLeft[0] - bottomRight[0]);
               _b2 = bottomRight[1] - (_m2 * bottomRight[0]);
               _rightEnd = floor441(((r1 - _b2) / _m2));
	       cout<< "m2: " << bottomLeft[1] << " - " << bottomRight[1] << "  /  " << bottomLeft[0] << " - " << bottomRight[0] << endl;
               cout << "b2 : " << bottomRight[1] << " - " <<  _m2 << " * " <<  bottomRight[0] << endl;
               cout << "rightEnd: " <<  r1 << " - " <<  _b2 << " / " <<  _m2 << endl;
           }
	   else{
               cout << "not case one or 2" << endl;
           }
           cout << "r1 : " << r1 << endl;
           cout << "case : "  << side << endl;
           cout << "left end : " << _leftEnd << endl;
           cout << "lright end : " << _rightEnd << endl;
           for(c = _leftEnd; c <= _rightEnd; c++){
                       cout << "inside this loop" <<endl;
               if(c >= 0 && r1 >= 0 && r1 < 1344 && c < 1786){
                   unsigned char *buf = (unsigned char *) image->GetScalarPointer(c,r1,0);   
                   buf[0] = triangles[i].color[0];
                   buf[1] = triangles[i].color[1];
                   buf[2] = triangles[i].color[2];
                   cout << "PIXEL DEPOSITED 1, c  :" << c << ", r : "  << r1 << endl;
               }
           } 
       }
       for(int r1 = midLine+1; r1 <= rowMax; r1++){
           if(side == 0 || side == 1){
	       m1 = (topVert[1] - bottomLeft[1]) / (topVert[0] - bottomLeft[0]);
	       cout<< "m1: " << topVert[1] << " - " << bottomLeft[1] << "  /  " << topVert[0] << " - " << bottomLeft[0] << endl;
               b1 = bottomLeft[1] - (m1 * bottomLeft[0]);
               cout << "b1 : " << bottomLeft[1] << " - " <<  m1 << " * " <<  bottomLeft[0] << endl;
               leftEnd = ceil441(((r1-b1) / m1));
               cout << "leftEnd: " <<  r1 << " - " <<  b1 << " / " <<  m1 << endl;
               
               m2 = (topVert[1] - bottomRight[1])/(topVert[0] - bottomRight[0]);
               b2 = bottomRight[1] - (m2 * bottomRight[0]);
	       rightEnd = floor441(((r1-b2) /m2));
           }

           cout << "r1 : " << r1 << endl;
           cout << "left end : " << leftEnd << endl;
           cout << "lright end : " << rightEnd << endl;
           for(c = leftEnd; c <= rightEnd; c++){
               if(c >= 0 && r1 >= 0 && r1 < 1344 && c < 1786){
                   unsigned char *buf = (unsigned char *) image->GetScalarPointer(c,r1,0);   
                   buf[0] = triangles[i].color[0];
                   buf[1] = triangles[i].color[1];
                   buf[2] = triangles[i].color[2];

                   cout << "PIXEL DEPOSITED 2, :" << c << ", "  << r1 << endl;
                   if( c == 1069 && r1 == 1270){
                   std::cout.clear();
                   cout << "PIXEL DEPOSITED 2, :" << c << ", "  << r1 << endl << "triangle :" << i << endl;
                   std::cout.setstate(std::ios_base::failbit);
                   }
               }
           }
       }
}
        

   // YOUR CODE GOES HERE TO DEPOSIT TRIANGLES INTO PIXELS USING THE SCANLINE ALGORITHM

   WriteImage(image, "duck");
}
