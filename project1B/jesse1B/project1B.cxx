#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <math.h>
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
   std::vector<Triangle> rv(100);

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

   return rv;
}

//set vertexs



int main()
{
   vtkImageData *image = NewImage(1000, 1000);
   unsigned char *buffer = 
     (unsigned char *) image->GetScalarPointer(0,0,0);
   int npixels = 1000*1000;
   for (int i = 0 ; i < npixels*3 ; i++)
       buffer[i] = 0;

   std::vector<Triangle> triangles = GetTriangles();
   
   Screen screen;
   screen.buffer = buffer;
   screen.width = 1000;
   screen.height = 1000;
   int numT = 100;
   double y1,m1,x1,b1,y2,m2,x2,b2,m3,b3;
   double rowMin, rowMax, leftEnd, rightEnd;
   int triangleCase;
   double bottomLeft[2], bottomRight[2], topVert[2];
   for(int i = 0; i < numT; i++){
       cout << "v1 (x,y) : (" << triangles[i].X[0] << ", "  << triangles[i].Y[0] << ")" <<endl;
       cout << "v2 (x,y) : (" << triangles[i].X[1] << ", "  << triangles[i].Y[1] << ")" <<endl;
       cout << "v3 (x,y) : (" << triangles[i].X[2] << ", "  << triangles[i].Y[2] << ")" <<endl;

       cout << i << endl;
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
       cerr << "top vertex:  x: " << std::setprecision(6) << topVert[0]<< " y: "  << topVert[1] <<std::setprecision(6) <<  endl;
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
       
       rowMin =  ceil441(bottomLeft[1]);
       rowMax = floor441(topVert[1]);
       cerr << "Height of rownMin: " << std::setprecision(6) << rowMin << endl<< " rowMax: "  << rowMax << std::setprecision(6) <<  endl;
       int c;
       if(rowMax > 999){ rowMax = 999;}

       for(int r = rowMin; r <= rowMax; r++){
           if(topVert[0] != bottomLeft[0]){

	       m1 = (topVert[1] - bottomLeft[1]) / (topVert[0] - bottomLeft[0]);
               b1 = bottomLeft[1] - (m1 * bottomLeft[0]);
           
               //cerr << "r :  " << r << endl;
               //cerr << "m1: " << std::setprecision(6) << m1 << " b1: " <<std::setprecision(6) <<b1 << endl;
               leftEnd = (r-b1) / m1;
           }else{
               leftEnd = topVert[0];
           }
           
	   cerr << "left end: " << std::setprecision(6) << leftEnd << endl;
	   if(topVert[0] != bottomRight[0]){
               m2 = (topVert[1] - bottomRight[1])/(topVert[0] - bottomRight[0]);
               b2 = bottomRight[1] - (m2 * bottomRight[0]);
	       rightEnd = (r-b2) /m2;
           }else{
               rightEnd = topVert[0];
           }
           
           leftEnd = ceil441(leftEnd);
	   if(leftEnd < 0){
               leftEnd = 0;
           }
           cerr << "right end: " << std::setprecision(6) << rightEnd << endl;
	   if(isnan(rightEnd) || rightEnd < 0){ rightEnd = 999;}
	   int max = floor441(rightEnd+1);
	   cerr << "max : " << std::setprecision(6) << max << endl;
          
           for(c = leftEnd; c < max; c++){
               
	      // cerr << "R :  "<< std::setprecision(6) << r << endl;
	       //cerr << "C :  "<< std::setprecision(6) << c << endl;
	       unsigned char *buf = (unsigned char *) image->GetScalarPointer(c,r,0);
	       buf[0] = triangles[i].color[0];               
               buf[1] = triangles[i].color[1];               
               buf[2] = triangles[i].color[2];               

           } 

       }
}
        

   // YOUR CODE GOES HERE TO DEPOSIT TRIANGLES INTO PIXELS USING THE SCANLINE ALGORITHM

   WriteImage(image, "allTriangles");
}
