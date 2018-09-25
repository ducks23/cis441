#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>

using std::cerr;
using std::endl;

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

int convert(int x){
   if(x == 0){
       return 0;
   }
   else if(x == 1){
       return 128;
   }
   else{ 
       return 255;
   }

}

int main()
{
   std::cerr << "In main!" << endl;
   vtkImageData *image = NewImage(1024, 1350);
   unsigned char *buffer = 
     (unsigned char *) image->GetScalarPointer(0,0,0);
   int i, j;

   int* dims = image->GetDimensions();
   cout << "Dims: " << "x: " << dims[0] << "y: " << dims[1] << "z: " << dims[2] << endl;
   int r,g,b, count;
   count = -1;
   for(int y = 0; y < dims[1]; y++){
       if((y % 50) == 0){
           count++;
           }
       for(int x = 0; x < dims[0]; x++){ 
           unsigned char *buffer = 
               (unsigned char *) image->GetScalarPointer(x,y,0);
           r = count / 9;
	   g = (count/3) % 3;
	   b = count % 3;  
	       //cout << "Red: " << r << " green: " << g << " blue: " << b << endl;
	   buffer[0] = convert(r);
           buffer[1] = convert(g);
	   buffer[2] = convert(b);
       }
   }
    
   WriteImage(image, "oneRedPixel");
}
