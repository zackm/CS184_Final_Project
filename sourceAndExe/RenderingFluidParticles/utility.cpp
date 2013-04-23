
#include "FluidSimulation.h"

inline void CrossProduct(double *a, double *b, double *c)
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

inline void Normalize(double *a)
{
    double length = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);

    if(!length) return;

    double inv_length = 1.0/length;

    a[0] *= inv_length;
    a[1] *= inv_length;
    a[2] *= inv_length;
}

inline double DotProduct(double *a, double *b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline double DotProduct4D(double *a, double *b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3];
}

inline void GetArea(double *normal, double &area)
{
    // assuming that "normal" is not normalized at this point 
    area = 0.5 * sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
}


inline double GetDistance(double *a, double *b)
{
    return sqrt( pow(a[0]-b[0], 2.0) + pow(a[1]-b[1], 2.0) + pow(a[2]-b[2], 2.0) ); 
}

inline double GetLength(double *a)
{
    return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]); 
}


inline void Swap(double &a, double &b)
{
    double t = a;  a = b;  b = t;
}

inline double Max(double a, double b)
{
    if(a > b) return a;
    else      return b;
}

inline double Min(double a, double b)
{
    if(a < b) return a;
    else      return b;
}


inline bool SolveLinearSystem(double (*matrix)[4], double *rhs, double *solution)
{  
    // perform gaussian elimination

    for(int i = 0; i <= 3; i++){
        double maxAbsoluteValue = -1.0;
        int    pivot_index;

        // choose maximum "matrix[j][i]" as a pivot
        for(int j = i; j <= 3; j++){
            if(fabs(matrix[j][i]) > maxAbsoluteValue){
                maxAbsoluteValue = fabs(matrix[j][i]);
                pivot_index      = j;
            }
        }

        // matrix is singular
        if(maxAbsoluteValue < 1.0e-6) return false;

        for(int j = i; j <= 3; j++) Swap(matrix[i][j], matrix[pivot_index][j]);
        Swap(rhs[i], rhs[pivot_index]);

        double scale = 1.0 / matrix[i][i];

        for(int j = i+1; j <= 3; j++){
            double pivot = -matrix[j][i]*scale;

            for(int k = 0;   k <= i;   k++) matrix[j][k] = 0.0;
            for(int k = i+1; k <= 3; k++){
                if(fabs(matrix[i][k]) > 1.0e-6)
                    matrix[j][k] += matrix[i][k] * pivot;
                else
                    break;
            }

            rhs[j] += rhs[i] * pivot;
        }
    }

   

//    cerr << "diag " << matrix[0][0] << " " << matrix[1][1] << " "  << matrix[2][2] << " " <<  matrix[3][3] << endl;
    

//    if(fabs(matrix[3][3]) < 1.0e-6) return false;

    for(int i = 3; i >= 0; i--){  
        solution[i] = 0.0;
        for(int j = i+1; j <= 3; j++) solution[i] += solution[j]*matrix[i][j];

        solution[i] = (rhs[i]-solution[i])/matrix[i][i];
    }


    return true;
}

inline void ComputeReflectionVector(double *normal, double *vector_in, double *reflection_vector)
{
    double scale = 2.0 * fabs( DotProduct(vector_in, normal) );

    for(int i = 0; i < 3; i++) reflection_vector[i] = vector_in[i] + scale*normal[i];
    Normalize(reflection_vector);
}

inline void ComputeRefractionVector(double *normal, double *vector_in, double *refraction_vector)
{
    double cos_angle = DotProduct(vector_in, normal);
    double angle_i, angle_r;

    if(cos_angle < 0.0){
        // from air to water
        angle_i = acos( -cos_angle );
        angle_r = asin( 0.7502*sqrt(1-cos_angle*cos_angle) );
    }else{
        // from water to air
        angle_r = acos(cos_angle);
        angle_i = asin( 1.333*sqrt(1-cos_angle*cos_angle) );
    }

    double rotation_axis[3];
    CrossProduct(normal, vector_in, rotation_axis);

    RotateAroundAxis(rotation_axis, angle_i-angle_r, vector_in, refraction_vector);

}



void RotateAroundAxis(double *rotation_axis, double theta, double *vec_in, double *vec_out)
{
    // rotate "vec_in" around "rotation_axis" using quarternion

    Normalize(rotation_axis);

    double common_factor = sin(theta*0.5);
    
    double a = cos(theta*0.5);
    double b = rotation_axis[0] * common_factor;
    double c = rotation_axis[1] * common_factor;
    double d = rotation_axis[2] * common_factor;
 
    double mat[9] = { a*a+b*b-c*c-d*d,     2*(b*c-a*d),      2*(b*d+a*c),
                          2*(b*c+a*d), a*a-b*b+c*c-d*d,      2*(c*d-a*b),
                          2*(b*d-a*c),     2*(c*d+a*b),  a*a-b*b-c*c+d*d };

    for(int i = 0; i < 3; i++) vec_out[i] = 0.0;

    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++) vec_out[i] += mat[i*3+j] * vec_in[j];
    }

}


bool ComputeRayBoxIntersection(double *org, double *dir, double *inv_dir, double *bbmin, double *bbmax, 
                               double &nearest_t, double *intersection, double *normal)
{
    nearest_t = 1.0e6;

    int nearest_object_index = -1;

    for(int i = 0; i < 6; i++){

        double t;

        switch(i){
        case 0:
            //if(org[0] > bbmin[0] || dir[0] < 0.0) continue;
            t = (bbmin[0]-org[0])*inv_dir[0];   
            break;
        case 1:
            //if(org[0] < bbmax[0] || dir[0] > 0.0) continue;            
            t = (bbmax[0]-org[0])*inv_dir[0]; 
            break;
        case 2: 
            //if(org[1] < bbmin[1] || dir[1] > 0.0) continue;     
            t = (bbmin[1]-org[1])*inv_dir[1]; 
            break;
        case 3: 
            //if(org[1] < bbm || dir[1] > 0.0) continue;     
            t = (bbmax[1]-org[1])*inv_dir[1]; 
            break;
        case 4: 
            //if(org[2] > -0.01 || dir[2] < 0.0) continue;
            t = (bbmin[2]-org[2])*inv_dir[2]; 
            break;
        case 5: 
            //if(org[2] < 0.19 || dir[2] > 0.0) continue;
            t = (bbmax[2]-org[2])*inv_dir[2];
            break;
        }

        if(t > 1.0e-6 && t < nearest_t){

            double intersect_1, intersect_2;

            switch(i){
            case 0:
            case 1:
                intersect_1 = org[1]+dir[1]*t;
                intersect_2 = org[2]+dir[2]*t;

                if(intersect_1 > bbmin[1] && intersect_1 < bbmax[1] && intersect_2 > bbmin[2] && intersect_2 < bbmax[2]){
                    nearest_t = t;
                    nearest_object_index = i;
                }
                break;
            case 2:
            case 3:
                intersect_1 = org[0]+dir[0]*t;
                intersect_2 = org[2]+dir[2]*t;

                if(intersect_1 > bbmin[0] && intersect_1 < bbmax[0] && intersect_2 > bbmin[2] && intersect_2 < bbmax[2]){
                    nearest_t = t;
                    nearest_object_index = i;
                }
                break;
            case 4:
            case 5:
                intersect_1 = org[0]+dir[0]*t;
                intersect_2 = org[1]+dir[1]*t;

                if(intersect_1 > bbmin[0] && intersect_1 < bbmax[0] && intersect_2 > bbmin[1] && intersect_2 < bbmax[1]){
                    nearest_t = t;
                    nearest_object_index = i;
                }
                break;
            }

        } // if(t > 0.0 && t < nearest_t){

    } // for(int i = 0; i < 6; i++){


    // ray does not intersect with this box
    if(nearest_object_index == -1) return false;


    for(int i = 0; i < 3; i++) intersection[i] = org[i] + dir[i]*nearest_t;

    switch(nearest_object_index){
    case 0:
        normal[0] = -1.0; normal[1] = 0.0; normal[2] = 0.0; 
        break;
    case 1:
        normal[0] =  1.0; normal[1] = 0.0; normal[2] = 0.0; 
        break;
    case 2:
        normal[0] = 0.0; normal[1] = -1.0; normal[2] = 0.0; 
        break;
    case 3:
        normal[0] = 0.0; normal[1] =  1.0; normal[2] = 0.0; 
        break;
    case 4:
        normal[0] = 0.0; normal[1] = 0.0; normal[2] = -1.0; 
        break;
    case 5:
        normal[0] = 0.0; normal[1] = 0.0; normal[2] = 1.0; 
        break;
    }

    return true;
}


double random0to1()
{
    static double invRAND_MAX = 1.0/RAND_MAX;
    return rand()*invRAND_MAX;
}
 

float* LoadPngData(WCHAR *filename, unsigned int &img_width, unsigned int &img_height)
{
    Bitmap bitmap(filename);

    img_width  = bitmap.GetWidth();
    img_height = bitmap.GetHeight();

    cerr << img_width << " " << img_height << endl;

    float *pixel_data = new float [img_width*img_height*3];

    for(unsigned int i = 0; i < img_height; i++){
        for(unsigned int j = 0; j < img_width; j++){
            Color pixelColor;
            bitmap.GetPixel(j, i, &pixelColor);

            static float scale = 1.0/255.0;
            pixel_data[((img_height-1-i)*img_width+j)*3+0] = pixelColor.GetR() * scale;
            pixel_data[((img_height-1-i)*img_width+j)*3+1] = pixelColor.GetG() * scale;
            pixel_data[((img_height-1-i)*img_width+j)*3+2] = pixelColor.GetB() * scale;
        }
    }

    return pixel_data;
}

int GetEncoderClsid(const WCHAR* format, CLSID* pClsid)
{
   UINT  num = 0;          // number of image encoders
   UINT  size = 0;         // size of the image encoder array in bytes

   ImageCodecInfo* pImageCodecInfo = NULL;

   GetImageEncodersSize(&num, &size);
   if(size == 0)
      return -1;  // Failure

   pImageCodecInfo = (ImageCodecInfo*)(malloc(size));
   if(pImageCodecInfo == NULL)
      return -1;  // Failure

   GetImageEncoders(num, size, pImageCodecInfo);

   for(UINT j = 0; j < num; ++j)
   {
      if( wcscmp(pImageCodecInfo[j].MimeType, format) == 0 )
      {
         *pClsid = pImageCodecInfo[j].Clsid;
         free(pImageCodecInfo);
         return j;  // Success
      }    
   }

   free(pImageCodecInfo);
   return -1;  // Failure
}


void SaveImage(wchar_t *filename, int w, int h, float *buffer)
{
    // assuming that GDI+ is already initialized. 
    // Otherwise, this function fails. 
    Bitmap output_img(w, h);

    for(int i = 0; i < h; i++){
        for(int j = 0; j < w; j++){
            int pixel_index = (i*w+j)*3;
            Color color(255, Min(buffer[pixel_index],   1.0)*255, 
                             Min(buffer[pixel_index+1], 1.0)*255, 
                             Min(buffer[pixel_index+2], 1.0)*255);
         
            output_img.SetPixel(j, h-1-i, color);
        }
    }

    CLSID encoderClsid;

    // Get the CLSID of the PNG encoder.
    GetEncoderClsid(L"image/png", &encoderClsid);

    output_img.Save(filename, &encoderClsid, NULL);
}

