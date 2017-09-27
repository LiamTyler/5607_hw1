#include "image.h"
// #include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <iomanip>

/**
 * Image
 **/
Image::Image (int width_, int height_){

    assert(width_ > 0);
    assert(height_ > 0);

    width           = width_;
    height          = height_;
    num_pixels      = width * height;
    sampling_method = IMAGE_SAMPLING_POINT;

    data.raw = new uint8_t[num_pixels*4];
    int b = 0; //which byte to write to
    for (int j = 0; j < height; j++){
        for (int i = 0; i < width; i++){
            data.raw[b++] = 0;
            data.raw[b++] = 0;
            data.raw[b++] = 0;
            data.raw[b++] = 0;
        }
    }

    assert(data.raw != NULL);
}

Image::Image (const Image& src){

    width           = src.width;
    height          = src.height;
    num_pixels      = width * height;
    sampling_method = IMAGE_SAMPLING_POINT;

    data.raw = new uint8_t[num_pixels*4];

    memcpy(data.raw, src.data.raw, 4*num_pixels);
    // *data.raw = *src.data.raw;
}

Image::Image (char* fname){

    int numComponents; //(e.g., Y, YA, RGB, or RGBA)
    data.raw = stbi_load(fname, &width, &height, &numComponents, 4);

    if (data.raw == NULL){
        printf("Error loading image: %s", fname);
        exit(-1);
    }


    num_pixels = width * height;
    sampling_method = IMAGE_SAMPLING_POINT;

}

Image::~Image (){
    delete data.raw;
    data.raw = NULL;
}

void Image::Write(char* fname){

    int lastc = strlen(fname);

    switch (fname[lastc-1]){
        case 'g': //jpeg (or jpg) or png
            if (fname[lastc-2] == 'p' || fname[lastc-2] == 'e') //jpeg or jpg
                stbi_write_jpg(fname, width, height, 4, data.raw, 95);  //95% jpeg quality
            else //png
                stbi_write_png(fname, width, height, 4, data.raw, width*4);
            break;
        case 'a': //tga (targa)
            stbi_write_tga(fname, width, height, 4, data.raw);
            break;
        case 'p': //bmp
        default:
            stbi_write_bmp(fname, width, height, 4, data.raw);
    }
}

void Image::AddNoise (double factor)
{
    int x,y;
    for (x = 0 ; x < Width() ; x++)
    {
        for (y = 0 ; y < Height() ; y++)
        {
            Pixel p = GetPixel(x, y);
            double rr = rand() % 256;
            double rg = rand() % 256;
            double rb = rand() % 256;
            double nr = p.r;
            double ng = p.g;
            double nb = p.b;
            nr = (1-factor)*nr + factor*rr;
            ng = (1-factor)*ng + factor*rg;
            nb = (1-factor)*nb + factor*rb;
            GetPixel(x,y) = Pixel(nr, ng, nb, p.a);
        }
    }
}

void Image::Brighten (double factor)
{
    int x,y;
    for (x = 0 ; x < Width() ; x++)
    {
        for (y = 0 ; y < Height() ; y++)
        {
            Pixel p = GetPixel(x, y);
            Pixel scaled_p = p*factor;
            GetPixel(x,y) = scaled_p;
        }
    }
}


void Image::ChangeContrast (double factor)
{
    int x,y;
    for (x = 0 ; x < Width() ; x++)
    {
        for (y = 0 ; y < Height() ; y++)
        {
            Pixel p = GetPixel(x, y);
            Component nr = ComponentClamp(128 + factor*(p.r - 128.0));
            Component ng = ComponentClamp(128 + factor*(p.g - 128.0));
            Component nb = ComponentClamp(128 + factor*(p.b - 128.0));
            GetPixel(x,y) = Pixel(nr, ng, nb, p.a);
        }
    }
}


void Image::ChangeSaturation(double factor)
{
    int x,y;
    for (x = 0 ; x < Width() ; x++)
    {
        for (y = 0 ; y < Height() ; y++)
        {
            Pixel p = GetPixel(x, y);
            Component l = p.Luminance();
            Component nr = ComponentClamp((1-factor)*l + factor*p.r);
            Component ng = ComponentClamp((1-factor)*l + factor*p.g);
            Component nb = ComponentClamp((1-factor)*l + factor*p.b);
            GetPixel(x,y) = Pixel(nr, ng, nb, p.a);
        }
    }
}


Image* Image::Crop(int x, int y, int w, int h)
{
    Image *ret = new Image(w, h);
    for (int h2 = 0; h2 < h; h2++) {
        for (int w2 = 0; w2 < w; w2++) {
            int i = width*(y+h2) + x + w2;
            if (i > width*height)
                ret->data.pixels[w*h2 + w2] = Pixel(0, 0, 0, 0);
            else
                ret->data.pixels[w*h2 + w2] = data.pixels[width*(y+h2) + x + w2];
        }
    }

    return ret;
}


void Image::ExtractChannel(int channel)
{
    for (int h = 0; h < height; h++) {
        for (int w = 0; w < width; w++) {
            Pixel p = GetPixel(w, h);
            if (channel == 0) {
                GetPixel(w, h) = Pixel(p.r, 0, 0, p.a);
            } else if (channel == 1) {
                GetPixel(w, h) = Pixel(0, p.g, 0, p.a);
            } else {
                GetPixel(w, h) = Pixel(0, 0, p.b, p.a);
            }
        }
    }
}


void Image::Quantize (int nbits)
{
    double d = 255.0 / (nbits - 1);
    double n = nbits - 1;
    for (int x = 0 ; x < Width() ; x++)
    {
        for (int y = 0 ; y < Height() ; y++)
        {
            Pixel p = GetPixel(x,y);
            int  r = 255 * round(p.r / d) / n;
            int  g = 255 * round(p.g / d) / n;
            int  b = 255 * round(p.b / d) / n;
            Pixel p2;
            p2.SetClamp(r, g, b, 255);
            GetPixel(x,y) = p2;
        }
    }
}

void Image::RandomDither (int nbits)
{
    double d = 255.0 / (nbits - 1);
    double n = nbits - 1;
    for (int x = 0 ; x < Width() ; x++)
    {
        for (int y = 0 ; y < Height() ; y++)
        {
            Pixel p = GetPixel(x,y);
            int  lr = 255 * std::floor(p.r / d) / n;
            int  hr = 255 * std::ceil(p.r / d) / n;
            int  lg = 255 * std::floor(p.r / d) / n;
            int  hg = 255 * std::ceil(p.r / d) / n;
            int  lb = 255 * std::floor(p.r / d) / n;
            int  hb = 255 * std::ceil(p.r / d) / n;
            int  r = (p.r < (lr + rand()%(hr-lr+1))) ? lr : hr;
            int  g = (p.g < (lg + rand()%(hg-lg+1))) ? lg : hg;
            int  b = (p.b < (lb + rand()%(hb-lb+1))) ? lb : hb;
            Pixel p2;
            p2.SetClamp(r, g, b, 255);
            GetPixel(x,y) = p2;
        }
    }
}


static int Bayer4[4][4] =
{
    {15,  7, 13,  5},
    { 3, 11,  1,  9},
    {12,  4, 14,  6},
    { 0,  8,  2, 10}
};


void Image::OrderedDither(int nbits)
{
    double d = 255.0 / (nbits - 1);
    double n = nbits - 1;
    for (int y = 0 ; y < Height(); y++) {
        for (int x = 0 ; x < Width(); x++) {
            int bay = Bayer4[x % 4][y % 4] / 16.0 * d;
            Pixel p = GetPixel(x,y);
            int  lr = 255 * std::floor(p.r / d) / n;
            int  hr = 255 * std::ceil(p.r / d) / n;
            int  lg = 255 * std::floor(p.r / d) / n;
            int  hg = 255 * std::ceil(p.r / d) / n;
            int  lb = 255 * std::floor(p.r / d) / n;
            int  hb = 255 * std::ceil(p.r / d) / n;
            int  r = (p.r < (lr + bay)) ? lr : hr;
            int  g = (p.g < (lg + bay)) ? lg : hg;
            int  b = (p.b < (lb + bay)) ? lb : hb;
            Pixel p2;
            p2.SetClamp(r, g, b, 255);
            GetPixel(x,y) = p2;
        }
    }
}

/* Error-diffusion parameters */
const double
ALPHA = 7.0 / 16.0,
      BETA  = 3.0 / 16.0,
      GAMMA = 5.0 / 16.0,
      DELTA = 1.0 / 16.0;


void Image::FloydSteinbergDither(int nbits)
{
    double d = 255.0 / (nbits - 1);
    double n = nbits - 1;
    for (int y = 0 ; y < Height() ; y++) {
        int er = 0, eg = 0, eb = 0;
        for (int x = 0 ; x < Width(); x++) {
            Pixel p = GetPixel(x,y);
            int  r = 255 * round((p.r+er) / d) / n;
            int  g = 255 * round((p.g+eg) / d) / n;
            int  b = 255 * round((p.b+eb) / d) / n;
            Pixel p2;
            p2.SetClamp(r, g, b, 255);
            GetPixel(x,y) = p2;
            er = p.r - p2.r;
            eg = p.g - p2.g;
            eb = p.b - p2.b;
            Pixel p3;
            if (x + 1 < width) {
                p3 = GetPixel(x+1, y);
                p3.SetClamp(p3.r + er*ALPHA, p3.g + eg*ALPHA, p3.b + eb*ALPHA, 255);
                GetPixel(x+1, y) = p3;
            }
            if (y + 1 < height) {
                if (x - 1 >= 0) {
                    p3 = GetPixel(x-1, y+1);
                    p3.SetClamp(p3.r + er*BETA, p3.g + eg*BETA, p3.b + eb*BETA, 255);
                    GetPixel(x-1, y+1) = p3;
                }
                p3 = GetPixel(x, y+1);
                p3.SetClamp(p3.r + er*GAMMA, p3.g + eg*GAMMA, p3.b + eb*GAMMA, 255);
                GetPixel(x, y+1) = p3;
                if (x + 1 < width) {
                    p3 = GetPixel(x+1, y+1);
                    p3.SetClamp(p3.r + er*DELTA, p3.g + eg*DELTA, p3.b + eb*DELTA, 255);
                    GetPixel(x+1, y+1) = p3;
                }
            }
        }
    }
}

Pixel Image::ApplyKernel(int x, int y, int size, double** kern) {
    double r = 0;
    double g = 0;
    double b = 0;
    int c = size / 2;
    for (int row = -c; row <= c; row++) {
        for (int col = -c; col <= c; col++) {
            int xx = std::min(width - 1, std::max(0, x + col));
            int yy = std::min(height - 1, std::max(0, y + row));
            Pixel p = GetPixel(xx, yy);
            r += p.r * kern[row+c][col+c];
            g += p.g * kern[row+c][col+c];
            b += p.b * kern[row+c][col+c];
        }
    }
    Pixel s;
    s.SetClamp(r,g,b,255);
    return s;
}

void Image::KernelFilter(int size, double** kern) {
    Image copy(*this);
    int x,y;
    for (x = 0; x < Width() ; x++)
    {
        for (y = 0; y < Height() ; y++)
        {
            copy.GetPixel(x,y) = ApplyKernel(x, y, size, kern);
        }
    }
    for (x = 0 ; x < Width() ; x++)
    {
        for (y = 0 ; y < Height() ; y++)
        {
            GetPixel(x,y) = copy.GetPixel(x,y);
        }
    }
}

void Image::Blur(int n)
{
    int SIZE = n;
    double MAX = 5;
    double** kern = new double*[SIZE];
    for (int i = 0; i < SIZE; i++)
        kern[i] = new double[SIZE];

    double total = 0;
    double center = SIZE / 2;
    double CAP = std::sqrt(2*center*center);
    for (int r = 0; r < SIZE; r++) {
        for (int c = 0; c < SIZE; c++) {
            double dx = -center + c;
            double dy = -center + r;
            kern[r][c] = MAX*(1 - (std::sqrt(dx*dx + dy*dy)/CAP));
            total += kern[r][c];
        }
    }
    for (int r = 0; r < SIZE; r++) {
        for (int c = 0; c < SIZE; c++) {
            kern[r][c] /= total;
            // std::cout << std::setprecision(4) << std::setw(8) <<
            //     std::setfill(' ') << kern[r][c] << " ";
        }
        // std::cout << std::endl;
    }
    KernelFilter(SIZE, kern);
}

void Image::Sharpen(int n)
{
    Image copy(*this);
    copy.Blur(n);
    for (int x = 0 ; x < Width() ; x++)
    {
        for (int y = 0 ; y < Height() ; y++)
        {
            Pixel cp = copy.GetPixel(x,y);
            Pixel op = GetPixel(x,y);
            Pixel p;
            p.SetClamp(2*op.r-cp.r, 2*op.g-cp.g, 2*op.b-cp.b, 255);
            GetPixel(x,y) = p;
        }
    }
}

void Image::EdgeDetect()
{
    int SIZE = 3;
    double** kern = new double*[SIZE];
    for (int i = 0; i < SIZE; i++)
        kern[i] = new double[SIZE];

    for (int r = 0; r < SIZE; r++) {
        for (int c = 0; c < SIZE; c++) {
            kern[r][c] = -1;
        }
    }
    kern[1][1] = 8;
    KernelFilter(SIZE, kern);
}

Image* Image::Scale(double sx, double sy)
{
    // SetSamplingMethod(IMAGE_SAMPLING_POINT);
    // SetSamplingMethod(IMAGE_SAMPLING_BILINEAR);
    int new_width = width*sx;
    int new_height = height*sy;
    Image *img = new Image(new_width, new_height);
    for (int y = 0 ; y < new_height; y++) {
        for (int x = 0 ; x < new_width; x++) {
            double u = x / (double) new_width;
            double v = y / (double) new_height;
            img->GetPixel(x,y) = Sample(u,v);
        }
    }
    return img;
}

typedef struct Point {
    Point() {
        x = 0;
        y = 0;
    }
    Point(double xx, double yy) {
        x = xx;
        y = yy;
    }
    double x;
    double y;
    friend Point operator-(Point &a, Point& b) {
        return Point(a.x - b.x, a.y - b.y);
    }
    friend Point operator*(double x, Point& a) {
        return Point(x * a.x, x * a.y);
    }
} Point;

Point rotate(double t, Point& p) {
    Point ret;
    double c = std::cos(t);
    double s = std::sin(t);
    ret.x = c*p.x - s*p.y;
    ret.y = s*p.x + c*p.y;
    return ret;
}

Image* Image::Rotate(double angle)
{
    using namespace std;
    angle = angle * M_PI / 180.0;
    Point center = Point(width / 2.0, height / 2.0);
    Point ur = center;
    Point lr = center;
    lr.y *= -1;
    ur = rotate(angle, ur);
    lr = rotate(angle, lr);
    int new_width = 2*ceil(max(fabs(ur.x), fabs(lr.x)));
    int new_height = 2*ceil(max(fabs(ur.y), fabs(lr.y)));
    Image* n = new Image(new_width, new_height);

    Point new_center = Point(new_width / 2, new_height / 2);
    for (int x = 0 ; x < new_width ; x++)
    {
        for (int y = 0 ; y < new_height ; y++)
        {
            Point p(x - new_center.x, y - new_center.y);
            p = rotate(-angle, p);
            p.x += center.x;
            p.y += center.y;
            double u = (p.x) / (double) width;
            double v = (p.y) / (double) height;
            n->GetPixel(x,y) = Sample(u,v);
        }
    }
    /*
       center = Point(new_width / 2, new_height / 2);
       for (int x = 0 ; x < new_width ; x++)
       {
       for (int y = 0 ; y < new_height ; y++)
       {
       Point p(x, -y);
       p.x = p.x - center.x;
       p.y = p.y + center.y;
       p = rotate(-angle, p);
       double u = (p.x + center.x) / (double) width;
       double v = (-1*(p.y - center.y)) / (double) height;
       n->GetPixel(x,y) = Sample(u,v);
       }
       }
       */
    return n;
}

void Image::Fun()
{
    Image img(width, height);
    int x,y;
    double cx = width / 2;
    double cy = height / 2;
    double miny = 100;
    double maxy = -100;
    for (y = 0; y < Height() ; y++) {
        for (x = 0; x < Width() ; x++) {
            double nx = (x*2.0) / width - 1;
            double ny = (y*2.0) / height - 1;
            double r = std::sqrt(nx*nx + ny*ny);
            if (r <= 1.0) {
                double theta = std::acos(-ny);
                double phi = std::asin(nx / std::sin(theta));
                if (theta == 0)
                    phi = 0;
                double xx = (phi + M_PI / 2) / M_PI;
                double yy = theta / M_PI;
                if (yy < miny)
                    miny = yy;
                if (yy > maxy)
                    maxy = yy;
                std::cout << "x: " << xx << ", y: " << yy << std::endl;
                img.GetPixel(x, y) = Sample(xx, yy);
            } else { 
                img.GetPixel(x,y) = Pixel(0,0,0,255);
            }
        }
    }
    std::cout << "miny: " << miny << ", maxy: " << maxy << std::endl;
    for (x = 0 ; x < Width() ; x++)
        for (y = 0 ; y < Height() ; y++)
            GetPixel(x,y) = img.GetPixel(x,y);

}

/**
 * Image Sample
 **/
void Image::SetSamplingMethod(int method)
{
    assert((method >= 0) && (method < IMAGE_N_SAMPLING_METHODS));
    sampling_method = method;
}


Pixel Image::Sample (double u, double v){
    if (u < 0 || u > 1 || v < 0 || v > 1)
        return Pixel(0, 0, 0, 0);

    Pixel p;
    double w = width*u;
    double h = height*v;
    switch(sampling_method) {
        case 0:
            {
                int x = std::min(width-1, (int) round(w));
                int y = std::min(height-1, (int) round(h));
                p = GetPixel(x,y);
            }
            break;
        case 1:
            {
                int x = std::floor(w);
                int y = std::floor(h);
                if (x == (width - 1))
                    x--;
                if (y == (height - 1))
                    y--;

                /*
                Pixel p1 = GetPixel(x,y);
                Pixel p2 = GetPixel(x+1,y);
                Pixel p3 = GetPixel(x,y+1);
                Pixel p4 = GetPixel(x+1,y+1);
                double dx = w - x;
                double dy = h - y;
                double fr = dy*(dx*p1.r + (1-dx)*p2.r) + (1-dy)*(dx*p3.r + (1-dx)*p4.r);
                double fg = dy*(dx*p1.g + (1-dx)*p2.g) + (1-dy)*(dx*p3.g + (1-dx)*p4.g);
                double fb = dy*(dx*p1.b + (1-dx)*p2.b) + (1-dy)*(dx*p3.b + (1-dx)*p4.b);
                p = Pixel(fr, fg, fb, 255);
                p.Clamp();
                return p;
                */

                // The LERP factors for x and y
                double tx = w - x;
                double tx2 = 1 - tx;
                double ty = h - y;
                double ty2 = 1-ty;
                Pixel tmp1 = GetPixel(x, y);
                Pixel tmp2 = GetPixel(x+1, y);
                // LERP horizontally on the top two points
                double p1r = tx*tmp1.r + tx2*tmp2.r;
                double p1g = tx*tmp1.g + tx2*tmp2.g;
                double p1b = tx*tmp1.b + tx2*tmp2.b;
                tmp1 = GetPixel(x, y+1);
                tmp2 = GetPixel(x+1, y+1);
                // LERP horizontally on the bottom two points
                double p2r = tx*tmp1.r + tx2*tmp2.r;
                double p2g = tx*tmp1.g + tx2*tmp2.g;
                double p2b = tx*tmp1.b + tx2*tmp2.b;
                // LERP vertically
                p = Pixel(ty*p1r + ty2*p2r, ty*p1g+ty2*p2g, ty*p1b+ty2*p1b, 255); 
                p.Clamp();
            }

            break;
        case IMAGE_SAMPLING_GAUSSIAN:
            {
                // Set up the kernel
                int SIZE = 5;
                double MAX = 5;
                double** kern = new double*[SIZE];
                for (int i = 0; i < SIZE; i++)
                    kern[i] = new double[SIZE];

                double total = 0;
                double center = SIZE / 2;
                double CAP = std::sqrt(2*center*center);
                for (int r = 0; r < SIZE; r++) {
                    for (int c = 0; c < SIZE; c++) {
                        double dx = -center + c;
                        double dy = -center + r;
                        kern[r][c] = MAX*(1 - (std::sqrt(dx*dx + dy*dy)/CAP));
                        total += kern[r][c];
                    }
                }
                for (int r = 0; r < SIZE; r++) {
                    for (int c = 0; c < SIZE; c++) {
                        kern[r][c] /= total;
                    }
                }
                int x = std::min(width-1, (int) round(w));
                int y = std::min(height-1, (int) round(h));
                p = ApplyKernel(x, y, SIZE, kern);
            }
            break;
    }

    return p;
}
