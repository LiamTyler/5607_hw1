#include "image.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <algorithm>
#include <iostream>
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
    auto f = [&](int c) {
        double dx = factor*(rand() / (float) RAND_MAX * 510 - 255);
        return ComponentClamp(c + dx);
    };

    for (x = 0 ; x < Width() ; x++)
    {
        for (y = 0 ; y < Height() ; y++)
        {
            Pixel p = GetPixel(x, y);
            Component nr = f(p.r);
            Component ng = f(p.g);
            Component nb = f(p.b);
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
            Component nr = ComponentClamp(128 + factor*(p.r - 128));
            Component ng = ComponentClamp(128 + factor*(p.g - 128));
            Component nb = ComponentClamp(128 + factor*(p.b - 128));
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
    for (int x = 0 ; x < Width() ; x++)
    {
        for (int y = 0 ; y < Height() ; y++)
        {
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
    for (int x = 0 ; x < Width(); x++)
    {
        int er = 0, eg = 0, eb = 0;
        for (int y = 0 ; y < Height() ; y++)
        {
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
    return NULL;
}

Image* Image::Rotate(double angle)
{
    /* WORK HERE */
    return NULL;
}

void Image::Fun()
{
    /* WORK HERE */
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
    /* WORK HERE */
    return Pixel();
}
