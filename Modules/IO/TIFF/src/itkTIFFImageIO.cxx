/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifdef _MSC_VER
#pragma warning( disable : 4611 )
#endif

#include "itkTIFFImageIO.h"
#include "itkRGBPixel.h"
#include "itkRGBAPixel.h"
#include "itksys/SystemTools.hxx"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/stat.h>

#include "itk_tiff.h"

namespace itk
{
class TIFFReaderInternal
{
public:
  TIFFReaderInternal();
  int Initialize();

  void Clean();

  int CanRead();

  int Open(const char *filename);

  TIFF *         m_Image;
  bool           m_IsOpen;
  unsigned int   m_Width;
  unsigned int   m_Height;
  unsigned short m_NumberOfPages;
  unsigned short m_CurrentPage;
  unsigned short m_SamplesPerPixel;
  unsigned short m_Compression;
  unsigned short m_BitsPerSample;
  unsigned short m_Photometrics;
  bool           m_HasValidPhotometricInterpretation;
  unsigned short m_PlanarConfig;
  unsigned short m_Orientation;
  uint32         m_TileDepth;
  unsigned int   m_TileRows;
  unsigned int   m_TileColumns;
  unsigned int   m_TileWidth;
  unsigned int   m_TileHeight;
  unsigned short m_NumberOfTiles;
  unsigned int   m_SubFiles;
  unsigned int   m_IgnoredSubFiles;
  unsigned int   m_ResolutionUnit;
  float          m_XResolution;
  float          m_YResolution;
  short          m_SampleFormat;
};

int TIFFReaderInternal::Open(const char *filename)
{
  this->Clean();
  struct stat fs;
  if ( stat(filename, &fs) )
    {
    return 0;
    }

  this->m_Image = TIFFOpen(filename, "r");
  if ( !this->m_Image )
    {
    this->Clean();
    return 0;
    }
  if ( !this->Initialize() )
    {
    this->Clean();
    return 0;
    }

  this->m_IsOpen = true;
  return 1;
}

void TIFFReaderInternal::Clean()
{
  if ( this->m_Image )
    {
    TIFFClose(this->m_Image);
    }
  this->m_Image = NULL;
  this->m_Width = 0;
  this->m_Height = 0;
  this->m_SamplesPerPixel = 0;
  this->m_Compression = 0;
  this->m_BitsPerSample = 0;
  this->m_Photometrics = 0;
  this->m_HasValidPhotometricInterpretation = false;
  this->m_PlanarConfig = 0;
  this->m_TileDepth = 0;
  this->m_CurrentPage = 0;
  this->m_NumberOfPages = 0;
  this->m_NumberOfTiles = 0;
  this->m_TileRows = 0;
  this->m_TileColumns = 0;
  this->m_TileWidth = 0;
  this->m_TileHeight = 0;
  this->m_XResolution = 1;
  this->m_YResolution = 1;
  this->m_SubFiles = 0;
  this->m_IgnoredSubFiles = 0;
  this->m_SampleFormat = 1;
  this->m_ResolutionUnit = 1; // none
  this->m_IsOpen = false;
}

TIFFReaderInternal::TIFFReaderInternal()
{
  this->m_Image           = NULL;
  this->Clean();
}

int TIFFReaderInternal::Initialize()
{
  if ( this->m_Image )
    {
    if ( !TIFFGetField(this->m_Image, TIFFTAG_IMAGEWIDTH, &this->m_Width)
         || !TIFFGetField(this->m_Image, TIFFTAG_IMAGELENGTH, &this->m_Height) )
      {
      return 0;
      }

    // Get the resolution in each direction
    TIFFGetField(this->m_Image,
                 TIFFTAG_XRESOLUTION, &this->m_XResolution);
    TIFFGetField(this->m_Image,
                 TIFFTAG_YRESOLUTION, &this->m_YResolution);
    TIFFGetField(this->m_Image,
                 TIFFTAG_RESOLUTIONUNIT, &this->m_ResolutionUnit);

    // Check the number of pages. First by looking at the number of directories
    this->m_NumberOfPages = TIFFNumberOfDirectories(this->m_Image);

    if ( this->m_NumberOfPages == 0 )
      {
      if ( !TIFFGetField(this->m_Image,
                         TIFFTAG_PAGENUMBER, &this->m_CurrentPage,
                         &this->m_NumberOfPages) )
        {
        // Check the Image Description tag to know the number of images
        // This is used by ImageJ
        char **description = new char *[255];
        if ( TIFFGetField(this->m_Image, TIFFTAG_IMAGEDESCRIPTION, description) )
          {
          // look for the number of images
          std::string desc = description[0];
          int         pos = desc.find("images=");
          int         pos2 = desc.find("\n");
          if ( ( pos != -1 ) && ( pos2 != -1 ) )
            {
            this->m_NumberOfPages = atoi( desc.substr(pos + 7, pos2 - pos - 7).c_str() );
            }
          }
        }
      }

    // If the number of pages is still zero we look if the image is tiled
    if ( this->m_NumberOfPages == 0 && TIFFIsTiled(this->m_Image) )
      {
      this->m_NumberOfTiles = TIFFNumberOfTiles(this->m_Image);

      if ( !TIFFGetField(this->m_Image, TIFFTAG_TILEWIDTH, &this->m_TileWidth)
           || !TIFFGetField(this->m_Image, TIFFTAG_TILELENGTH, &this->m_TileHeight) )
        {
        itkGenericExceptionMacro(
          << "Cannot read tile width and tile length from file");
        }
      else
        {
        this->m_TileRows = this->m_Height / this->m_TileHeight;
        this->m_TileColumns = this->m_Width / this->m_TileWidth;
        }
      }

    // Checking if the TIFF contains subfiles
    if ( this->m_NumberOfPages > 1 )
      {
      this->m_SubFiles = 0;
      this->m_IgnoredSubFiles = 0;

      for ( unsigned int page = 0; page < this->m_NumberOfPages; page++ )
        {
        int32 subfiletype = 6;
        if ( TIFFGetField(this->m_Image, TIFFTAG_SUBFILETYPE, &subfiletype) )
          {
          if ( subfiletype == 0 )
            {
            this->m_SubFiles += 1;
            }
          // ignored flags
          else if ( subfiletype & FILETYPE_REDUCEDIMAGE
                    || subfiletype & FILETYPE_MASK )
            {
            ++this->m_IgnoredSubFiles;
            }

          }
        TIFFReadDirectory(this->m_Image);
        }

      // Set the directory to the first image, and reads it
      TIFFSetDirectory(this->m_Image, 0);
      }

    TIFFGetFieldDefaulted(this->m_Image, TIFFTAG_ORIENTATION,
                          &this->m_Orientation);
    TIFFGetFieldDefaulted(this->m_Image, TIFFTAG_SAMPLESPERPIXEL,
                          &this->m_SamplesPerPixel);
    TIFFGetFieldDefaulted(this->m_Image, TIFFTAG_COMPRESSION, &this->m_Compression);
    TIFFGetFieldDefaulted(this->m_Image, TIFFTAG_BITSPERSAMPLE,
                          &this->m_BitsPerSample);
    TIFFGetFieldDefaulted(this->m_Image, TIFFTAG_PLANARCONFIG, &this->m_PlanarConfig);
    TIFFGetFieldDefaulted(this->m_Image, TIFFTAG_SAMPLEFORMAT, &this->m_SampleFormat);

    // If TIFFGetField returns false, there's no Photometric Interpretation
    // set for this image, but that's a required field so we set a warning flag.
    // (Because the "Photometrics" field is an enum, we can't rely on setting
    // this->m_Photometrics to some signal value.)
    if ( TIFFGetField(this->m_Image, TIFFTAG_PHOTOMETRIC, &this->m_Photometrics) )
      {
      this->m_HasValidPhotometricInterpretation = true;
      }
    else
      {
      this->m_HasValidPhotometricInterpretation = false;
      }
    if ( !TIFFGetField(this->m_Image, TIFFTAG_TILEDEPTH, &this->m_TileDepth) )
      {
      this->m_TileDepth = 0;
      }
    }

  return 1;
}

int TIFFReaderInternal::CanRead()
{
  return ( this->m_Image && ( this->m_Width > 0 ) && ( this->m_Height > 0 )
           && ( this->m_SamplesPerPixel > 0 )
           && ( this->m_Compression == COMPRESSION_NONE
                || this->m_Compression == COMPRESSION_PACKBITS
                || this->m_Compression == COMPRESSION_LZW
                )
           && ( this->m_HasValidPhotometricInterpretation )
           && ( this->m_Photometrics == PHOTOMETRIC_RGB
                || this->m_Photometrics == PHOTOMETRIC_MINISWHITE
                || this->m_Photometrics == PHOTOMETRIC_MINISBLACK
                || this->m_Photometrics == PHOTOMETRIC_PALETTE )
           && ( this->m_PlanarConfig == PLANARCONFIG_CONTIG )
           && ( !this->m_TileDepth )
           && ( this->m_BitsPerSample == 8 || this->m_BitsPerSample == 16 ) );
}

bool TIFFImageIO::CanReadFile(const char *file)
{
  // First check the extension
  std::string filename = file;

  if (  filename == "" )
    {
    itkDebugMacro(<< "No filename specified.");
    return false;
    }

  // Now check if this is a valid TIFF image
  TIFFErrorHandler save = TIFFSetErrorHandler(0);
  int              res = m_InternalImage->Open(file);
  if ( res )
    {
    TIFFSetErrorHandler(save);
    return true;
    }
  m_InternalImage->Clean();
  TIFFSetErrorHandler(save);
  return false;
}

/** To Support Zeiss images that contains only 2 samples per pixel but
 *  are actually RGB images */
void TIFFImageIO::ReadTwoSamplesPerPixelImage(void *out,
                                              unsigned int width,
                                              unsigned int height)
{
  unsigned int isize = TIFFScanlineSize(m_InternalImage->m_Image);
  unsigned int cc;
  int          row;
  tdata_t      buf = _TIFFmalloc(isize);

  int inc = 1;

  if ( m_ComponentType == UCHAR )
    {
    unsigned char *image;
    if ( m_InternalImage->m_PlanarConfig == PLANARCONFIG_CONTIG )
      {
      for ( row = 0; row < (int)height; row++ )
        {
        if ( TIFFReadScanline(m_InternalImage->m_Image, buf, row, 0) <= 0 )
          {
          itkExceptionMacro(<< "Problem reading the row: " << row);
          break;
          }

        if ( m_InternalImage->m_Orientation == ORIENTATION_TOPLEFT )
          {
          image = reinterpret_cast< unsigned char * >( out ) + row * width * inc;
          }
        else
          {
          image = reinterpret_cast< unsigned char * >( out ) + width * inc * ( height - ( row + 1 ) );
          }

        for ( cc = 0; cc < isize;
              cc += m_InternalImage->m_SamplesPerPixel )
          {
          inc = this->EvaluateImageAt(image,
                                      static_cast< unsigned char * >( buf )
                                      + cc);
          image += inc;
          }
        }
      }
    else if ( m_InternalImage->m_PlanarConfig == PLANARCONFIG_SEPARATE )
      {
      uint32 s;
      uint32 nsamples = 0;
      TIFFGetField(m_InternalImage->m_Image, TIFFTAG_SAMPLESPERPIXEL, &nsamples);
      for ( s = 0; s < nsamples; s++ )
        {
        for ( row = 0; row < (int)height; row++ )
          {
          if ( TIFFReadScanline(m_InternalImage->m_Image, buf, row, s) <= 0 )
            {
            itkExceptionMacro(<< "Problem reading the row: " << row);
            break;
            }

          inc = 3;

          if ( m_InternalImage->m_Orientation == ORIENTATION_TOPLEFT )
            {
            image = reinterpret_cast< unsigned char * >( out ) + row * width * inc;
            }
          else
            {
            image = reinterpret_cast< unsigned char * >( out ) + width * inc * ( height - ( row + 1 ) );
            }

          // We translate the output pixel to be on the right RGB
          image += s;
          for ( cc = 0; cc < isize;
                cc += 1 )
            {
            ( *image ) = *( static_cast< unsigned char * >( buf ) + cc );
            inc = 3;
            image += inc;
            }
          }
        }
      }
    }
  else if ( m_ComponentType == USHORT )
    {
    isize /= 2;
    unsigned short *image;
    if ( m_InternalImage->m_PlanarConfig == PLANARCONFIG_CONTIG )
      {
      for ( row = 0; row < (int)height; row++ )
        {
        if ( TIFFReadScanline(m_InternalImage->m_Image, buf, row, 0) <= 0 )
          {
          itkExceptionMacro(<< "Problem reading the row: " << row);
          break;
          }

        if ( m_InternalImage->m_Orientation == ORIENTATION_TOPLEFT )
          {
          image = reinterpret_cast< unsigned short * >( out ) + row * width * inc;
          }
        else
          {
          image = reinterpret_cast< unsigned short * >( out ) + width * inc * ( height - ( row + 1 ) );
          }

        for ( cc = 0; cc < isize;
              cc += m_InternalImage->m_SamplesPerPixel )
          {
          inc = this->EvaluateImageAt(image,
                                      static_cast< unsigned short * >( buf )
                                      + cc);
          image += inc;
          }
        }
      }
    else if ( m_InternalImage->m_PlanarConfig == PLANARCONFIG_SEPARATE )
      {
      uint32 s, nsamples;
      TIFFGetField(m_InternalImage->m_Image, TIFFTAG_SAMPLESPERPIXEL, &nsamples);
      for ( s = 0; s < nsamples; s++ )
        {
        for ( row = 0; row < (int)height; row++ )
          {
          if ( TIFFReadScanline(m_InternalImage->m_Image, buf, row, s) <= 0 )
            {
            itkExceptionMacro(<< "Problem reading the row: " << row);
            break;
            }

          if ( m_InternalImage->m_Orientation == ORIENTATION_TOPLEFT )
            {
            image = reinterpret_cast< unsigned short * >( out ) + row * width * inc;
            }
          else
            {
            image = reinterpret_cast< unsigned short * >( out ) + width * inc * ( height - ( row + 1 ) );
            }
          // We translate the output pixel to be on the right RGB
          image += s;
          for ( cc = 0; cc < isize;
                cc += 1 )
            {
            ( *image ) = *( static_cast< unsigned short * >( buf ) + cc );
            inc = 3;
            image += inc;
            }
          }
        }
      }
    }
  _TIFFfree(buf);
}

void TIFFImageIO::ReadGenericImage(void *out,
                                   unsigned int width,
                                   unsigned int height)
{
  unsigned int isize = TIFFScanlineSize(m_InternalImage->m_Image);
  unsigned int cc;
  int          row, inc;
  tdata_t      buf = _TIFFmalloc(isize);

  // It is necessary to re-initialize the colors for eachread so
  // that the colormap remains valid.
  this->InitializeColors();

  if ( m_InternalImage->m_PlanarConfig != PLANARCONFIG_CONTIG )
    {
    itkExceptionMacro(<< "This reader can only do PLANARCONFIG_CONTIG");
    return;
    }

  switch ( this->GetFormat() )
    {
    case TIFFImageIO::GRAYSCALE:
    case TIFFImageIO::PALETTE_GRAYSCALE:
      inc = 1;
      break;
    case TIFFImageIO::RGB_:
      inc = m_InternalImage->m_SamplesPerPixel;
      break;
    case TIFFImageIO::PALETTE_RGB:
      inc = 3;
      break;
    default:
      inc = 1;
      break;
    }

  if ( m_ComponentType == UCHAR )
    {
    unsigned char *image;
    if ( m_InternalImage->m_PlanarConfig == PLANARCONFIG_CONTIG )
      {
      for ( row = 0; row < (int)height; row++ )
        {
        if ( TIFFReadScanline(m_InternalImage->m_Image, buf, row, 0) <= 0 )
          {
          itkExceptionMacro(<< "Problem reading the row: " << row);
          break;
          }

        if ( m_InternalImage->m_Orientation == ORIENTATION_TOPLEFT )
          {
          image = reinterpret_cast< unsigned char * >( out ) + row * width * inc;
          }
        else
          {
          image = reinterpret_cast< unsigned char * >( out ) + width * inc * ( height - ( row + 1 ) );
          }

        for ( cc = 0; cc < isize;
              cc += m_InternalImage->m_SamplesPerPixel )
          {
          inc = this->EvaluateImageAt(image,
                                      static_cast< unsigned char * >( buf )
                                      + cc);
          image += inc;
          }
        }
      }
    else if ( m_InternalImage->m_PlanarConfig == PLANARCONFIG_SEPARATE )
      {
      uint32 s;
      uint32 nsamples;
      TIFFGetField(m_InternalImage->m_Image, TIFFTAG_SAMPLESPERPIXEL, &nsamples);
      for ( s = 0; s < nsamples; s++ )
        {
        for ( row = 0; row < (int)height; row++ )
          {
          if ( TIFFReadScanline(m_InternalImage->m_Image, buf, row, s) <= 0 )
            {
            itkExceptionMacro(<< "Problem reading the row: " << row);
            break;
            }

          inc = 3;
          if ( m_InternalImage->m_Orientation == ORIENTATION_TOPLEFT )
            {
            image = reinterpret_cast< unsigned char * >( out ) + row * width * inc;
            }
          else
            {
            image = reinterpret_cast< unsigned char * >( out ) + width * inc * ( height - ( row + 1 ) );
            }

          for ( cc = 0; cc < isize;
                cc += m_InternalImage->m_SamplesPerPixel )
            {
            inc = this->EvaluateImageAt(image,
                                        static_cast< unsigned char * >( buf )
                                        + cc);
            image += inc;
            }
          }
        }
      }
    }
  else if ( m_ComponentType == CHAR )
    {
    char *image;
    if ( m_InternalImage->m_PlanarConfig == PLANARCONFIG_CONTIG )
      {
      for ( row = 0; row < (int)height; row++ )
        {
        if ( TIFFReadScanline(m_InternalImage->m_Image, buf, row, 0) <= 0 )
          {
          itkExceptionMacro(<< "Problem reading the row: " << row);
          break;
          }

        if ( m_InternalImage->m_Orientation == ORIENTATION_TOPLEFT )
          {
          image = reinterpret_cast< char * >( out ) + row * width * inc;
          }
        else
          {
          image = reinterpret_cast< char * >( out ) + width * inc * ( height - ( row + 1 ) );
          }

        for ( cc = 0; cc < isize;
              cc += m_InternalImage->m_SamplesPerPixel )
          {
          inc = this->EvaluateImageAt(image,
                                      static_cast< char * >( buf )
                                      + cc);
          image += inc;
          }
        }
      }
    else if ( m_InternalImage->m_PlanarConfig == PLANARCONFIG_SEPARATE )
      {
      uint32 s;
      uint32 nsamples;
      TIFFGetField(m_InternalImage->m_Image, TIFFTAG_SAMPLESPERPIXEL, &nsamples);
      for ( s = 0; s < nsamples; s++ )
        {
        for ( row = 0; row < (int)height; row++ )
          {
          if ( TIFFReadScanline(m_InternalImage->m_Image, buf, row, s) <= 0 )
            {
            itkExceptionMacro(<< "Problem reading the row: " << row);
            break;
            }

          inc = 3;
          if ( m_InternalImage->m_Orientation == ORIENTATION_TOPLEFT )
            {
            image = reinterpret_cast< char * >( out ) + row * width * inc;
            }
          else
            {
            image = reinterpret_cast< char * >( out ) + width * inc * ( height - ( row + 1 ) );
            }

          for ( cc = 0; cc < isize;
                cc += m_InternalImage->m_SamplesPerPixel )
            {
            inc = this->EvaluateImageAt(image,
                                        static_cast< char * >( buf )
                                        + cc);
            image += inc;
            }
          }
        }
      }
    }
  else if ( m_ComponentType == USHORT )
    {
    isize /= 2;
    unsigned short *image;
    if ( m_InternalImage->m_PlanarConfig == PLANARCONFIG_CONTIG )
      {
      for ( row = 0; row < (int)height; row++ )
        {
        if ( TIFFReadScanline(m_InternalImage->m_Image, buf, row, 0) <= 0 )
          {
          itkExceptionMacro(<< "Problem reading the row: " << row);
          break;
          }

        if ( m_InternalImage->m_Orientation == ORIENTATION_TOPLEFT )
          {
          image = reinterpret_cast< unsigned short * >( out ) + row * width * inc;
          }
        else
          {
          image = reinterpret_cast< unsigned short * >( out ) + width * inc * ( height - ( row + 1 ) );
          }

        for ( cc = 0; cc < isize;
              cc += m_InternalImage->m_SamplesPerPixel )
          {
          inc = this->EvaluateImageAt(image,
                                      static_cast< unsigned short * >( buf )
                                      + cc);
          image += inc;
          }
        }
      }
    else if ( m_InternalImage->m_PlanarConfig == PLANARCONFIG_SEPARATE )
      {
      uint32 s, nsamples;
      TIFFGetField(m_InternalImage->m_Image, TIFFTAG_SAMPLESPERPIXEL, &nsamples);
      for ( s = 0; s < nsamples; s++ )
        {
        for ( row = 0; row < (int)height; row++ )
          {
          if ( TIFFReadScanline(m_InternalImage->m_Image, buf, row, s) <= 0 )
            {
            itkExceptionMacro(<< "Problem reading the row: " << row);
            break;
            }

          if ( m_InternalImage->m_Orientation == ORIENTATION_TOPLEFT )
            {
            image = reinterpret_cast< unsigned short * >( out ) + row * width * inc;
            }
          else
            {
            image = reinterpret_cast< unsigned short * >( out ) + width * inc * ( height - ( row + 1 ) );
            }
          for ( cc = 0; cc < isize;
                cc += m_InternalImage->m_SamplesPerPixel )
            {
            inc = this->EvaluateImageAt(image,
                                        static_cast< unsigned short * >( buf )
                                        + cc);
            image += inc;
            }
          }
        }
      }
    }
  // Short type
  else if ( m_ComponentType == SHORT )
    {
    isize /= 2;
    short *image;
    if ( m_InternalImage->m_PlanarConfig == PLANARCONFIG_CONTIG )
      {
      for ( row = 0; row < (int)height; row++ )
        {
        if ( TIFFReadScanline(m_InternalImage->m_Image, buf, row, 0) <= 0 )
          {
          itkExceptionMacro(<< "Problem reading the row: " << row);
          break;
          }

        if ( m_InternalImage->m_Orientation == ORIENTATION_TOPLEFT )
          {
          image = reinterpret_cast< short * >( out ) + row * width * inc;
          }
        else
          {
          image = reinterpret_cast< short * >( out ) + width * inc * ( height - ( row + 1 ) );
          }

        for ( cc = 0; cc < isize;
              cc += m_InternalImage->m_SamplesPerPixel )
          {
          inc = this->EvaluateImageAt(image,
                                      static_cast< short * >( buf )
                                      + cc);
          image += inc;
          }
        }
      }
    else if ( m_InternalImage->m_PlanarConfig == PLANARCONFIG_SEPARATE )
      {
      uint32 s, nsamples;
      TIFFGetField(m_InternalImage->m_Image, TIFFTAG_SAMPLESPERPIXEL, &nsamples);
      for ( s = 0; s < nsamples; s++ )
        {
        for ( row = 0; row < (int)height; row++ )
          {
          if ( TIFFReadScanline(m_InternalImage->m_Image, buf, row, s) <= 0 )
            {
            itkExceptionMacro(<< "Problem reading the row: " << row);
            break;
            }

          if ( m_InternalImage->m_Orientation == ORIENTATION_TOPLEFT )
            {
            image = reinterpret_cast< short * >( out ) + row * width * inc;
            }
          else
            {
            image = reinterpret_cast< short * >( out ) + width * inc * ( height - ( row + 1 ) );
            }
          for ( cc = 0; cc < isize;
                cc += m_InternalImage->m_SamplesPerPixel )
            {
            inc = this->EvaluateImageAt(image,
                                        static_cast< short * >( buf )
                                        + cc);
            image += inc;
            }
          }
        }
      }
    }
  _TIFFfree(buf);
}

int TIFFImageIO::EvaluateImageAt(void *out, void *in)
{
  unsigned char *image = (unsigned char *)out;
  unsigned char *source = (unsigned char *)in;

  int            increment;
  unsigned short red, green, blue, alpha;

  switch ( this->GetFormat() )
    {
    case TIFFImageIO::GRAYSCALE:
      if ( m_InternalImage->m_Photometrics ==
           PHOTOMETRIC_MINISBLACK )
        {
        if ( m_ComponentType == USHORT )
          {
          unsigned short *image_us = (unsigned short *)out;
          unsigned short *source_us = (unsigned short *)in;
          *image_us = *source_us;
          }
        else if ( m_ComponentType == SHORT )
          {
          short *image_us = (short *)out;
          short *source_us = (short *)in;
          *image_us = *source_us;
          }
        else if ( m_ComponentType == CHAR )
          {
          char *image_us = (char *)out;
          char *source_us = (char *)in;
          *image_us = *source_us;
          }
        else
          {
          *image = *source;
          }
        }
      else
        {
        *image = ~( *source );
        }
      increment = 1;
      break;
    case TIFFImageIO::PALETTE_GRAYSCALE:
      this->GetColor(*source, &red, &green, &blue);
      *image = static_cast< unsigned char >( red >> 8 );
      increment = 1;
      break;
    case TIFFImageIO::RGB_:
      if ( m_ComponentType == USHORT )
        {
        unsigned short *image_us = (unsigned short *)out;
        unsigned short *source_us = (unsigned short *)in;

        red   = *( source_us );
        green = *( source_us + 1 );
        blue  = *( source_us + 2 );
        *( image_us )   = red;
        *( image_us + 1 ) = green;
        *( image_us + 2 ) = blue;
        if ( m_InternalImage->m_SamplesPerPixel == 4 )
          {
          alpha = *( source_us + 3 );
          *( image_us + 3 ) = 65535 - alpha;
          }
        }
      else
        {
        red   = *( source );
        green = *( source + 1 );
        blue  = *( source + 2 );
        *( image )   = red;
        *( image + 1 ) = green;
        *( image + 2 ) = blue;
        if ( m_InternalImage->m_SamplesPerPixel == 4 )
          {
          alpha = *( source + 3 );
          *( image + 3 ) = 255 - alpha;
          }
        }
      increment = m_InternalImage->m_SamplesPerPixel;
      break;
    case TIFFImageIO::PALETTE_RGB:
      if ( m_ComponentType == USHORT )
        {
        unsigned short *image_us = (unsigned short *)out;
        unsigned short *source_us = (unsigned short *)in;
        this->GetColor(*source_us, &red, &green, &blue);
        *( image_us )   = red << 8;
        *( image_us + 1 ) = green << 8;
        *( image_us + 2 ) = blue << 8;
        }
      else if ( m_ComponentType == SHORT )
        {
        short *image_us = (short *)out;
        short *source_us = (short *)in;
        this->GetColor(*source_us, &red, &green, &blue);
        *( image_us )   = red << 8;
        *( image_us + 1 ) = green << 8;
        *( image_us + 2 ) = blue << 8;
        }
      else if ( m_ComponentType == CHAR )
        {
        this->GetColor(*source, &red, &green, &blue);
        *( image )   = static_cast< char >( red >> 8 );
        *( image + 1 ) = static_cast< char >( green >> 8 );
        *( image + 2 ) = static_cast< char >( blue >> 8 );
        }
      else
        {
        this->GetColor(*source, &red, &green, &blue);
        *( image )   = static_cast< unsigned char >( red >> 8 );
        *( image + 1 ) = static_cast< unsigned char >( green >> 8 );
        *( image + 2 ) = static_cast< unsigned char >( blue >> 8 );
        }
      increment = 3;
      break;
    default:
      return 0;
    }

  return increment;
}

void TIFFImageIO::GetColor(int index, unsigned short *red,
                           unsigned short *green, unsigned short *blue)
{
  *red   = 0;
  *green = 0;
  *blue  = 0;
  if ( index < 0 )
    {
    itkExceptionMacro(<< "Color index has to be greater than 0");
    return;
    }
  if ( m_TotalColors > 0
       && m_ColorRed && m_ColorGreen && m_ColorBlue )
    {
    if ( index >= m_TotalColors )
      {
      itkExceptionMacro(<< "Color index has to be less than number of colors ("
                        << m_TotalColors << ")");
      return;
      }
    *red   = *( m_ColorRed   + index );
    *green = *( m_ColorGreen + index );
    *blue  = *( m_ColorBlue  + index );
    return;
    }

  unsigned short photometric;

  if ( !TIFFGetField(m_InternalImage->m_Image, TIFFTAG_PHOTOMETRIC, &photometric) )
    {
    if ( m_InternalImage->m_Photometrics != PHOTOMETRIC_PALETTE )
      {
      itkExceptionMacro(<< "You can only access colors for palette images");
      return;
      }
    }

  unsigned short *red_orig, *green_orig, *blue_orig;

  switch ( m_InternalImage->m_BitsPerSample )
    {
    case 1:
    case 2:
    case 4:
    case 8:
    case 16:
      break;
    default:
      itkExceptionMacro(<<  "Sorry, can not handle image with "
                        << m_InternalImage->m_BitsPerSample
                        << "-bit samples");
      return;
    }
  if ( !TIFFGetField(m_InternalImage->m_Image, TIFFTAG_COLORMAP,
                     &red_orig, &green_orig, &blue_orig) )
    {
    itkExceptionMacro(<< "Missing required \"Colormap\" tag");
    return;
    }
  m_TotalColors = ( 1L << m_InternalImage->m_BitsPerSample );

  if ( index >= m_TotalColors )
    {
    itkExceptionMacro(<< "Color index has to be less than number of colors ("
                      << m_TotalColors << ")");
    return;
    }
  m_ColorRed   =   red_orig;
  m_ColorGreen = green_orig;
  m_ColorBlue  =  blue_orig;

  *red   = *( red_orig   + index );
  *green = *( green_orig + index );
  *blue  = *( blue_orig  + index );
}

unsigned int TIFFImageIO::GetFormat()
{
  unsigned int cc;

  if ( m_ImageFormat != TIFFImageIO::NOFORMAT )
    {
    return m_ImageFormat;
    }

  switch ( m_InternalImage->m_Photometrics )
    {
    case PHOTOMETRIC_RGB:
    case PHOTOMETRIC_YCBCR:
      m_ImageFormat = TIFFImageIO::RGB_;
      return m_ImageFormat;
    case PHOTOMETRIC_MINISWHITE:
    case PHOTOMETRIC_MINISBLACK:
      m_ImageFormat = TIFFImageIO::GRAYSCALE;
      return m_ImageFormat;
    case PHOTOMETRIC_PALETTE:
      for ( cc = 0; cc < 256; cc++ )
        {
        unsigned short red, green, blue;
        this->GetColor(cc, &red, &green, &blue);
        if ( red != green || red != blue )
          {
          m_ImageFormat = TIFFImageIO::PALETTE_RGB;
          return m_ImageFormat;
          }
        }
      m_ImageFormat = TIFFImageIO::PALETTE_GRAYSCALE;
      return m_ImageFormat;
    }
  m_ImageFormat = TIFFImageIO::OTHER;
  return m_ImageFormat;
}

/** Read a tiled tiff */
void TIFFImageIO::ReadTiles(void *buffer)
{
  unsigned char *volume = reinterpret_cast< unsigned char * >( buffer );

  for ( unsigned int col = 0; col < m_InternalImage->m_Width; col += m_InternalImage->m_TileWidth )
    {
    for ( unsigned int row = 0; row < m_InternalImage->m_Height; row += m_InternalImage->m_TileHeight )
      {
      unsigned char *tempImage;
      tempImage =
        new unsigned char[m_InternalImage->m_TileWidth * m_InternalImage->m_TileHeight
                          * m_InternalImage->m_SamplesPerPixel];

      if ( TIFFReadTile(m_InternalImage->m_Image, tempImage, col, row, 0, 0) < 0 )
        {
        itkExceptionMacro(<< "Cannot read tile : " << row << "," << col << " from file");
        if ( tempImage != buffer )
          {
          delete[] tempImage;
          }

        return;
        }

      unsigned int xx, yy;
      for ( yy = 0; yy < m_InternalImage->m_TileHeight; yy++ )
        {
        for ( xx = 0; xx <  m_InternalImage->m_TileWidth; xx++ )
          {
          for ( unsigned int i = 0; i < m_InternalImage->m_SamplesPerPixel; i++ )
            {
            *volume = *( tempImage++ );
            volume++;
            }
          }
        }
      }
    }
}

/** Read a multipage tiff */
void TIFFImageIO::ReadVolume(void *buffer)
{
  int width  = m_InternalImage->m_Width;
  int height = m_InternalImage->m_Height;

  for ( unsigned int page = 0; page < m_InternalImage->m_NumberOfPages; page++ )
    {
    if ( m_InternalImage->m_IgnoredSubFiles > 0 )
      {
      int32 subfiletype = 6;
      if ( TIFFGetField(m_InternalImage->m_Image, TIFFTAG_SUBFILETYPE, &subfiletype) )
        {
        if ( subfiletype & FILETYPE_REDUCEDIMAGE
             || subfiletype & FILETYPE_MASK )
          {
          // skip subfile
          TIFFReadDirectory(m_InternalImage->m_Image);
          continue;
          }
        }
      }

    // It is necessary to re-initialize the colors for each page so
    // that the colormap is reset in the GetColor method.  This is
    // also true in the case that each slice has a different colormap.
    this->InitializeColors();

    // if we have a Zeiss image meaning that the SamplesPerPixel is 2
    if ( m_InternalImage->m_SamplesPerPixel == 2 )
      {
      if ( m_ComponentType == USHORT )
        {
        unsigned short *volume = reinterpret_cast< unsigned short * >( buffer );
        volume += width * height * m_InternalImage->m_SamplesPerPixel * page;
        this->ReadTwoSamplesPerPixelImage(volume, width, height);
        }
      else if ( m_ComponentType == SHORT )
        {
        short *volume = reinterpret_cast< short * >( buffer );
        volume += width * height * m_InternalImage->m_SamplesPerPixel * page;
        this->ReadTwoSamplesPerPixelImage(volume, width, height);
        }
      else if ( m_ComponentType == CHAR )
        {
        char *volume = reinterpret_cast< char * >( buffer );
        volume += width * height * m_InternalImage->m_SamplesPerPixel * page;
        this->ReadTwoSamplesPerPixelImage(volume, width, height);
        }
      else
        {
        unsigned char *volume = reinterpret_cast< unsigned char * >( buffer );
        volume += width * height * m_InternalImage->m_SamplesPerPixel * page;
        this->ReadTwoSamplesPerPixelImage(volume, width, height);
        }
      break;
      }
    else if ( !m_InternalImage->CanRead() )
      {
      uint32 *tempImage;
      tempImage = new uint32[width * height];

      if ( !TIFFReadRGBAImage(m_InternalImage->m_Image,
                              width, height,
                              tempImage, 1) )
        {
        itkExceptionMacro(<< "Cannot read TIFF image or as a TIFF RGBA image");
        if ( tempImage != buffer )
          {
          delete[] tempImage;
          }
        return;
        }
      int     xx, yy;
      uint32 *ssimage;

      if ( m_ComponentType == USHORT )
        {
        unsigned short *fimage = (unsigned short *)buffer;
        fimage += width * height * 4 * page;
        for ( yy = 0; yy < height; yy++ )
          {
          ssimage = tempImage + ( height - yy - 1 ) * width;
          for ( xx = 0; xx < width; xx++ )
            {
            unsigned short red   = static_cast< unsigned short >( TIFFGetR(*ssimage) );
            unsigned short green = static_cast< unsigned short >( TIFFGetG(*ssimage) );
            unsigned short blue  = static_cast< unsigned short >( TIFFGetB(*ssimage) );
            unsigned short alpha = static_cast< unsigned short >( TIFFGetA(*ssimage) );

            *( fimage  ) = red;
            *( fimage + 1 ) = green;
            *( fimage + 2 ) = blue;
            *( fimage + 3 ) = alpha;
            fimage += 4;
            ssimage++;
            }
          }
        }
      else if ( m_ComponentType == SHORT )
        {
        short *fimage = (short *)buffer;
        fimage += width * height * 4 * page;
        for ( yy = 0; yy < height; yy++ )
          {
          ssimage = tempImage + ( height - yy - 1 ) * width;
          for ( xx = 0; xx < width; xx++ )
            {
            short red   = static_cast< short >( TIFFGetR(*ssimage) );
            short green = static_cast< short >( TIFFGetG(*ssimage) );
            short blue  = static_cast< short >( TIFFGetB(*ssimage) );
            short alpha = static_cast< short >( TIFFGetA(*ssimage) );

            *( fimage  ) = red;
            *( fimage + 1 ) = green;
            *( fimage + 2 ) = blue;
            *( fimage + 3 ) = alpha;
            fimage += 4;
            ssimage++;
            }
          }
        }
      else if ( m_ComponentType == CHAR )
        {
        char *fimage = (char *)buffer;
        fimage += width * height * 4 * page;
        for ( yy = 0; yy < height; yy++ )
          {
          ssimage = tempImage + ( height - yy - 1 ) * width;
          for ( xx = 0; xx < width; xx++ )
            {
            char red   = static_cast< char >( TIFFGetR(*ssimage) );
            char green = static_cast< char >( TIFFGetG(*ssimage) );
            char blue  = static_cast< char >( TIFFGetB(*ssimage) );
            char alpha = static_cast< char >( TIFFGetA(*ssimage) );

            *( fimage  ) = red;
            *( fimage + 1 ) = green;
            *( fimage + 2 ) = blue;
            *( fimage + 3 ) = alpha;
            fimage += 4;
            ssimage++;
            }
          }
        }
      else
        {
        unsigned char *fimage = (unsigned char *)buffer;
        fimage += width * height * 4 * page / 2;
        for ( yy = 0; yy < height; yy++ )
          {
          ssimage = tempImage + ( height - yy - 1 ) * width;
          for ( xx = 0; xx < width; xx++ )
            {
            unsigned char red   = static_cast< unsigned char >( TIFFGetR(*ssimage) );
            unsigned char green = static_cast< unsigned char >( TIFFGetG(*ssimage) );
            unsigned char blue  = static_cast< unsigned char >( TIFFGetB(*ssimage) );
            unsigned char alpha = static_cast< unsigned char >( TIFFGetA(*ssimage) );

            *( fimage  ) = red;
            *( fimage + 1 ) = green;
            *( fimage + 2 ) = blue;
            *( fimage + 3 ) = alpha;
            fimage += 4;
            ssimage++;
            }
          }
        }
      if ( tempImage != 0 && tempImage != buffer )
        {
        delete[] tempImage;
        }
      }
    else
      {
      unsigned int format = this->GetFormat();

      switch ( format )
        {
        case TIFFImageIO::GRAYSCALE:
        case TIFFImageIO::RGB_:
        case TIFFImageIO::PALETTE_GRAYSCALE:
          if ( m_ComponentType == USHORT )
            {
            unsigned short *volume = reinterpret_cast< unsigned short * >( buffer );
            volume += width * height * m_InternalImage->m_SamplesPerPixel * page;
            this->ReadGenericImage(volume, width, height);
            }
          else if ( m_ComponentType == SHORT )
            {
            short *volume = reinterpret_cast< short * >( buffer );
            volume += width * height * m_InternalImage->m_SamplesPerPixel * page;
            this->ReadGenericImage(volume, width, height);
            }
          else if ( m_ComponentType == CHAR )
            {
            char *volume = reinterpret_cast< char * >( buffer );
            volume += width * height * m_InternalImage->m_SamplesPerPixel * page;
            this->ReadGenericImage(volume, width, height);
            }
          else
            {
            unsigned char *volume = reinterpret_cast< unsigned char * >( buffer );
            volume += width * height * m_InternalImage->m_SamplesPerPixel * page;
            this->ReadGenericImage(volume, width, height);
            }
          break;
        case TIFFImageIO::PALETTE_RGB:
          // This differs from PALLETTE_GRAYSCALE only in that the
          // volume is incremented by 3 times more since the colormap
          // consists of RGB.
          if ( m_ComponentType == USHORT )
            {
            unsigned short *volume = reinterpret_cast< unsigned short * >( buffer );
            volume += width * height * m_InternalImage->m_SamplesPerPixel * page * 3;
            this->ReadGenericImage(volume, width, height);
            }
          else if ( m_ComponentType == SHORT )
            {
            short *volume = reinterpret_cast< short * >( buffer );
            volume += width * height * m_InternalImage->m_SamplesPerPixel * page * 3;
            this->ReadGenericImage(volume, width, height);
            }
          else if ( m_ComponentType == CHAR )
            {
            char *volume = reinterpret_cast< char * >( buffer );
            volume += width * height * m_InternalImage->m_SamplesPerPixel * page * 3;
            this->ReadGenericImage(volume, width, height);
            }
          else
            {
            unsigned char *volume = reinterpret_cast< unsigned char * >( buffer );
            volume += width * height * m_InternalImage->m_SamplesPerPixel * page * 3;
            this->ReadGenericImage(volume, width, height);
            }
          break;
        default:
          return;
        }
      }
    TIFFReadDirectory(m_InternalImage->m_Image);
    }
}

void TIFFImageIO::Read(void *buffer)
{

  // re-open the file if it was closed
  if ( !m_InternalImage->m_IsOpen )
    {
    if ( !this->CanReadFile( m_FileName.c_str() ) )
      {
      itkExceptionMacro(<< "Cannot open file " << this->m_FileName << "!");
      return;
      }
    }

  if ( m_InternalImage->m_Compression == COMPRESSION_OJPEG )
    {
    itkExceptionMacro(<< "This reader cannot read old JPEG compression");
    return;
    }

  // The IO region should be of dimensions 3 otherwise we read only the first
  // page
  if ( m_InternalImage->m_NumberOfPages > 0 && this->GetIORegion().GetImageDimension() > 2 )
    {
    this->ReadVolume(buffer);
    m_InternalImage->Clean();
    return;
    }

  if ( m_InternalImage->m_NumberOfTiles > 0 && this->GetIORegion().GetImageDimension() > 2 )
    {
    this->ReadTiles(buffer);
    m_InternalImage->Clean();
    return;
    }

  int width  = m_InternalImage->m_Width;
  int height = m_InternalImage->m_Height;

  if ( !m_InternalImage->CanRead() )
    {
    uint32 *tempImage;
    tempImage = new uint32[width * height];

    if ( !TIFFReadRGBAImage(m_InternalImage->m_Image,
                            width, height,
                            tempImage, 1) )
      {
      itkExceptionMacro(<< "Cannot read TIFF image or as a TIFF RGBA image");
      if ( tempImage != buffer )
        {
        delete[] tempImage;
        }

      m_InternalImage->Clean();
      return;
      }
    int            xx, yy;
    uint32 *       ssimage;
    unsigned char *fimage = (unsigned char *)buffer;

    for ( yy = 0; yy < height; yy++ )
      {
      ssimage = tempImage + ( height - yy - 1 ) * width;
      for ( xx = 0; xx < width; xx++ )
        {
        unsigned char red   = static_cast< unsigned char >( TIFFGetR(*ssimage) );
        unsigned char green = static_cast< unsigned char >( TIFFGetG(*ssimage) );
        unsigned char blue  = static_cast< unsigned char >( TIFFGetB(*ssimage) );
        unsigned char alpha = static_cast< unsigned char >( TIFFGetA(*ssimage) );

        *( fimage  ) = red;
        *( fimage + 1 ) = green;
        *( fimage + 2 ) = blue;
        *( fimage + 3 ) = alpha;
        fimage += 4;

        ssimage++;
        }
      }

    if ( tempImage != 0 && tempImage != buffer )
      {
      delete[] tempImage;
      }
    m_InternalImage->Clean();
    return;
    }

  unsigned int format = this->GetFormat();

  switch ( format )
    {
    case TIFFImageIO::GRAYSCALE:
    case TIFFImageIO::RGB_:
    case TIFFImageIO::PALETTE_RGB:
    case TIFFImageIO::PALETTE_GRAYSCALE:
      this->ReadGenericImage(buffer, width, height);
      break;
    default:
      return;
    }
  m_InternalImage->Clean();
}

TIFFImageIO::TIFFImageIO()
{
  this->SetNumberOfDimensions(2);
  m_PixelType = SCALAR;
  m_ComponentType = UCHAR;

  this->InitializeColors();
  m_InternalImage = new TIFFReaderInternal;

  m_Spacing[0] = 1.0;
  m_Spacing[1] = 1.0;

  m_Origin[0] = 0.0;
  m_Origin[1] = 0.0;

  m_Compression = TIFFImageIO::PackBits;

  this->AddSupportedWriteExtension(".tif");
  this->AddSupportedWriteExtension(".TIF");
  this->AddSupportedWriteExtension(".tiff");
  this->AddSupportedWriteExtension(".TIFF");

  this->AddSupportedReadExtension(".tif");
  this->AddSupportedReadExtension(".TIF");
  this->AddSupportedReadExtension(".tiff");
  this->AddSupportedReadExtension(".TIFF");
}

TIFFImageIO::~TIFFImageIO()
{
  m_InternalImage->Clean();
  delete m_InternalImage;
}

void TIFFImageIO::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Compression: " << m_Compression << "\n";
}

void TIFFImageIO::InitializeColors()
{
  m_ColorRed    = 0;
  m_ColorGreen  = 0;
  m_ColorBlue   = 0;
  m_TotalColors = -1;
  m_ImageFormat = TIFFImageIO::NOFORMAT;
}

void TIFFImageIO::ReadImageInformation()
{
  // If the internal image was not open we open it.
  // This is usually done when the user sets the ImageIO manually
  if ( !m_InternalImage->m_IsOpen )
    {
    if ( !this->CanReadFile( m_FileName.c_str() ) )
      {
      itkExceptionMacro(<< "Cannot open file " << this->m_FileName << "!");
      return;
      }
    }

  m_Spacing[0] = 1.0;
  m_Spacing[1] = 1.0;

  // If we have some spacing information we use it
  if ( m_InternalImage->m_ResolutionUnit > 0
       && m_InternalImage->m_XResolution > 0
       && m_InternalImage->m_YResolution > 0
       )
    {
    if ( m_InternalImage->m_ResolutionUnit == 2 ) // inches
      {
      m_Spacing[0] = 25.4 / m_InternalImage->m_XResolution;
      m_Spacing[1] = 25.4 / m_InternalImage->m_YResolution;
      }
    else if ( m_InternalImage->m_ResolutionUnit == 3 ) // cm
      {
      m_Spacing[0] = 10.0 / m_InternalImage->m_XResolution;
      m_Spacing[1] = 10.0 / m_InternalImage->m_YResolution;
      }
    }

  m_Origin[0] = 0.0;
  m_Origin[1] = 0.0;

  int width  = m_InternalImage->m_Width;
  int height = m_InternalImage->m_Height;

  m_Dimensions[0] = width;
  m_Dimensions[1] = height;

  switch ( this->GetFormat() )
    {
    case TIFFImageIO::GRAYSCALE:
    case TIFFImageIO::PALETTE_GRAYSCALE:
      this->SetNumberOfComponents(1);
      this->SetPixelType(SCALAR);
      break;
    case TIFFImageIO::RGB_:
      this->SetNumberOfComponents(m_InternalImage->m_SamplesPerPixel);
      this->SetPixelType(RGB);
      break;
    case TIFFImageIO::PALETTE_RGB:
      this->SetNumberOfComponents(3);
      this->SetPixelType(RGB);
      break;
    default:
      this->SetNumberOfComponents(4);
      this->SetPixelType(RGBA);
    }

  if ( !m_InternalImage->CanRead() )
    {
    this->SetNumberOfComponents(4);
    this->SetPixelType(RGBA);
    }

  if ( m_InternalImage->m_BitsPerSample <= 8 )
    {
    if ( m_InternalImage->m_SampleFormat == 2 )
      {
      m_ComponentType = CHAR;
      }
    else
      {
      m_ComponentType = UCHAR;
      }
    }
  else
    {
    if ( m_InternalImage->m_SampleFormat == 2 )
      {
      m_ComponentType = SHORT;
      }
    else
      {
      m_ComponentType = USHORT;
      }
    }

  // We check if we have a Zeiss image.
  // Meaning that the SamplesPerPixel is 2 but the image should be treated as
  // an RGB image.
  if ( m_InternalImage->m_SamplesPerPixel == 2 )
    {
    this->SetNumberOfComponents(3);
    this->SetPixelType(RGB);
    }

  // if the tiff file is multi-pages
  if ( m_InternalImage->m_NumberOfPages - m_InternalImage->m_IgnoredSubFiles > 1 )
    {
    this->SetNumberOfDimensions(3);
    if ( m_InternalImage->m_SubFiles > 0 )
      {
      m_Dimensions[2] = m_InternalImage->m_SubFiles;
      }
    else
      {
      m_Dimensions[2] = m_InternalImage->m_NumberOfPages - m_InternalImage->m_IgnoredSubFiles;
      }
    m_Spacing[2] = 1.0;
    m_Origin[2] = 0.0;
    }

  // if the tiff is tiled
  if ( m_InternalImage->m_NumberOfTiles > 1 )
    {
    this->SetNumberOfDimensions(3);
    m_Dimensions[0] = m_InternalImage->m_TileWidth;
    m_Dimensions[1] = m_InternalImage->m_TileHeight;
    m_Dimensions[2] = m_InternalImage->m_NumberOfTiles;
    m_Spacing[2] = 1.0;
    m_Origin[2] = 0.0;
    }

  return;
}

bool TIFFImageIO::CanWriteFile(const char *name)
{
  std::string filename = name;

  if ( filename == "" )
    {
    return false;
    }

  std::string::size_type TIFFPos = filename.rfind(".TIFF");
  if ( ( TIFFPos != std::string::npos )
       && ( TIFFPos == filename.length() - 5 ) )
    {
    return true;
    }

  TIFFPos = filename.rfind(".tiff");
  if ( ( TIFFPos != std::string::npos )
       && ( TIFFPos == filename.length() - 5 ) )
    {
    return true;
    }

  TIFFPos = filename.rfind(".tif");
  if ( ( TIFFPos != std::string::npos )
       && ( TIFFPos == filename.length() - 4 ) )
    {
    return true;
    }

  TIFFPos = filename.rfind(".TIF");
  if ( ( TIFFPos != std::string::npos )
       && ( TIFFPos == filename.length() - 4 ) )
    {
    return true;
    }

  return false;
}

void TIFFImageIO::WriteImageInformation()
{}

void TIFFImageIO::Write(const void *buffer)
{
  if ( m_NumberOfDimensions == 2 || m_NumberOfDimensions == 3 )
    {
    this->InternalWrite(buffer);
    }
  else
    {
    itkExceptionMacro(<< "TIFF Writer can only write 2-d or 3-d images");
    }
}

class TIFFWriterIO
{
public:
  // Writing file no reading
  static tsize_t TIFFRead(thandle_t, tdata_t, tsize_t) { return 0; }

  // Write data
  static tsize_t TIFFWrite(thandle_t fd, tdata_t buf, tsize_t size)
  {
    std::ostream *out = reinterpret_cast< std::ostream * >( fd );

    out->write(static_cast< char * >( buf ), size);
    return out->fail() ? static_cast< tsize_t >( 0 ) : size;
  }

  static toff_t TIFFSeek(thandle_t fd, toff_t off, int whence)
  {
    std::ostream *out = reinterpret_cast< std::ostream * >( fd );

    switch ( whence )
      {
      case SEEK_SET:
        out->seekp(off, std::ios::beg);
        break;
      case SEEK_END:
        out->seekp(off, std::ios::end);
        break;
      case SEEK_CUR:
        out->seekp(off, std::ios::cur);
        break;
      default:
        return static_cast< toff_t >( out->tellp() );
      }
    return static_cast< toff_t >( out->tellp() );
  }

  static toff_t TIFFSize(thandle_t fd)
  {
    std::ostream *out = reinterpret_cast< std::ostream * >( fd );

    out->seekp(0, std::ios::end);
    return static_cast< toff_t >( out->tellp() );
  }

  static int TIFFMapFile(thandle_t, tdata_t *, toff_t *) { return ( 0 ); }
  static void TIFFUnmapFile(thandle_t, tdata_t, toff_t) {}
};

void TIFFImageIO::InternalWrite(const void *buffer)
{
  char *outPtr = (char *)buffer;

  unsigned int width, height, page, pages = 1;

  width =  m_Dimensions[0];
  height = m_Dimensions[1];
  if ( m_NumberOfDimensions == 3 )
    {
    pages = m_Dimensions[2];
    }

  int    scomponents = this->GetNumberOfComponents();
  float  resolution_x = static_cast< float >( m_Spacing[0] != 0.0 ? 25.4 / m_Spacing[0] : 0.0);
  float  resolution_y = static_cast< float >( m_Spacing[1] != 0.0 ? 25.4 / m_Spacing[1] : 0.0);
  uint32 rowsperstrip = ( uint32 ) - 1;
  int    bps;

  switch ( this->GetComponentType() )
    {
    case UCHAR:
      bps = 8;
      break;
    case CHAR:
      bps = 8;
      break;
    case USHORT:
      bps = 16;
      break;
    case SHORT:
      bps = 16;
      break;
    default:
      itkExceptionMacro(
        << "TIFF supports unsigned/signed char and unsigned/signed short");
    }

  int predictor;

  TIFF *tif = TIFFOpen(m_FileName.c_str(), "w");
  if ( !tif )
    {
    itkExceptionMacro( "Error while trying to open file for writing: "
                       << this->GetFileName()
                       << std::endl
                       << "Reason: "
                       << itksys::SystemTools::GetLastSystemError() );
    }

  if ( this->GetComponentType() == SHORT
       || this->GetComponentType() == CHAR )
    {
    TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_INT);
    }

  uint32 w = width;
  uint32 h = height;

  if ( m_NumberOfDimensions == 3 )
    {
    TIFFCreateDirectory(tif);
    }
  for ( page = 0; page < pages; page++ )
    {
    TIFFSetDirectory(tif, page);
    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, w);
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, h);
    TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, scomponents);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, bps); // Fix for stype
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    if ( this->GetComponentType() == SHORT
         || this->GetComponentType() == CHAR )
      {
      TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_INT);
      }
    TIFFSetField(tif, TIFFTAG_SOFTWARE, "InsightToolkit");

    if ( scomponents > 3 )
      {
      // if number of scalar components is greater than 3, that means we assume
      // there is alpha.
      uint16  extra_samples = scomponents - 3;
      uint16 *sample_info = new uint16[scomponents - 3];
      sample_info[0] = EXTRASAMPLE_ASSOCALPHA;
      int cc;
      for ( cc = 1; cc < scomponents - 3; cc++ )
        {
        sample_info[cc] = EXTRASAMPLE_UNSPECIFIED;
        }
      TIFFSetField(tif, TIFFTAG_EXTRASAMPLES, extra_samples,
                   sample_info);
      delete[] sample_info;
      }

    int compression;

    if ( m_UseCompression )
      {
      switch ( m_Compression )
        {
        case TIFFImageIO::PackBits:
          compression = COMPRESSION_PACKBITS; break;
        case TIFFImageIO::JPEG:
          compression = COMPRESSION_JPEG; break;
        case TIFFImageIO::Deflate:
          compression = COMPRESSION_DEFLATE; break;
        case TIFFImageIO::LZW:
          compression = COMPRESSION_LZW; break;
        default:
          compression = COMPRESSION_NONE;
        }
      }
    else
      {
      compression = COMPRESSION_NONE;
      }

    TIFFSetField(tif, TIFFTAG_COMPRESSION, compression); // Fix for compression

    uint16 photometric = ( scomponents == 1 ) ? PHOTOMETRIC_MINISBLACK : PHOTOMETRIC_RGB;

    if ( compression == COMPRESSION_JPEG )
      {
      TIFFSetField(tif, TIFFTAG_JPEGQUALITY, 75); // Parameter
      TIFFSetField(tif, TIFFTAG_JPEGCOLORMODE, JPEGCOLORMODE_RGB);
      photometric = PHOTOMETRIC_YCBCR;
      }
    else if ( compression == COMPRESSION_LZW )
      {
      predictor = 2;
      TIFFSetField(tif, TIFFTAG_PREDICTOR, predictor);
      itkDebugMacro(<< "LZW compression is patented outside US so it is disabled");
      }
    else if ( compression == COMPRESSION_DEFLATE )
      {
      predictor = 2;
      TIFFSetField(tif, TIFFTAG_PREDICTOR, predictor);
      }

    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, photometric); // Fix for scomponents

    TIFFSetField( tif,
                  TIFFTAG_ROWSPERSTRIP,
                  TIFFDefaultStripSize(tif, rowsperstrip) );

    if ( resolution_x > 0 && resolution_y > 0 )
     {
     TIFFSetField(tif, TIFFTAG_XRESOLUTION, resolution_x);
     TIFFSetField(tif, TIFFTAG_YRESOLUTION, resolution_y);
     TIFFSetField(tif, TIFFTAG_RESOLUTIONUNIT, RESUNIT_INCH);
     }

    if ( m_NumberOfDimensions == 3 )
      {
      // We are writing single page of the multipage file
      TIFFSetField(tif, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
      // Set the page number
      TIFFSetField(tif, TIFFTAG_PAGENUMBER, page, pages);
      }
    int rowLength; // in bytes

    switch ( this->GetComponentType() )
      {
      case UCHAR:
        rowLength = sizeof( unsigned char );
        break;
      case USHORT:
        rowLength = sizeof( unsigned short );
        break;
      case CHAR:
        rowLength = sizeof( char );
        break;
      case SHORT:
        rowLength = sizeof( short );
        break;
      default:
        itkExceptionMacro(
          << "TIFF supports unsigned/signed char and unsigned/signed short");
      }

    rowLength *= this->GetNumberOfComponents();
    rowLength *= width;

    int row = 0;
    for ( unsigned int idx2 = 0; idx2 < height; idx2++ )
      {
      if ( TIFFWriteScanline(tif, const_cast< char * >( outPtr ), row, 0) < 0 )
        {
        itkExceptionMacro(<< "TIFFImageIO: error out of disk space");
        break;
        }
      outPtr += rowLength;
      row++;
      }

    if ( m_NumberOfDimensions == 3 )
      {
      TIFFWriteDirectory(tif);
      }
    }
  TIFFClose(tif);
}

bool TIFFImageIO::CanFindTIFFTag(unsigned int t)
{
  // m_InternalImage needs to be valid
  if ( !m_InternalImage )
    {
    itkExceptionMacro(<< "Need to call CanReadFile before");
    return false;
    }

  ttag_t               tag = t; // 32bits integer
  const TIFFFieldInfo *fld = TIFFFieldWithTag(m_InternalImage->m_Image, tag);
  if ( fld == NULL )
    {
    return false;
    }
  return true;
}

void * TIFFImageIO::ReadRawByteFromTag(unsigned int t, short & value_count)
{
  // m_InternalImage needs to be valid
  if ( !m_InternalImage )
    {
    itkExceptionMacro(<< "Need to call CanReadFile before");
    return NULL;
    }
  ttag_t               tag = t;
  void *               raw_data = NULL;
  const TIFFFieldInfo *fld = TIFFFieldWithTag(m_InternalImage->m_Image, tag);
  if ( fld == NULL )
    {
    itkExceptionMacro(<< "fld is NULL");
    return NULL;
    }
  else
    {
    if ( fld->field_passcount )
      {
      if ( TIFFGetField(m_InternalImage->m_Image, tag, &value_count, &raw_data) != 1 )
        {
        itkExceptionMacro(<< "Tag cannot be found");
        return NULL;
        }
      else
        {
        if ( fld->field_type != TIFF_BYTE )
          {
          itkExceptionMacro(<< "Tag is not of type TIFF_BYTE");
          return NULL;
          }
        }
      }
    }
  return raw_data;
}
} // end namespace itk
