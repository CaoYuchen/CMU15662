#include "texture.h"
#include "color.h"

#include <assert.h>
#include <iostream>
#include <algorithm>

using namespace std;

namespace CMU462 {

inline void uint8_to_float( float dst[4], unsigned char* src ) {
  uint8_t* src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
  dst[3] = src_uint8[3] / 255.f;
}

inline void float_to_uint8( unsigned char* dst, float src[4] ) {
  uint8_t* dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[0])));
  dst_uint8[1] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[1])));
  dst_uint8[2] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[2])));
  dst_uint8[3] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[3])));
}

void Sampler2DImp::generate_mips(Texture& tex, int startLevel) {

  // NOTE: 
  // This starter code allocates the mip levels and generates a level 
  // map by filling each level with placeholder data in the form of a 
  // color that differs from its neighbours'. You should instead fill
  // with the correct data!

  // Task 7: Implement this

  // check start level
  if ( startLevel >= tex.mipmap.size() ) {
    std::cerr << "Invalid start level"; 
  }

  // allocate sublevels
  int baseWidth  = tex.mipmap[startLevel].width;
  int baseHeight = tex.mipmap[startLevel].height;
  int numSubLevels = (int)(log2f( (float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  tex.mipmap.resize(startLevel + numSubLevels + 1);

  int width  = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel& level = tex.mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width  = max( 1, width  / 2); assert(width  > 0);
    height = max( 1, height / 2); assert(height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(4 * width * height);

  }

  // fill all 0 sub levels with interchanging colors (JUST AS A PLACEHOLDER)
  // Color colors[3] = { Color(1,0,0,1), Color(0,1,0,1), Color(0,0,1,1) };
  // for(size_t i = 1; i < tex.mipmap.size(); ++i) {

  //   Color c = colors[i % 3];
  //   MipLevel& mip = tex.mipmap[i];

  //   for(size_t i = 0; i < 4 * mip.width * mip.height; i += 4) {
  //     float_to_uint8( &mip.texels[i], &c.r );
  //   }
  // }

  for(size_t index = 1; index < tex.mipmap.size(); ++index) {
    MipLevel& preMip = tex.mipmap[index - 1];
    MipLevel& curMip = tex.mipmap[index];

    for (size_t x = 0; x < curMip.width; ++x) {
      size_t preX = x * 2;
      for (size_t y = 0;  y < curMip.height; ++y) {
        size_t preY = y * 2;
        uint16_t r = 0, g = 0, b = 0, a = 0;
        for (int i = 0; i < 2; ++i) {
          for (int j = 0; j < 2; ++j) {
            r += preMip.texels[4 * (preX + i + (preY + j) * preMip.width)];
            g += preMip.texels[4 * (preX + i + (preY + j) * preMip.width) + 1];
            b += preMip.texels[4 * (preX + i + (preY + j) * preMip.width) + 2];
            a += preMip.texels[4 * (preX + i + (preY + j) * preMip.width) + 3];
          }
        }
        curMip.texels[4 * (x + y * curMip.width)] =  r / 4;
        curMip.texels[4 * (x + y * curMip.width) + 1] =  g / 4;
        curMip.texels[4 * (x + y * curMip.width) + 2] =  b / 4;
        curMip.texels[4 * (x + y * curMip.width) + 3] =  a / 4;
      }
    }
  }

}

Color Sampler2DImp::sample_nearest(Texture& tex, 
                                   float u, float v, 
                                   int level) {

  // Task 6: Implement nearest neighbour interpolation
  if(level > kMaxMipLevels && level < 0)
  {
    return Color(1,0,1,1);
  }
  size_t x = round(u * tex.mipmap[level].width - 0.5f);
  size_t y = round(v * tex.mipmap[level].height - 0.5f);

  Color c;
  uint8_to_float(&c.r, &tex.mipmap[level].texels[4*(x+y*tex.mipmap[level].width)]);

  // return magenta for invalid level
  return c;
}

Color Sampler2DImp::sample_bilinear(Texture& tex, 
                                    float u, float v, 
                                    int level) {
  
  // Task 6: Implement bilinear filtering
  if(level > kMaxMipLevels && level < 0)
  {
    return Color(1,0,1,1);
  }

  MipLevel& mip = tex.mipmap[level];
  size_t x = abs(floor(u * mip.width - 0.5f));
  size_t y = abs(floor(v * mip.height - 0.5f));

  float u_ratio = u * mip.width - 0.5f - x;
  float v_ratio = v * mip.height - 0.5f - y;

  Color c1, c2, c3, c4;
  uint8_to_float(&c1.r, &mip.texels[4 * (x + mip.width * y)]);
  uint8_to_float(&c2.r, &mip.texels[4 * (x + 1 + mip.width * y)]);
  uint8_to_float(&c3.r, &mip.texels[4 * (x + mip.width * (y + 1))]);
  uint8_to_float(&c4.r, &mip.texels[4 * (x + 1 + mip.width * (y + 1))]);
  // return magenta for invalid level

  Color c = (c1 * (1 - u_ratio) + c2 * u_ratio) * (1 - v_ratio) +
            (c3 * (1 - u_ratio) + c4 * u_ratio) * v_ratio;
  return c;

}

Color Sampler2DImp::sample_trilinear(Texture& tex, 
                                     float u, float v, 
                                     float u_scale, float v_scale) {

  // Task 7: Implement trilinear filtering
  float d = max(0.f, log2f(max(u_scale, v_scale)));
  int k0 = floor(d), k1 = k0 + 1;
  float a = k1 - d, b = 1 - a;

  Color c0 = sample_bilinear(tex, u, v, k0);
  Color c1 = sample_bilinear(tex, u, v, k1);

  return a * c0 + b * c1;
  // return magenta for invalid level

}

} // namespace CMU462
