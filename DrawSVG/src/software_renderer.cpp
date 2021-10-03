#include "software_renderer.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "triangulation.h"

using namespace std;

namespace CMU462 {


// Implements SoftwareRenderer //

void SoftwareRendererImp::draw_svg( SVG& svg ) {

  // set top level transformation
  transformation = svg_2_screen;

  // draw all elements
  for ( size_t i = 0; i < svg.elements.size(); ++i ) {
    draw_element(svg.elements[i]);
  }

  // draw canvas outline
  Vector2D a = transform(Vector2D(    0    ,     0    )); a.x--; a.y--;
  Vector2D b = transform(Vector2D(svg.width,     0    )); b.x++; b.y--;
  Vector2D c = transform(Vector2D(    0    ,svg.height)); c.x--; c.y++;
  Vector2D d = transform(Vector2D(svg.width,svg.height)); d.x++; d.y++;

  rasterize_line(a.x, a.y, b.x, b.y, Color::Black);
  rasterize_line(a.x, a.y, c.x, c.y, Color::Black);
  rasterize_line(d.x, d.y, b.x, b.y, Color::Black);
  rasterize_line(d.x, d.y, c.x, c.y, Color::Black);

  // resolve and send to render target
  resolve();

}

void SoftwareRendererImp::set_sample_rate( size_t sample_rate ) {

  // Task 4: 
  // You may want to modify this for supersampling support
  this->sample_rate = sample_rate;
  this->sample_w = this->target_w * this->sample_rate;
  this->sample_h = this->target_h * this->sample_rate;
  this->sample_buffer.resize(4 * this->sample_w * this->sample_h);

}

void SoftwareRendererImp::set_render_target( unsigned char* render_target,
                                             size_t width, size_t height ) {

  // Task 4: 
  // You may want to modify this for supersampling support
  this->render_target = render_target;
  this->target_w = width;
  this->target_h = height;
  this->sample_w = target_w * this->sample_rate;
  this->sample_h = target_h * this->sample_rate;
  this->sample_buffer.resize(4 * this->sample_w * this->sample_h);
}

void SoftwareRendererImp::draw_element( SVGElement* element ) {

  // Task 5 (part 1):
  // Modify this to implement the transformation stack

  // push transformation matrix
  Matrix3x3 transform_save = transformation;

  // set object transformation
  transformation = transformation * element->transform;

  switch(element->type) {
    case POINT:
      draw_point(static_cast<Point&>(*element));
      break;
    case LINE:
      draw_line(static_cast<Line&>(*element));
      break;
    case POLYLINE:
      draw_polyline(static_cast<Polyline&>(*element));
      break;
    case RECT:
      draw_rect(static_cast<Rect&>(*element));
      break;
    case POLYGON:
      draw_polygon(static_cast<Polygon&>(*element));
      break;
    case ELLIPSE:
      draw_ellipse(static_cast<Ellipse&>(*element));
      break;
    case IMAGE:
      draw_image(static_cast<Image&>(*element));
      break;
    case GROUP:
      draw_group(static_cast<Group&>(*element));
      break;
    default:
      break;
  }
  
  transformation = transform_save;
}


// Primitive Drawing //

void SoftwareRendererImp::draw_point( Point& point ) {

  Vector2D p = transform(point.position);
  rasterize_point( p.x, p.y, point.style.fillColor );

}

void SoftwareRendererImp::draw_line( Line& line ) { 

  Vector2D p0 = transform(line.from);
  Vector2D p1 = transform(line.to);
  rasterize_line( p0.x, p0.y, p1.x, p1.y, line.style.strokeColor );

}

void SoftwareRendererImp::draw_polyline( Polyline& polyline ) {

  Color c = polyline.style.strokeColor;

  if( c.a != 0 ) {
    int nPoints = polyline.points.size();
    for( int i = 0; i < nPoints - 1; i++ ) {
      Vector2D p0 = transform(polyline.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polyline.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_rect( Rect& rect ) {

  Color c;
  
  // draw as two triangles
  float x = rect.position.x;
  float y = rect.position.y;
  float w = rect.dimension.x;
  float h = rect.dimension.y;

  Vector2D p0 = transform(Vector2D(   x   ,   y   ));
  Vector2D p1 = transform(Vector2D( x + w ,   y   ));
  Vector2D p2 = transform(Vector2D(   x   , y + h ));
  Vector2D p3 = transform(Vector2D( x + w , y + h ));
  
  // draw fill
  c = rect.style.fillColor;
  if (c.a != 0 ) {
    rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    rasterize_triangle( p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c );
  }

  // draw outline
  c = rect.style.strokeColor;
  if( c.a != 0 ) {
    rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    rasterize_line( p1.x, p1.y, p3.x, p3.y, c );
    rasterize_line( p3.x, p3.y, p2.x, p2.y, c );
    rasterize_line( p2.x, p2.y, p0.x, p0.y, c );
  }

}

void SoftwareRendererImp::draw_polygon( Polygon& polygon ) {

  Color c;

  // draw fill
  c = polygon.style.fillColor;
  if( c.a != 0 ) {

    // triangulate
    vector<Vector2D> triangles;
    triangulate( polygon, triangles );

    // draw as triangles
    for (size_t i = 0; i < triangles.size(); i += 3) {
      Vector2D p0 = transform(triangles[i + 0]);
      Vector2D p1 = transform(triangles[i + 1]);
      Vector2D p2 = transform(triangles[i + 2]);
      rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    }
  }

  // draw outline
  c = polygon.style.strokeColor;
  if( c.a != 0 ) {
    int nPoints = polygon.points.size();
    for( int i = 0; i < nPoints; i++ ) {
      Vector2D p0 = transform(polygon.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polygon.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_ellipse( Ellipse& ellipse ) {

  // Extra credit 

}

void SoftwareRendererImp::draw_image( Image& image ) {

  Vector2D p0 = transform(image.position);
  Vector2D p1 = transform(image.position + image.dimension);

  rasterize_image( p0.x, p0.y, p1.x, p1.y, image.tex );
}

void SoftwareRendererImp::draw_group( Group& group ) {

  for ( size_t i = 0; i < group.elements.size(); ++i ) {
    draw_element(group.elements[i]);
  }

}

// Rasterization //

// The input arguments in the rasterization functions 
// below are all defined in screen space coordinates

void SoftwareRendererImp::rasterize_point( float x, float y, Color color ) {

  // fill in the nearest pixel
  int sx = (int) floor(x);
  int sy = (int) floor(y);

  // check bounds
  if ( sx < 0 || sx >= target_w ) return;
  if ( sy < 0 || sy >= target_h ) return;

  // fill sample - NOT doing alpha blending!
  // render_target[4 * (sx + sy * target_w)    ] = (uint8_t) (color.r * 255);
  // render_target[4 * (sx + sy * target_w) + 1] = (uint8_t) (color.g * 255);
  // render_target[4 * (sx + sy * target_w) + 2] = (uint8_t) (color.b * 255);
  // render_target[4 * (sx + sy * target_w) + 3] = (uint8_t) (color.a * 255);


  sx *= sample_rate;
  sy *= sample_rate;

  for (int i = sx; i < sx + sample_rate; ++i)
  {
    for (int j = sy; j < sy + sample_rate; ++j)
    {
      size_t pxy = 4 * (i + j * sample_w);
      sample_buffer[pxy] = (uint8_t)(color.r * 255);
      sample_buffer[pxy + 1] = (uint8_t)(color.g * 255);
      sample_buffer[pxy + 2] = (uint8_t)(color.b * 255);
      sample_buffer[pxy + 3] = (uint8_t)(color.a * 255);
    }
  }

  fill_sample(sx,sy,color);

}

void SoftwareRendererImp::rasterize_line( float x0, float y0,
                                          float x1, float y1,
                                          Color color) {

  // Task 2: 
  // Implement line rasterization
  x0 *= sample_rate;
  x1 *= sample_rate;
  y0 *= sample_rate;
  y1 *= sample_rate;
  // The Bresenham Line-Drawing Algorithm
  float m = (y1 - y0) / (x1 - x0);
  float e = 0;
  float x,y;
  float ystep = 1;
  float xstep = 1;
 
  if(-1 < m && m < 1) {  // if the slope is between (-1,1)
    if(x1 < x0) {
      float temp = x0;
      x0 = x1;
      x1 = temp;
      temp = y0;
      y0 = y1;
      y1 = temp;
    }
    y = y0;
    for(x = x0; x < x1; x++) {
      rasterize_point(x, y, color);
      if(y1>y0) { // if slope is positive
        if(e + m < 0.5) {
          e = e + m;
        } 
        else {
          y = y + ystep;
          e = e + m - ystep;
        }
      } 
      else {
        if(-0.5 < e + m) {
          e = e + m;
        } 
        else {
          y = y + ystep * -1; // if slope is negative, let ystep be -1
          e = e + m - ystep * -1;
        }
      }
    }
  } 
  else { // if the slope is outside (-1,1)
    m = (x1 - x0) / (y1 - y0);
    if(y1 < y0) {
      float temp = x0;
      x0 = x1;
      x1 = temp;
      temp = y0;
      y0 = y1;
      y1 = temp;
    }
    x = x0;
    for(y = y0; y < y1; y++) {
      rasterize_point(x, y, color);
      if(x1 > x0) { // if slope is positive
        if(e + m < 0.5) {
          e = e + m;
        } 
        else {
          x = x + xstep;
          e = e + m - xstep;
        }
      } 
      else {
        if(-0.5 < e + m) {
          e = e + m;
        } 
        else {
          x = x + xstep * -1;
          e = e + m - xstep * -1;
        }
      }
    }
  }
}

void SoftwareRendererImp::rasterize_triangle( float x0, float y0,
                                              float x1, float y1,
                                              float x2, float y2,
                                              Color color ) {
  // Task 3: 
  // Implement triangle rasterization 
  x0 *= sample_rate;
  y0 *= sample_rate;
  x1 *= sample_rate;
  y1 *= sample_rate;
  x2 *= sample_rate;
  y2 *= sample_rate;
  // draw edges
  rasterize_line(x0, y0, x1, y1, color);
  rasterize_line(x1, y1, x2, y2, color);
  rasterize_line(x2, y2, x0, y0, color);

  // The Barycentric Algorithm
  // fill the triangle
  float xMin = min(x0,min(x1,x2));
  float yMin = min(y0,min(y1,y2));
  float xHigh = max(x0,max(x1,x2));
  float yHigh = max(y0,max(y1,y2));

  float x, y;
  for(x = xMin; x < xHigh; x += 0.5) {
    for(y = yMin; y < yHigh; y += 0.5) {

      float area_half1=(x1-x0)*(y-y0)-(y1-y0)*(x-x0);
      float area_half2=(x-x0)*(y2-y0)-(y-y0)*(x2-x0);
      float area_whole=(x1-x0)*(y2-y0)-(y1-y0)*(x2-x0);
      float a1 = area_half1/area_whole;
      float a2 = area_half2/area_whole;
      float a3 = a1+a2;
      // don't need to consider area1/area == 0, area2/area == 0 situation cause already draw the line. 
      if((a1>0 || a1==0) && (a2>0 || a2==0) && (a3<1 || a3==1)){
        rasterize_point(x, y, color);
      }
    }
  }

}

void SoftwareRendererImp::rasterize_image( float x0, float y0,
                                           float x1, float y1,
                                           Texture& tex ) {
  // Task 6: 
  // Implement image rasterization
  x0 *= sample_rate;
  x1 *= sample_rate;
  y0 *= sample_rate;
  y1 *= sample_rate;
  Sampler2D* sampler = this->sampler;

  float u_step = 1 / (x1 - x0);
  float v_step = 1 / (y1 - y0);
  Color c;
  for(float x = x0; x < x1; x++) {
    for(float y = y0; y < y1; y++) {
      float u = (x - x0) / (y1 - y0);
      float v = (y - y0) / (y1 - y0);

      // c = sampler->sample_nearest(tex, u, v, 0);
      c = sampler->sample_bilinear(tex, u, v, 0);
      // c = sampler->sample_trilinear(tex, u, v, u_step, v_step);
      rasterize_point(x, y, c);
    }
  }
}

// resolve samples to render target
void SoftwareRendererImp::resolve( void ) {

  // Task 4: 
  // Implement supersampling
  // You may also need to modify other functions marked with "Task 4".
  for(size_t x=0; x< sample_w; x+=sample_rate)
  {
    for (size_t y = 0; y < sample_h; y += sample_rate)
    {
      uint16_t r = 0, g = 0, b = 0, a = 0;
      for (size_t i = 0; i < sample_rate; ++i)
      {
        for (size_t j = 0; j < sample_rate; ++j)
        {
          r += sample_buffer[4 * (x + i + (y + j) * sample_w)    ];
          g += sample_buffer[4 * (x + i + (y + j) * sample_w) + 1];
          b += sample_buffer[4 * (x + i + (y + j) * sample_w) + 2];
          a += sample_buffer[4 * (x + i + (y + j) * sample_w) + 3];
        }
      }

      r /= sample_rate*sample_rate;
      g /= sample_rate*sample_rate;
      b /= sample_rate*sample_rate;
      a /= sample_rate*sample_rate;

      size_t sx = x / sample_rate;
      size_t sy = y / sample_rate;

      render_target[4 * (sx + sy * target_w)    ] = (uint8_t)(r);
      render_target[4 * (sx + sy * target_w) + 1] = (uint8_t)(g);
      render_target[4 * (sx + sy * target_w) + 2] = (uint8_t)(b);
      render_target[4 * (sx + sy * target_w) + 3] = (uint8_t)(a);
    }
  }

  return;

}

void SoftwareRendererImp::fill_sample(int sx, int sy, const Color &color)
{
  // check bounds
  if (sx < 0 || sx >= this->sample_w)
    return;
  if (sy < 0 || sy >= this->sample_h)
    return;

  // task 8
  // Simple alpha compositing - premultipled method
  float Er = color.r, Eg = color.g, Eb = color.b, Ea = color.a;
  float Cr = sample_buffer[4 * (sx + sy * sample_w)    ]/255.0f;
  float Cg = sample_buffer[4 * (sx + sy * sample_w)  +1]/255.0f;
  float Cb = sample_buffer[4 * (sx + sy * sample_w)  +2]/255.0f;
  float Ca = sample_buffer[4 * (sx + sy * sample_w)  +3]/255.0f;

  float Ca_new = 1.0f - (1.0f - Ea) * (1.0f - Ca);
  float Cr_new = (1.0f - Ea) * Cr + Er;
  float Cg_new = (1.0f - Ea) * Cg + Eg;
  float Cb_new = (1.0f - Ea) * Cb + Eb;

  sample_buffer[4 * (sx + sy * sample_w)    ] = (uint8_t) (Cr_new * 255);
  sample_buffer[4 * (sx + sy * sample_w) + 1] = (uint8_t) (Cg_new * 255);
  sample_buffer[4 * (sx + sy * sample_w) + 2] = (uint8_t) (Cb_new * 255);
  sample_buffer[4 * (sx + sy * sample_w) + 3] = (uint8_t) (Ca_new * 255);
}

} // namespace CMU462
