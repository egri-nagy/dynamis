#include "png_util.h"

#include <stdlib.h>

#define ERROR -1
#define OK 0

#define BIT_DEPTH 8

/* write an indexed 8bpp png file */
int write_indexed_png(char *file_name, png_byte *image, png_uint_32 width, png_uint_32 height, png_byte *palette, png_byte palette_length){
  int i;
  FILE *fp;
  png_structp png_ptr;
  png_infop info_ptr;
  png_colorp png_palette;
  png_textp text_ptr;
  png_uint_32 k;
  png_bytep *row_pointers;

  /* open the file */
  fp = fopen(file_name, "wb");
  if (fp == NULL){
    return (ERROR);
  }
  
  /* Create and initialize the png_struct with the desired error handler
   * functions. 
   */
  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  
  if (png_ptr == NULL){
    fclose(fp);
    return (ERROR);
  }
  
  /* Allocate/initialize the image information data. */
  info_ptr = png_create_info_struct(png_ptr);
  if (info_ptr == NULL){
    fclose(fp);
    png_destroy_write_struct(&png_ptr,  (png_infopp)NULL);
    return (ERROR);
  }
  
   /* Set error handling.    */
  if (setjmp(png_jmpbuf(png_ptr))){
    /* If we get here, we had a problem reading the file */
    fclose(fp);
    png_destroy_write_struct(&png_ptr, &info_ptr);
    return (ERROR);
  }
  
   /* set up the output control if you are using standard C streams */
   png_init_io(png_ptr, fp);

   /* the image header */
   png_set_IHDR(png_ptr, info_ptr, width, height, BIT_DEPTH, PNG_COLOR_TYPE_PALETTE,
      PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

   /* set the palette if there is one.  REQUIRED for indexed-color images */
   png_palette = (png_colorp)png_malloc(png_ptr, palette_length * sizeof (png_color));
   /* ... set palette colors ... */
   for (i = 0; i < palette_length; i++){
     png_palette[i].red = palette[i*3];
     png_palette[i].green = palette[i*3+1];
     png_palette[i].blue = palette[i*3+2];
   }
   png_set_PLTE(png_ptr, info_ptr, png_palette, palette_length);

   
   /* Optionally write comments into the image */
   text_ptr = (png_textp) calloc(3, sizeof(png_text));
  if (text_ptr == NULL){
    fclose(fp);
    return (ERROR);
  }

  
   text_ptr[0].key = "Title";
   text_ptr[0].text = "Dynamis Output";
   text_ptr[0].compression = PNG_TEXT_COMPRESSION_NONE;
   text_ptr[1].key = "Author";
   text_ptr[1].text = "sirna";
   text_ptr[1].compression = PNG_TEXT_COMPRESSION_NONE;
   text_ptr[2].key = "Description";
   text_ptr[2].text = "Generated output picture from dynamis.";
   text_ptr[2].compression = PNG_TEXT_COMPRESSION_zTXt;
   
   png_set_text(png_ptr, info_ptr, text_ptr, 3);

   /* Write the file header information. */
   png_write_info(png_ptr, info_ptr);

   //creating scanlines
   row_pointers = (png_bytep*) png_malloc(png_ptr, height * sizeof(png_bytep*));
   for (k = 0; k < height; k++){
     row_pointers[k] = image + k*width;
   }
     
   png_write_image(png_ptr, row_pointers);
   
   /* It is REQUIRED to call this to finish writing the rest of the file */
   png_write_end(png_ptr, info_ptr);
   

   png_free(png_ptr, png_palette);
   palette=NULL;

   /* clean up after the write, and free any memory allocated */
   png_destroy_write_struct(&png_ptr, &info_ptr);

   free(text_ptr);

   /* close the file */
   fclose(fp);

   /* that's it */
   return (OK);
}

/* //SIMPLE EXAMPLE OF USAGE: RED-BLUE PATTERN 
  int main(){
  int i,j;
  png_byte pal[6] = {255,0,0,0,0,255};
  png_byte *image = (png_byte*) calloc(10*10, sizeof(png_byte));
  for (i = 0; i < 10; i++){
  for (j = 0; j < 10; j++){
  image[10*j+i] = (i+j)%2;
  }
  } 
  write_indexed_png("x.png",image, 10, 10, pal, 2);
  }
*/
