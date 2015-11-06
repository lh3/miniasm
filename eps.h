#ifndef EPS_H_
#define EPS_H_

#include <stdio.h>

#define EPS FILE
#define EPSPTR  FILE *
#define eps_open(s) fopen((s),"w+")
#define eps_close(fp) fclose(fp)

#define eps_header(fp,x,y,linewidth) { \
	fprintf(fp,"%%!PS-Adobe-3.0 EPSF-3.0\n"); \
	fprintf(fp,"%%%%BoundingBox:"); \
	fprintf(fp," 1 1 %g %g\n\n",(float)(x),(float)(y)); \
	fprintf(fp,"/C { dup 255 and 255 div exch dup -8 bitshift 255 and 255 div 3 1 roll -16 bitshift 255 and 255 div 3 1 roll setrgbcolor } bind def\n"); \
	fprintf(fp,"/L { 4 2 roll moveto lineto } bind def\n"); \
	fprintf(fp,"/LX { dup 4 -1 roll exch moveto lineto } bind def\n"); \
	fprintf(fp,"/LY { dup 4 -1 roll moveto exch lineto } bind def\n"); \
	fprintf(fp,"/LS { 3 1 roll moveto show } bind def\n"); \
	fprintf(fp,"/MS { dup stringwidth pop 2 div 4 -1 roll exch sub 3 -1 roll moveto show } bind def\n"); \
	fprintf(fp,"/RS { dup stringwidth pop 4 -1 roll exch sub 3 -1 roll moveto show } bind def\n"); \
	fprintf(fp,"/B { 4 copy 3 1 roll exch 6 2 roll 8 -2 roll moveto lineto lineto lineto closepath } bind def\n");\
	fprintf(fp,"%g setlinewidth\n\n",linewidth);\
}
#define eps_font(fp,f,s) do { \
	fprintf(fp,"/FS %d def\n",s); \
	fprintf(fp,"/FS4 FS 4 div def\n"); \
	fprintf(fp,"/%s findfont FS scalefont setfont\n\n",f); \
  } while (0)

#define eps_bottom(fp) fprintf(fp,"stroke showpage\n")
#define eps_color(fp,col) fprintf(fp,"stroke %d C\n",col)
#define eps_gray(fp,gray) fprintf(fp, "%g setgray\n",(float)gray)
#define eps_linewidth(fp, lw) fprintf(fp, "%g setlinewidth\n", (float)(lw))
#define eps_line(fp,x1,y1,x2,y2) fprintf(fp,"%g %g %g %g L\n",(float)(x1),(float)(y1),(float)(x2),(float)(y2))
#define eps_linex(fp,x1,x2,y) fprintf(fp,"%g %g %g LX\n",(float)(x1),(float)(x2),(float)(y))
#define eps_liney(fp,y1,y2,x) fprintf(fp,"%g %g %g LY\n",(float)(y1),(float)(y2),(float)(x))
#define eps_Lstr(fp,x,y,s) fprintf(fp,"%g %g (%s) LS\n",(float)(x),(float)(y),s)
#define eps_Mstr(fp,x,y,s) fprintf(fp,"%g %g (%s) MS\n",(float)(x),(float)(y),s)
#define eps_Rstr(fp,x,y,s) fprintf(fp,"%g %g (%s) RS\n",(float)(x),(float)(y),s)
#define eps_Lstr4(fp,x,y,s) fprintf(fp,"%g %g FS4 add (%s) LS\n",(float)(x),(float)(y),s)
#define eps_Rstr4(fp,x,y,s) fprintf(fp,"%g %g FS4 add (%s) RS\n",(float)(x),(float)(y),s)
#define eps_Lstr4s(fp,x,y,s) fprintf(fp,"%g %g FS4 sub (%s) LS\n",(float)(x),(float)(y),s)
#define eps_Rstr4s(fp,x,y,s) fprintf(fp,"%g %g FS4 sub (%s) RS\n",(float)(x),(float)(y),s)
#define eps_box(fp,x1,y1,x2,y2) fprintf(fp,"%g %g %g %g B\n",(float)(x1),(float)(y1),(float)(x2),(float)(y2))
#define eps_fill(fp) fprintf(fp,"fill\n")
#define eps_stroke(fp) fprintf(fp,"stroke\n")

#endif
