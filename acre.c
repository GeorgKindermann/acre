/*  acre 0.1a, Programm to show the influence area of points and their density
    Copyright (C) 2006  Georg Kindermann

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

Post: Georg Kindermann, Bundesforschungszentrum für Wald, Institut für Waldwachstum und Waldbau
      Seckendorff-Gudent-Weg 8, 1131 Wien, Austria
Home: https://github.comGeorgKindermann/acre */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

unsigned long *nr; /*Point Nr.*/
int *ba; /*Pointtype*/
double *w; /*Pointimportance*/
double *x,*y; /*x and y Koordiantes*/
double *max; /*Max influence radius*/
unsigned long n,nn; /*Number of points(n) and arraysize(nn)*/

float *area; /*standing area of tree*/
unsigned long *ysort; /*y influence sorted trees index*/
unsigned long nysort; /*number of first unused tree in ysort*/
unsigned long *stab; /*selected tab with y influence for the x-row*/
unsigned long nstab; /*number of stab*/
unsigned long *atab; /*selected stab with x influence*/
unsigned long natab; /*number of atab*/
char *bord; /*Register if Tree is at the border or not*/

double xmin=0,xmax=0,ymin=0,ymax=0;
/*the outer corners of the area (-x nr nr nr nr)*/
int res=400,xres,yres; /*the reslution of the picture (-r <nr>)*/
float RES=0; /*the pixelsize (-R <nr>)*/
double maa=-1.7E+308,mia=1.7E+308;  /*mia..minimum (maa..maximum) activation of the pic*/
float *pic; /*the register for the calculated picturedata*/
unsigned long *picnr; /*tree Nr.*/

char iname[255]; /*name of the inputfile*/
char oname[251]; /*name of the outputfile*/
char aname[255]; /*for outpuname and extension*/
char mdl=2; /*modell number (-m nr)*/
unsigned char selec = 95; /*selection for the saving (-s <nr> only 4 bit)
			      1 textoutput 1
			      2 sw-pic 1
			      4 grec pic 1
			      8 color pic 1
			     16 border 1
			     32 text to stdout 0
			     64 autoclac min/max x/y 1
			    128 autocalc influecradius 0*/
unsigned char sout = 5; /*selection for the textoutput (-t <nr> only 4 bit)
			      1 pixel(0)/area(1) DEF 1
			      2 drop bordertrees 0
			      4 show header 1
			      8 replace area with NA if bordertree 0
			     16 unused
			     32 unused
			     64 unused
			     128 unused*/
double infl = 0; /*influence factor*/

FILE *fp;

int main (int argc, char **argv)
{
  int i,reg;
  void dload (void);
  void mpic (int mdl);
  void gsave (void);
  void csave (void);
  void bsave (void);
  void sorty (void);
  void border (void);
  void tsave (void);
  void help (void);
  void vers (void);
  
  for(i=1;i<argc;++i) {
    if(! strncmp(argv[i],"-",1)) {
      if(! strcmp(argv[i],"-m")) {if(argc>i+1){++i; mdl = atoi(argv[i]);}}
      if(! strcmp(argv[i],"-s")) {if(argc>i+1){++i; reg = atoi(argv[i]);
      reg &= 15; selec &= 240; selec |= reg;}}
      if(! strcmp(argv[i],"-t")) {if(argc>i+1){++i; reg = atoi(argv[i]);
      reg &= 15; sout &= 240; sout |= reg;}}
      if(! strcmp(argv[i],"-o")) {selec |= 32;}
      if(! strcmp(argv[i],"-r")) {if(argc>i+1){++i; res = atoi(argv[i]);
      if (res <= 0) {exit(EXIT_FAILURE);}}}
      if(! strcmp(argv[i],"-R")) {if(argc>i+1){++i; RES = atof(argv[i]);
      if (RES <= 0) {exit(EXIT_FAILURE);}}}
      if(! strcmp(argv[i],"-x")) {if(argc>i+4){++i; xmin=atof(argv[i]); ++i;
      ymax=atof(argv[i]); ++i; xmax=atof(argv[i]); ++i; ymin=atof(argv[i]);
      ymax *= -1; ymin *= -1; selec &= 191;}}
      if(! strcmp(argv[i],"-i")) {if(argc>i+1){++i;infl = atof(argv[i]);
      if (infl == 0) {selec &= 127;} else {selec |= 128;}}}
      if(! strcmp(argv[i],"-b")) {selec &= 239;}
      if(! strcmp(argv[i],"-h")) {help();exit(0);}
      if(! strcmp(argv[i],"--help")) {help();exit(0);}
      if(! strcmp(argv[i],"-V")) {vers();exit(0);}
      }
    else if(strlen(iname)) {if(!strlen(oname)) {strcpy(oname, argv[i]);}}
    else {strcpy(iname, argv[i]);}
  }
  if(!strlen(oname)) {help();exit(EXIT_FAILURE);}
  if(mdl == 1) {selec &= 228;}
  if(! selec & 15) {return(0);} /*nothing to do*/

  dload();
  sorty();
  mpic(mdl);
  if(selec & 16 && selec & 14){border();}/*add borders for g/c-save if wanted*/
  if(selec & 4) {gsave();}
  if(selec & 8) {csave();}
  if(selec & 2) {bsave();}
  /*bsave at last of the pictur-saves because it destroys the other data*/
  if(selec & 1) {tsave();}

  return(0);
}

void sorty (void) /*sorts ysort for the min y influence*/
{
  int i;
  int ycmp(const void *aa, const void *ba);
  qsort(ysort, n, sizeof(unsigned long), ycmp);
}

int ycmp(const void *aa, const void *ba)
{
  int *a = (int *)aa;
  int *b = (int *)ba;
  if(y[*a]-max[*a] > y[*b]-max[*b]) return (1);
  if(y[*a]-max[*a] < y[*b]-max[*b]) return (-1);
  else return(0);
}

void sortx (void) /*sorts ysort for the min y influence*/
{
  int i;
  int xcmp(const void *aa, const void *ba);
  qsort(stab, nstab, sizeof(unsigned long), xcmp);
}

int xcmp(const void *aa, const void *ba)
{
  int *a = (int *)aa;
  int *b = (int *)ba;
  if(x[*a]-max[*a] > x[*b]-max[*b]) return (1);
  if(x[*a]-max[*a] < x[*b]-max[*b]) return (-1);
  else return(0);
}

void mpic (int mdl) /*makes the picture
		      mdl..calcmode
		      1..wzp, 2..circle, 3..hyperbola, 4..poygon
		    */
{
  int px,py; /*x and y coordinate of the pixel*/
  int b; /*counter for the treenumber*/
  double f,g; /*f..picturresizefactor and later treedistance
		g..theoreticaly minimum activation in the pic*/
  void sortx (void);
  unsigned long i; /*index*/
  unsigned long batab; /*beginning number of atab*/
  unsigned long fatab; /*flip-register for stab*/
  float spic; /*storage for the actual picture activation*/
  unsigned long spicnr; /*storage for the actual picture number*/
  unsigned long *xbord, ybord; /*Border-finding registers*/

  if(RES>0) { /*set the resolution if the pixelsize is set*/
    res = (xmax-xmin)/RES;
    if((ymax-ymin)/RES>res) {res = (ymax-ymin)/RES;}
  }

  if(xmax-xmin <= 0 || ymax-ymin <= 0) {exit(EXIT_FAILURE);}
  f=res/(xmax-xmin);
  if(f > res/(ymax-ymin)){
    f = res/(ymax-ymin); yres = res;
    xres = ceil(res*((xmax-xmin)/(ymax-ymin)));}
  else {xres = res; yres = ceil(res*((ymax-ymin)/(xmax-xmin)));}

  if((xbord = (unsigned long *) calloc(xres, sizeof(unsigned long))) == NULL)
    {exit(EXIT_FAILURE);} /*resize the xpictureborder storage register*/
  
  if(selec & 14) { /*only for picturoutput*/	  
    if((pic = (float *) malloc(xres*yres* sizeof(float))) == NULL)
      {exit(EXIT_FAILURE);} /*resize the picturedata storage register*/
    if((picnr = (unsigned long *) malloc(xres*yres* sizeof(unsigned long))) == 
       NULL) {exit(EXIT_FAILURE);} /*resize the picturedata storage register*/
  }

  switch (mdl) { /*search the minimum activation*/
  case 1 : g = 0; break;
  case 2 : g = (max[0]/w[0])*(-1); break;
  case 3 : g = w[atab[0]]-max[0]; break;
  case 4 : g = -max[0]; break;
  default: g=0;
  }
  for(b=0;b<n;++b) { /*transform the coordinates and search min activation*/
    x[b] = (x[b]-xmin)*f;
    y[b] = (y[b]-ymin)*f;
    max[b] *= f;
    switch (mdl) {
    case 2 : if(g > (max[b]/w[b])*(-1)) {g=(max[b]/w[b])*(-1);}
      break;
    case 3 : if(g > w[atab[b]]-max[b]) {g = w[atab[b]]-max[b];}
      break;
    case 4 : if(g > (-max[b])) {g = -max[b];}
      break;
    default: {}
    }
  }

  nysort = 0;

  for(px=0;px<xres;++px) {  /*set the last tree number to no tree*/
    xbord[px] = n+1;
  }

  for ( py = 0; py < yres; ++py )
    {
      i = 0;
      for(b=0;b<nstab;++b) { /*drop trees without influence in py*/
	if(y[stab[b]]+max[stab[b]] < py) {continue;}
	stab[i] = stab[b];++i;
      }
      nstab = i;
      for(;nysort<n;++nysort) { /*select trees with influence in py*/
	if(y[ysort[nysort]]-max[ysort[nysort]] > py) {break;}
	stab[nstab] = ysort[nysort];++nstab;
      }
      sortx();
      natab = 0;
      batab = 0;
      ybord = n+1; /*set the last tree number to no tree*/
      for(px=0;px<xres;++px)
	{
	  spic = g;/*init the pixel with the minimum activation*/
	  spicnr = n;/*init the pixel with no tree*/
	  for(;natab<nstab;++natab) { /*select trees with influence in px*/
	    if(x[stab[natab]]-max[stab[natab]] > px) {break;}
	    atab[natab] = stab[natab];
	    /*atab is only therefor, that stab is not so much unsorted.
	      maybe this don't increase the speed!!*/
	  }
	  for(b=batab;b<natab;++b) {
	    /*some models can work without the sqrt but now not released*/
	    f = sqrt(((x[atab[b]]-px)*(x[atab[b]]-px))+
		     ((y[atab[b]]-py)*(y[atab[b]]-py)));
	    if(max[atab[b]]<f) {
	      if(x[atab[b]] < px) {
		if(b == batab) {++batab;}
		else {fatab = atab[b]; atab[b] = atab[batab];
		atab[batab] = fatab; ++batab;}
	      }
	      continue;
	    }
	    /*for speed you can reduche only to one modell*/
	    switch (mdl) {
	    case 1 : ++spic; break;
	    case 2 :
	      if((f/w[atab[b]])*(-1) > spic) {spic = (f/w[atab[b]])*(-1);
	      spicnr = atab[b];}
	      break;
	    case 3 :
	      if(w[atab[b]]-f > spic) {spic = w[atab[b]]-f;
	      spicnr = atab[b];}
	      break;
	    case 4 : if(-f > spic) {spic = -f; spicnr = atab[b];}
	      break;
	    default: {}
	    }
	  }

	  if(xbord[px] == n || ybord == n) {bord[spicnr] = 1;} /*bordertree*/
	  if(spicnr == n) {
	    if(xbord[px] <= n) {bord[xbord[px]] = 1;}
	    if(ybord <= n) {bord[ybord] = 1;}
	  }
	  xbord[px] = spicnr; ybord = spicnr;

	  ++area[spicnr]; /*increase the area of tree x*/
	  if(spic > maa) {maa = spic;} /*max activation*/
	  if(spic < mia) {mia = spic;} /*min activation*/
	  if(selec & 14) { /*only for picturoutput*/	  
	    pic[py*xres+px] = spic;
	    picnr[py*xres+px] = spicnr;
	  }
	}
    }
}

void border (void) /*add a border of the standings*/
{
  unsigned long i;
  unsigned char m, full;

  for(i=0;i<xres;++i) { /*top-bottom-border controlling*/
    if(picnr[i] == n) {pic[i] = mia;}
    if(picnr[(xres*yres)-xres+i] == n) {pic[(xres*yres)-xres+i] = mia;}
  }
  for(i=0;i<yres;++i) { /*left-right-border controlling*/
    if(picnr[i*xres] == n) {pic[i*xres] = mia;}
    if(picnr[(i*xres)+xres-1] == n) {pic[(i*xres)+xres-1] = mia;}
  }
  for(i=0;i<xres*yres-xres;++i) {
    if(picnr[i] != picnr[i+1] || picnr[i] != picnr[i+xres]) {pic[i] = mia;}
  }
}

void tsave (void) /*text-save*/
{
  unsigned long i;
  double m,a;

  if(sout & 1) { /*calculate the real area not the pixels*/
    if(xres < 2) {exit(EXIT_FAILURE);}
    if(yres < 2) {exit(EXIT_FAILURE);}
    m = ((xmax-xmin)/(xres-1));m *= ((ymax-ymin)/(yres-1));
    for(i=0;i<=n;++i) {area[i] *= m;}
  }
  
  if(selec & 32) { /*text to stdout and not to a file*/
    if(sout & 4) {
      printf("#min: %lf max: %lf\n", mia, maa);
      printf("#min:x/y %lf / %lf max:x/y %lf / %lf\n#res:x/y %d / %d\n",
	      xmin, ymax*(-1), xmax, ymin*(-1), xres, yres);
      printf("#outputselection: %d\n", sout);
      printf("#treenr\tnr of pixels\t[border]");
      printf("\nfree\t%.2f\t[1]", area[n]);
    }
    else {
      if(sout & 10) {printf("nr\tarea");}
      else {printf("nr\tarea\tborder");}
    }
    for(i=0;i<n;++i) {
      if(sout & 2) { /*drop the bordertrees*/
	if(bord[i] == 0) {
	  printf("\n%d\t%.2f", nr[i], area[i]);
	}
      }
      else {
	if(sout & 8) { /*replace area with NA*/
	  if(bord[i]) {
	    printf("\n%d\tNA", nr[i]);
	  } else {printf("\n%d\t%.2f", nr[i], area[i]);}
	} else {printf("\n%d\t%.2f\t%d", nr[i], area[i], bord[i]);}
      }
    }
  }
  else {
    strcpy(aname, oname);
    strcat(aname, ".txt");
    fp = fopen(aname, "wt");
    if(sout & 4) {
      fprintf(fp,"#min: %lf max: %lf\n", mia, maa);
      fprintf(fp,"#min:x/y %lf / %lf max:x/y %lf / %lf\n#res:x/y %d / %d\n",
	      xmin, ymax*(-1), xmax, ymin*(-1), xres, yres);
      fprintf(fp,"#outputselection: %d\n", sout);
      fprintf(fp,"#treenr\tnr of pixels\t[border]");
      fprintf(fp,"\nfree\t%.2f\t[1]", area[n]);
    }
    else {
      if(sout & 10) {fprintf(fp,"nr\tarea");}
      else {fprintf(fp,"nr\tarea\tborder");}
    }
    for(i=0;i<n;++i) {
      if(sout & 2) { /*drop the bordertrees*/
	if(bord[i] == 0) {
	  fprintf(fp,"\n%d\t%.2f", nr[i], area[i]);
	}
      }
      else {
	if(sout & 8) { /*replace area with NA*/
	  if(bord[i]) {
	    fprintf(fp,"\n%d\tNA", nr[i]);
	  } else {fprintf(fp,"\n%d\t%.2f", nr[i], area[i]);}
	} else {fprintf(fp,"\n%d\t%.2f\t%d", nr[i], area[i], bord[i]);}
      }
    }
    fclose(fp);
  }
}

void bsave (void) /*black-white-picture save*/
{
  unsigned long i,j;
  unsigned char m, full;

  strcpy(aname, oname);
  strcat(aname, ".pbm");
  fp = fopen(aname, "wb");
  fprintf(fp,"P4\n#min: %lf max: %lf\n", mia, maa);
  fprintf(fp,"#min:x/y %lf / %lf max:x/y %lf / %lf\n%d %d\n",
	  xmin, ymax*(-1), xmax, ymin*(-1), xres, yres);
  for(i=0;i<xres;++i) { /*top-bottom-border controlling*/
    if(picnr[i] < n) {pic[i] = 1;}
    else {pic[i] = 0;}
    if(picnr[(xres*yres)-xres+i] < n) {pic[(xres*yres)-xres+i] = 1;}
    else {pic[(xres*yres)-xres+i] = 0;}
  }
  for(i=0;i<yres;++i) { /*left-right-border controlling*/
    if(picnr[i*xres] < n) {pic[i*xres] = 1;}
    else {pic[i*xres] = 0;}
    if(picnr[(i*xres)+xres-1] < n) {pic[(i*xres)+xres-1] = 1;}
    else {pic[(i*xres)+xres-1] = 0;}
  }
  for(i=0;i<xres*yres-xres;++i) {
    if(picnr[i] == picnr[i+1] && picnr[i] == picnr[i+xres]) {pic[i] = 0;}
    else {pic[i] = 1;}
  }
  for(j=0;j<yres;++j) {
    full = 0;
    for(i=0;i<xres;++i) {
      ++full;
      m <<= 1;
      m |= (char)pic[(xres*j)+i];
      if(full == 8) {putc(m,fp); full = 0;}
    }
    if(full) {m <<= 8-full; putc(m,fp);}
  }
  fclose(fp);
}

void gsave (void) /*greyscale-picture save*/
{
  unsigned long i;
  double m;

  strcpy(aname, oname);
  strcat(aname, ".pgm");
  fp = fopen(aname, "wb");
  fprintf(fp,"P5\n#min: %lf max: %lf\n", mia, maa);
  fprintf(fp,"#min:x/y %lf / %lf max:x/y %lf / %lf\n%d %d 255\n",
	  xmin, ymax*(-1), xmax, ymin*(-1), xres, yres);
  m = 255/(maa-mia);
  for(i=0;i<xres*yres;++i) {
    putc((pic[i]-mia)*m,fp);
  }
  fclose(fp);
}

void csave (void) /*color-picture save*/
{
  unsigned long i;
  double m;

  strcpy(aname, oname);
  strcat(aname, ".ppm");
  fp = fopen(aname, "wb");
  fprintf(fp,"P6\n#min: %lf max: %lf\n", mia, maa);
  fprintf(fp,"#min:x/y %lf / %lf max:x/y %lf / %lf\n%d %d 255\n",
	  xmin, ymax*(-1), xmax, ymin*(-1), xres, yres);
  m = 255/(maa-mia);
  for(i=0;i<xres*yres;++i) {
    if(ba[picnr[i]] < n) {
      switch(ba[picnr[i]]) {
      case 1 : putc((pic[i]-mia)*m,fp);putc(0,fp);putc(0,fp); break;
      case 2 : putc(0,fp);putc((pic[i]-mia)*m,fp);putc(0,fp); break;
      case 3 : putc(0,fp);putc(0,fp);putc((pic[i]-mia)*m,fp); break;
      case 4 : putc((pic[i]-mia)*m,fp);putc((pic[i]-mia)*m,fp);putc(0,fp);
	break;
      case 5 : putc((pic[i]-mia)*m,fp);putc(0,fp);putc((pic[i]-mia)*m,fp);
	break;
      case 6 : putc(0,fp);putc((pic[i]-mia)*m,fp);putc((pic[i]-mia)*m,fp);
	break;
      default: putc((pic[i]-mia)*m,fp);putc((pic[i]-mia)*m,fp);
	putc((pic[i]-mia)*m,fp);
      }
    } else {putc(0,fp);putc(0,fp);putc(0,fp);}
  }
  fclose(fp);
}

void dload (void) /*loads the stand-data*/
     /*File format for inputdata:
       nr ba w x y max
       nr..treenumber (u.long = 0..4294967295)
       ba..species (integer = -32768 to 32767 bzw. -2147483648 to 2147483647)
       w..growing-factor (double = 2.3E-308 to 1.7E+308, 15 digits)
       x..x-cordinate (double)
       y..y-cordinate (double)
       max ..maximum influenceradius (double)*/
{
  unsigned long i;
  nn = 1000; /*first size of the data-size (autoexpanding on loading)*/

  if((nr = (unsigned long *) malloc(nn * sizeof(unsigned long))) == NULL)
    {exit(EXIT_FAILURE);}
  if((ba = (int *) malloc(nn * sizeof(int))) == NULL) {exit(EXIT_FAILURE);}
  if((w = (double *) malloc(nn * sizeof(double))) == NULL)
    {exit(EXIT_FAILURE);}
  if((x = (double *) malloc(nn * sizeof(double))) == NULL)
    {exit(EXIT_FAILURE);}
  if((y = (double *) malloc(nn * sizeof(double))) == NULL)
    {exit(EXIT_FAILURE);}
  if((max = (double *) malloc(nn * sizeof(double))) == NULL)
    {exit(EXIT_FAILURE);}

  if((fp=fopen(iname, "rt"))==NULL)
    {(void)puts("File not found"); exit(EXIT_FAILURE);}

  n = 0;
  while((fscanf(fp,"%lu %d %lf %lf %lf %lf\n",
		&nr[n],&ba[n],&w[n],&x[n],&y[n],&max[n])) != EOF) {
    y[n] *= -1;/*flip the y axe that x-Axe is left and y-Axe is up in the pic*/
    if(selec & 128) { /*calculate the maximum influence radius*/
      max[n] = w[n]/infl;}
    if(selec & 64) { /*calculate max/min x/y*/
      if(xmin > x[n]-max[n]) {xmin = x[n]-max[n];}
      if(xmax < x[n]+max[n]) {xmax = x[n]+max[n];}
      if(ymin > y[n]-max[n]) {ymin = y[n]-max[n];}
      if(ymax < y[n]+max[n]) {ymax = y[n]+max[n];}
    }
    ++n;
    if(n==nn) {
      nn += 1000;
      if((nr = (unsigned long *) realloc(nr, nn * sizeof(unsigned long))) == 
	 NULL) {exit(EXIT_FAILURE);}
      if((ba = (int *) realloc(ba, nn * sizeof(int))) == NULL)
	{exit(EXIT_FAILURE);}
      if((w = (double *) realloc(w, nn * sizeof(double))) == NULL)
	{exit(EXIT_FAILURE);}
      if((x = (double *) realloc(x, nn * sizeof(double))) == NULL)
	{exit(EXIT_FAILURE);}
      if((y = (double *) realloc(y, nn * sizeof(double))) == NULL)
	{exit(EXIT_FAILURE);}
      if((max = (double *) realloc(max, nn * sizeof(double))) == NULL)
	{exit(EXIT_FAILURE);}
    }
  }

  if((area = (float *) calloc(n+1, sizeof(float))) == NULL)
    {exit(EXIT_FAILURE);}
  if((ysort = (unsigned long *) malloc(n * sizeof(unsigned long))) == NULL)
    {exit(EXIT_FAILURE);}
  if((stab = (unsigned long *) malloc(n * sizeof(unsigned long))) == NULL)
    {exit(EXIT_FAILURE);}
  if((atab = (unsigned long *) malloc(n * sizeof(unsigned long))) == NULL)
    {exit(EXIT_FAILURE);}
  if((bord = (char *) calloc(n, sizeof(char))) == NULL)
    {exit(EXIT_FAILURE);}
  for(i=0;i<n;++i) {ysort[i] = i;}

  fclose(fp);
}

void vers()
{
  puts("acre 0.1a, Copyright (C) 2006 Georg Kindermann");
  puts("acre comes with ABSOLUTELY NO WARRANTY; for details type `more copyright'.");
  puts("This is free software, and you are welcome to redistribute it");
  puts("under certain conditions; type `more copyright' for details.");

  puts("If you find acre useful, or you have some proplems with it, or you have some");
  puts("idea how it could work better, or ...   write to");
  puts("Post: Georg Kindermann, Institut fuer Waldwachstumsforschung - BOKU");
  puts("      Peter Jordanstrasse 82, A-1190 Wien, Austria");
  puts("Home: https://sourceforge.net/projects/acre");
  puts("Mail: georg.kindermann@boku.ac.at\n");
  puts("To cite acre in publications use: Georg Kindermann (1998),");
  puts(" Die Fl\"achenanteile der Baumarten (The areashare of tree-species),");
  puts(" Masterthesis, Istitut of forest growth and yield,");
  puts(" University of applied life science (Universit\"at f\"ur Bodenkultur),");
  puts(" Vienna, Austria.");
  puts("I have invested some time and effort in creating acre, please cite");
  puts("it when using it.");
}

void help()
{
puts("\nSyntax: acre INFILE OUTFILE(S) [OPTIONS]");
puts("Format of the INFILE:");
puts("            *)tree-number  *)species     *)growingsize");
puts("            *)X-position   *)Y-position  *)distance of influence");
puts("Format of the text-OUTFILE:");
puts("            [first lines with # ..info of the calculation]");
puts("            *)tree-number  *)pixels/area [*)bordertree]");
puts("options:");
puts(" -m <nr>  ..modell: 1..layers(wzp,ccf..), 2..circle (DEF), 3..hyperbola");
puts("                    4..polygon");
puts(" -s <nr>  ..selection: 1..textoutput, 2..s/wpic, 4..greypic, 8..colpic");
puts("             you can add the numbers for multi-output (DEF=15)");
puts(" -o       ..output of text to stdout (DEF not set)");
puts(" -r <nr>  ..resolution: set the maximum resolution of pic (DEF=400)");
puts(" -R <nr>  ..pixelsize: set the pixelsize (overrides the settings of -r)");
puts(" -x <minx> <miny> <maxx> <maxy> ..set the border to x/y");
puts("                                  (DEF..not set -> auto fit)");
puts(" -i <nr>  ..influence: set the max influenceradius to growingfactor/<nr>");
puts("            if <nr> = 0 ... use the radius from the inputdatafile (DEF=0)");
puts(" -b <nr>  ..border(not mod 1): set the pic-borders off (DEF borders on)");
puts(" -t <nr>  ..asciioutput: 1..pixels(0)/area(1), 2..drop bordertrees,");
puts("                         4..show header, 8..area of bordertrees as NA");
puts("             you can add the numbers (DEF=5)");
puts(" -h <nr>  ..help: display this and exit");
puts(" -V <nr>  ..display Version, copyright info and exit");
puts("example: acre data.txt result");
}
