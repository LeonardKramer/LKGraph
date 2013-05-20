#include "LKGraph.h"
#include <math.h>

double _Xc(LKGraph &G, double x);
double _Xc_withclipping(LKGraph &G, double x);
double _YC(LKGraph &G, double x);
double _YC_withclipping(LKGraph &G, double x);


double (*p_Xc)(LKGraph &G, double x);
double (*p_Yc)(LKGraph &G, double x);
#define XC(X) (p_Xc(*this,X))
#define YC(X) (p_Yc(*this,X))

double _Xc_withclipping(LKGraph &G, double x) { 
	double R = G.cxmin + (x-G.xmin)*(G.cxmax-G.cxmin)/(G.xmax-G.xmin);
	if (G.xclip) {
		R = (R > G.cxmin) ?  R : G.cxmin;
		R = (R < G.cxmax) ?  R : G.cxmax;
	}
	return R;
}
double _Xc(LKGraph &G, double x) { 
	return G.cxmin + (x-G.xmin)*(G.cxmax-G.cxmin)/(G.xmax-G.xmin);

}
double _YC(LKGraph &G, double y) { 
	return G.cymin + (y-G.ymin)*(G.cymax-G.cymin)/(G.ymax-G.ymin);
}
double _YC_withclipping(LKGraph &G, double y) { 
	double R=G.cymin + (y-G.ymin)*(G.cymax-G.cymin)/(G.ymax-G.ymin);
	if (G.yclip) {
		R = (R > G.cymin) ? G.cymin : R;
		R = (R < G.cymax) ? G.cymax : R;
	}
	return R;
}

void LKGraph::set_xclip(int n) {
	xclip = n;
	p_Xc = xclip ? _Xc_withclipping: _Xc ;
}
void LKGraph::set_yclip(int n) {
	yclip = n;
	p_Yc = yclip ? _YC_withclipping  :  _YC ;
}
void LKGraph::map(){
	if (!noerase) erase();
	set_xy(xmin,xmax,ymin,ymax);
	frame();
	cairo_set_dash(cr,dashes[linetype],ndash[linetype],0.00);
	cairo_set_source(cr,foreground);
	cairo_set_line_width(cr,linewidth);

	for (int i=0;i<this->numMapNodes;i++) {
		if (this->mapArray[i].pen && (( this->mapArray[i].y-.5) >xmin )&& ((this->mapArray[i].y+.5) < xmax) ) 
			cairo_line_to(cr,XC(this->mapArray[i].y),YC(this->mapArray[i].x));
		else {
			cairo_stroke(cr);
			cairo_move_to(cr,XC(this->mapArray[i].y),YC(this->mapArray[i].x));
		}
	}

	cairo_stroke(cr);

}

void (*draw_o_move[2])(cairo_t *cr,double x, double y);

struct SymbolData {
	int pen;
	double x; double y;
};
const double d=8.4853;
static SymbolData Box[]     = { 0,3,3, 1,-6,0,  1,0,-6,   1,6,0,   1,0,6,  -1,0,0};
static SymbolData Diamond[] = { 0,3,0, 1,-3,3,  1,-3,-3,  1,3,-3,  1,3,3,  -1,0,0};
static SymbolData Ex[]		= { 0,3,3, 1,-6,-6, 0,0,6,    1,6,-6,  -1,0,0 };
static SymbolData Cross[]   = { 0,-3,0,  1,6,0,  0,-3,-3,  1,0,6,  -1,0,0};
static SymbolData Dash[]    = { 0,-3,0,  1,6,0, -1, 0, 0};
static SymbolData Aster[]   = { 0,-3,0,  1,6,0,  0,-3,-3,  1,0,6, 0,3,0, 1,-6,-6, 0,0,6,    1,6,-6,  -1,0,0 };
static SymbolData Dot[] = {0,-1,0, 1, 1,1, 1, 1,-1,  1,-1,-1, 1,-1,1, -1,0,0};
int MaxNumPsym;
SymbolData *Symbol[64];
class StartUp {
public:
	StartUp() {
		draw_o_move[0]=cairo_rel_move_to;
		draw_o_move[1]=cairo_rel_line_to;
		Symbol[0] = Cross;   // not used - psym = 0 signals none.
		Symbol[1] = Cross;
		Symbol[2] = Ex;
		Symbol[3] = Box;
		Symbol[4] = Diamond;
		Symbol[5] = Dash;
		Symbol[6] = Aster;
		Symbol[7] = Dot;
		MaxNumPsym = 8;
	}
} Start;

void DrawSymbol(cairo_t *cr, SymbolData *S) {
	for (int i=0;(S[i].pen>=0);i++) draw_o_move[S[i].pen](cr,S[i].x,S[i].y);
}
LKGraph::LKGraph(double _width, double _height) {
//	if (screen) cleanup();
	init(_width,_height);
	Initialized = 1;
}
LKGraph::LKGraph(void) {
	init(640,480);
	Initialized = 1;
}
void LKGraph::set_outputfilename(char *_pngfilename){
		strcpy(pngfilename,_pngfilename);
}
void LKGraph::erase() {
		unsigned int old = set_color(0xFFFFFFFF);
		cairo_rectangle(cr,0,0,width,height);
		cairo_set_source(cr,foreground);
		cairo_fill(cr);
		set_color(old);
}
void LKGraph::plot(double *x, double *y, int N) {
	if(!overplot) {
		if (xmin == 0.0 && xmax == 0.0) minmaxdata(x,N,&xmin,&xmax);
		if (ymin == 0.0 && ymax == 0.0) minmaxdata(y,N,&ymin,&ymax);
		if (ymin == ymax) {ymin=ymin-1.; ymax=ymax+1.0;}
		if (!noerase) erase();
		set_xy(xmin,xmax,ymin,ymax);
		frame();
		xaxis();
		yaxis();
		maintitle();
		set_color(color);
	}
	if (!nodata) oplot(x,y,N);
}

void LKGraph::oplot(double *x, double *y, int N) {
	cairo_set_dash(cr,dashes[linetype],ndash[linetype],0.00);
	cairo_set_line_width(cr,linewidth);
	cairo_set_source(cr,foreground);
	if (N==0) return;
	if (psym != 0) {
		for (int i=0;i<N;i++ ) {
			cairo_move_to(cr,XC(x[i]),YC(y[i])) ;
			DrawSymbol(cr, Symbol[psym]);
			if ((i % 4096) == 0) {
				if (i) printf("%d of %d points\n",i,N);
				fflush(stdout);		
				cairo_stroke(cr);
			}
		}
		cairo_stroke(cr);
	}
	if (linetype != 0) {
		cairo_line_to(cr,XC(x[0]),YC(y[0])) ;
		for (int i=1;i<N;i++ ) {
			cairo_line_to(cr,XC(x[i]),YC(y[i])) ;
			if ((i % 4096) == 0) {
				printf("%d of %d points\n",i,N);
				fflush(stdout);		
				cairo_stroke(cr);
			}
		}
		cairo_stroke(cr);
	}
}
unsigned int LKGraph::set_color(unsigned int hexcolor) {
	unsigned char *bgra = (unsigned char*)&hexcolor;
	double a = (double)bgra[3]/255.0;
	double r = (double)bgra[2]/255.0;
	double g = (double)bgra[1]/255.0;
	double b = (double)bgra[0]/255.0;
	foreground = cairo_pattern_create_rgba(r,g,b,a);
	unsigned int oldcolor = color;
	color = hexcolor;
	return oldcolor;
}
void LKGraph::set_charsize(double _charsize) {		
	charsize = _charsize;
	cairo_set_font_size(cr,charsize);
}



void LKGraph::plot_oi(double *x, double *y, int N) {
	double *lnx = new double[N];
	for (int i=0;i<N;i++) lnx[i] = log10(x[i]);
	plot(lnx,y,N);
	delete lnx;
}
void LKGraph::plot_oo(double *x, double *y, int N) {
	double *lnx = new double[N];
	double *lny = new double[N];
	for (int i=0;i<N;i++) {
		lnx[i] = (x[i] > 0) ? log10(x[i]) : -64.0;
		lny[i] = (y[i] > 0.) ? log10(y[i]) : -64.0;
	}

	logplot=1;
	plot(lnx,lny,N);
	logplot=0;
	delete lny, lnx;
}

void LKGraph::maintitle(void) {
	if (mtitle[0] != '\0') {
		double lwdth, lhght;
		double xwidth, xheight;
		c_text_extents("X",&xwidth,&xheight);
		c_text_extents(mtitle,&lwdth, &lhght);
		double cx = (XC(xmin)+XC(xmax)-lwdth)/2.0;
		double cy = 2*xheight;
		c_moveto(cx,YC(ymax)-cy);
		c_text(mtitle,0.0);
		cairo_stroke(cr);
	}
}
void LKGraph::xaxis(void ) {
	double cx;
	double lwdth, lhght;
	double xwidth, xheight;
	c_text_extents("X",&xwidth,&xheight);
	double dy= YC(0.05*(ymax-ymin));
	for (int i=1; i<Nxticks; i++) {
		cx = XC(xmin+i*(xmax-xmin)/Nxticks);
		double y0 = YC(ymin);
		double y1 = YC(ymin+ 0.03*(ymax-ymin));
		c_moveto(cx,y0);
		c_lineto(cx,y1);
	}
	for (int i=0;i<=Nxticks; i++) {
		double x =xmin+i*(xmax-xmin)/Nxticks;
		char xlabel[256];
		if (logplot) {
			sprintf(xlabel,"%G",pow(10,x));
		} else {
			sprintf(xlabel,"%G",x);
		}
		c_text_extents(xlabel,&lwdth, &lhght);
		cx = XC(x)-lwdth/2;
		double y1 = YC(ymin)+2.0*xheight;
		c_moveto(cx,y1);
		c_text(xlabel,0.0);
	}
	if (xtitle[0] != '\0') {
		c_text_extents(xtitle,&lwdth, &lhght);
		cx = (XC(xmin)+XC(xmax)-lwdth)/2.0;
		double cy = 4*xheight;
		c_moveto(cx,YC(ymin)+cy);
		c_text(xtitle,0.0);
	}
	cairo_stroke(cr);
}
void LKGraph::yaxis(void ) {
	double cy;
	double lwdth, lhght;
	double xwidth, xheight;
	c_text_extents("X",&xwidth,&xheight);

	for (int i=1; i<Nyticks; i++) {
		cy = YC(ymin+i*(ymax-ymin)/Nyticks);
		double cx = XC(xmin);
		double cx1 = XC(xmin) + 0.03*(cymin - cymax);
		c_moveto(cx,cy);
		c_lineto(cx1,cy);
	}
	double mincx=XC(xmin);
	for (int i=0;i<=Nyticks; i++) {
		double y =ymin+i*(ymax-ymin)/Nyticks;
		char ylabel[256];
		if (logplot) {
			sprintf(ylabel,"%G",pow(10,y));
		} else {
			sprintf(ylabel,"%G",y);
		}
		c_text_extents(ylabel,&lwdth, &lhght);
		double cx = XC(xmin)-lwdth-xwidth;
		if (cx < mincx) mincx = cx;
		double cy = YC(y) + lhght/2;
		c_moveto(cx,cy);
		c_text(ylabel,0.0);
	}
	if (ytitle[0] != '\0') {
		c_text_extents(ytitle,&lwdth, &lhght);
		cy = (YC(ymin)+YC(ymax)+lwdth)/2.0;
		double cx = mincx-2.*xheight;
		xyoutstring(cx,cy,ytitle,90.0);
		//c_moveto(cx,cy);
		//c_text(ytitle,90.0);
	}
	cairo_stroke(cr);
}
void LKGraph::xyoutstring(double _cx,double _cy,char *string, double angle) {
	double lwdth, lhght;
	double xwidth, xheight;
	c_text_extents(string,&lwdth, &lhght);
	c_moveto(_cx,_cy);
	c_text(string,angle);
}
void LKGraph::set_xticks(int _Nxticks) { Nxticks = _Nxticks;}
void LKGraph::set_yticks(int _Nyticks) { Nyticks = _Nyticks;}
void LKGraph::set_zticks(int _Nzticks) { Nzticks = _Nzticks;} 
void LKGraph::frame(void ) {
	double oldlinewidth = linewidth;
	linewidth = 1.5*linewidth;
	c_moveto(XC(xmin),YC(ymin));
	c_lineto(XC(xmin),YC(ymax));
	c_lineto(XC(xmax),YC(ymax));
	c_lineto(XC(xmax),YC(ymin));
	c_lineto(XC(xmin),YC(ymin));
	cairo_stroke(cr);
	linewidth = oldlinewidth;
}
void LKGraph::c_moveto(double x,double y) {  // hardware coordinate move.
		cx0=x;
		cy0=y;
	}
void LKGraph::c_lineto(double x, double y) {
		cairo_set_source(cr,foreground);
		cairo_set_line_width(cr,linewidth);
		cairo_move_to(cr,cx0,cy0);
		cairo_line_to(cr,x,y);
		cx0=x; cy0=y;
	}
void LKGraph::c_text_extents(char *string, double *width, double *height) {
		cairo_text_extents_t extents;
		cairo_text_extents(cr,string,&extents);
		*width = extents.width;
		*height = extents.height;
	}
void LKGraph::c_text(char *string,double angle){
		cairo_set_source(cr,foreground);
		cairo_move_to(cr,cx0,cy0);
		cairo_rotate(cr,-angle*pi/180.0);
		cairo_show_text(cr,string);
		cairo_stroke(cr);
		cairo_rotate(cr,angle*pi/180.0);
	}
void LKGraph::set_viewport(double fxmin, double fxmax, double fymin, double fymax) {
		// default viewport { 0.2, .95, .15, .90 }  %height - .75  & %width = 0.75
		cxmin = 0 + fxmin*width;
		cxmax = 0 + fxmax*width;
		cymin = height - fymin*height;
		cymax = height - fymax*height;
		double new_width_frac = (cxmax-cxmin)/width;
		double new_height_frac = (cymin-cymax)/height;
		double frac = ((new_width_frac < new_height_frac)? new_width_frac : new_height_frac)/0.75;
		
		set_charsize(14*frac);
	}
void LKGraph::set_xy(double _xmin,double _xmax, double _ymin, double _ymax) {
	xmin = _xmin;
	xmax = _xmax;
	ymin = _ymin;
	ymax = _ymax;
}

void LKGraph::init(double _width, double _height){
	logplot=0;   // selects linear plot (not logplot) by default.
	p_Yc=_YC;
	p_Xc=_Xc;
	overplot=0;
	xclip = 0;
	yclip = 0;
	Nzticks = 2;
	nodata = 0;
	noerase = 0;
	width = _width;
	height = _height;
	white = cairo_pattern_create_rgba(1.,1.,1.,1.);
	black = cairo_pattern_create_rgba(1.0,0.,0.,0.0);
	red = cairo_pattern_create_rgba(1.0,0.,0.,1.0);
	green = cairo_pattern_create_rgba(0.,1.0,0.,1.0);
	blue = cairo_pattern_create_rgba(0.,0.,1.0,1.0);
	color = 0xFF000000;
	background = white;
	foreground = black;
	cx0=0;
	cy0=0;
	image = NULL;
	screen = NULL;
	surface = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, (int)width, (int)height);
	if (surface < 0) {
		fprintf(stderr, "LKGraph: cairo_image_surface_create: format invalid.\n");
		exit(1);
	}
	cr=cairo_create(surface);
	if ( cairo_status(cr) == CAIRO_STATUS_NO_MEMORY) {
		fprintf(stderr, "LKGraph: cairo_create: no memory.\n");
		exit(1);
	}
	cairo_rectangle(cr,0,0,width,height);
	cairo_set_source(cr,background);
	cairo_fill(cr);

	cairo_select_font_face(cr,"serif",CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
	set_charsize(14);
	set_viewport(0.2,.95,.15, .90);
	set_xy(0.,0.0,0.0, 0.0);

	Nxticks = 5;
	Nyticks = 5;
	linewidth = 1;

	xtitle[0]='\0';
	ytitle[0]='\0';
	mtitle[0]='\0';
	linetype = 0;

	dashes[0] = dashes0;
	dashes[1] = dashes0;
	dashes[2] = dashes1;
	dashes[3] = dashes2;
	ndash[0]=ndash0;
	ndash[1]=ndash0;
	ndash[2]=ndash1;
	ndash[3]=ndash2;
	psym = 1;
	strcpy(pngfilename,"LKGraph.png");

	FILE *fp = fopen(mapdatafile,"rb");
	if (fp == NULL) {
		fprintf(stderr,"LKGraph.cpp: init() could not open %s\n",mapdatafile);
		exit (1);
	}
	fread(&numMapNodes,sizeof(numMapNodes),1,fp);
	mapArray = new MapData[numMapNodes];
	fread(mapArray,sizeof(mapArray[0]),numMapNodes,fp);
	fclose(fp);
}
void LKGraph::set_psym(int _psym) { psym = _psym % MaxNumPsym;};
void LKGraph::set_linetype(int _linetype) {
	linetype=_linetype;
}
void LKGraph::set_linewidth(double _linewidth) {linewidth = _linewidth;}
void LKGraph::set_xtitle(char *_xtitle) {
	strcpy(xtitle,_xtitle);
}
void LKGraph::set_mtitle(char *_mtitle) {
	strcpy(mtitle,_mtitle);
}
void LKGraph::set_ytitle(char *_ytitle) {
	strcpy(ytitle,_ytitle);
}
LKGraph::~LKGraph(){
	delete mapArray;
	if (screen != NULL) SDL_Quit();
	cleanup();
}
void LKGraph::show(unsigned int delay) {
	cairo_stroke(cr);
	if (surface<0) {
		fprintf(stderr, "LKGraph.show: surface invalid\n");
		exit(1);
	}
	cairo_surface_write_to_png (surface, pngfilename);
	image = IMG_Load(pngfilename);
	if (image != NULL) {
		if (screen == NULL) {
			SDL_Init(SDL_INIT_EVERYTHING);
			screen = SDL_SetVideoMode(image->clip_rect.w,image->clip_rect.h,32,SDL_SWSURFACE);
		}
		SDL_BlitSurface(image,NULL,screen,NULL);
		SDL_Flip(screen);
		SDL_Event event;
		int Done = 0;
		while (!Done) {
			if (delay == 0xFFFFFFFF) {
				SDL_WaitEvent(&event);
			} else {
				SDL_PollEvent(&event);
				SDL_Delay(delay);
				Done = 1;
			}
			switch(event.type) {
				case SDL_QUIT:
					Done = 1;
					break;
				default:
					break;
			}
		}
		SDL_FreeSurface(image);
	}
}
int LKGraph::inside(double Q, double X0, double X1, double X2, double X3) {
	double _min = X0;
	if (X1 < _min) _min = X1;
	if (X2 < _min) _min = X2;
	if (X3 < _min) _min = X3;
	double _max = X0;
	if (X1 > _max) _max = X1;
	if (X2 > _max) _max = X2;
	if (X3 > _max) _max = X3;
	return ((_min <= Q) && (Q < _max));
}
double c[2][2];
int between( double Z0, double Z, double Z1) {
	return (((Z0 <= Z) && (Z < Z1)) || ((Z1 < Z) && (Z <= Z0))); 
}
int LocMax(double **C,int j,int i) {
	return ((C[j][i] >= C[j][i+1]) && (C[j][i] >= C[j+1][i])  && (C[j][i] >= C[j][i-1])  && (C[j][i] >= C[j-1][i]) ) ;
}
int LocMin(double **C,int j,int i) {
	return ((C[j][i] <= C[j][i+1]) && (C[j][i] <= C[j+1][i])  && (C[j][i] <= C[j][i-1])  && (C[j][i] <= C[j-1][i]) ) ;
}
void  LKGraph::DrawContour(double x[2], double y[2], double **C, int i, int j, double Z) {
	double Z0=C[i][j];
	double Z1=C[i][j+1];
	double Z2=C[i+1][j+1];
	double Z3=C[i+1][j];
	double dx = x[1]-x[0];
	double dy = y[1]-y[0];
	if ( between(Z0,Z,Z1) ) {
		double x0;
		x0= x[0]+ (Z-Z0)*(dx)/(Z1-Z0);
		c_moveto(x0,y[0]); 
		if (between(Z1,Z,Z2)) {
			double y1= y[0]+ (Z-Z1)*(dy)/(Z2-Z1);
			c_lineto(x[1],y1);
		} else if ( between(Z0,Z,Z3) ) {
			double y1= y[0]+ (Z-Z0)*(dy)/(Z3-Z0);
			c_lineto(x[0],y1);
		} else if ( between(Z3,Z,Z2)) {
			double x1 = x[0]+(Z-Z3)*(dx)/(Z2-Z3);
			c_lineto(x1,y[1]);
		}
	} else if (between(Z3,Z,Z2) ) {
		double x1= x[0]+(Z-Z3)*(dx)/(Z2-Z3);
		c_moveto(x1,y[1]);
		if (between(Z1,Z,Z2)) {
			double y1= y[0]+ (Z-Z1)*(dy)/(Z2-Z1);
			c_lineto(x[1],y1);
		}else if (between(Z0,Z,Z3)) {
			double y1= y[0]+ (Z-Z0)*(dy)/(Z3-Z0);
			c_lineto(x[0],y1);
		}
	} else if (between(Z0,Z,Z3) ) {
		double y1= y[0]+ (Z-Z0)*(dy)/(Z3-Z0);
		c_moveto(x[0],y1);
		y1 = y[0]+(Z-Z1)*(dy)/(Z2-Z1);
		c_lineto(x[1],y1);
	}
}
void LKGraph::contour(double **C, int N, int M) {
	double *x = new double[M];
	double *y = new double[N];
	for (int i=0;i<M;i++) x[i]=i;
	for (int i=0;i<N;i++) y[i]=i;
	this->set_xy(x[0],x[M-1],y[0],y[N-1]);
	this->xaxis();
	this->yaxis();
	maintitle();
	frame();
	set_color(color);
	cairo_set_dash(cr,dashes[linetype],ndash[linetype],0.00);
	cairo_set_line_width(cr,linewidth);
	zmin = zmax = C[0][0];
	for (int i=0;i<N;i++) {
		for (int j=0;j<M;j++) {
			if (C[i][j] < zmin) zmin = C[i][j];
			if (C[i][j] > zmax) zmax = C[i][j];
		}
	}
	double *zb = new double[Nzticks];
	double cx[2],cy[2];
	for (int i=0;i<Nzticks; i++) zb[i]=zmin+ (i+1)*(zmax-zmin)/(Nzticks+1);
	for (int i=0;i<Nzticks;i++) {
		double ZB = zb[i];
		for (int k=0;k<(N-1);k++) {
			for (int j=0;j<(M-1);j++) 
				if (inside(ZB,C[k][j],C[k][j+1],C[k+1][j+1],C[k+1][j])) {
					cx[0]=XC(x[j]); cx[1]=XC(x[j+1]);
					cy[0]=YC(y[k]); cy[1]=YC(y[k+1]);
				DrawContour(cx,cy,C,k,j,ZB);
			}
		}
	}
	delete zb;

	for (int j=1; j<N-1;j++) {
		for (int i=1; i< M-1; i++) {
			if (LocMax(C,j,i)) {
				cairo_move_to(cr,XC(x[i]),YC(y[j])) ;
				DrawSymbol(cr, Symbol[1]);  // Plus
			} else if (LocMin(C,j,i)) {
				cairo_move_to(cr,XC(x[i]),YC(y[j])) ;
				DrawSymbol(cr, Symbol[5]);  // Dash
			}
		}
	}

	delete x, y;
}
