//Author : Raphael Mahiet
//Hire Me : collaborate@designlense.io
//Website : designlense.io

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>
#include <stdlib.h>


using namespace std;

struct Couleur{
	double r;
	double g;
	double b;
};

struct ImgHeader {
	int nCol;
	int nLine;
};


typedef vector<vector<Couleur> > Image;

string uneSeuleLigne(string s){
	ifstream f(s.c_str());
	string line;
	string allLines;

	while(getline(f,line))
		allLines += line + " ";
	f.close();
	return allLines;
}

int readHeader(ImgHeader &header, string sSource)
{
    int nIndex ;

    header.nCol=0 ;
    header.nLine=0 ;
    nIndex=2;
    if ((sSource[0]=='P') && (sSource[1]=='3'))
    {

          int k ;

          k=0 ;
          for (nIndex++,header.nCol=0;(sSource[nIndex] >='0') && (sSource[nIndex] <= '9') && (k<70); nIndex++,k++)
             {header.nCol=header.nCol*10+(sSource[nIndex]-0x30) ;
             } ;

          for (nIndex++,header.nLine=0;(sSource[nIndex] >='0') && (sSource[nIndex] <= '9') && (k<70); nIndex++,k++)
             {header.nLine=header.nLine*10+(sSource[nIndex]-0x30) ;
             } ;
          if (k==70) {cout  <<"erreur ligne" << endl ;  } ;
		  nIndex++;
          if ((sSource[nIndex]='2') && (sSource[nIndex+1]=='5') && (sSource[nIndex+2]=='5'))
            return (nIndex+3) ;
          else
            return (0);
    }
    else return (0) ;
}

Image readPPM(string s){
    int nIndex;
    int nLineMax;
    int nColMax;
    Image img;
	int i;
	int j;
	int k;

    string  sSource ;
    struct ImgHeader header ;

	cout << s << endl;
    sSource = uneSeuleLigne(s) ;
    if (nIndex=readHeader(header,sSource))
    {

		nLineMax=header.nLine ;
		nColMax=header.nCol ;
		img.resize(nLineMax) ;
		for (i=0 ; i < nLineMax ; i++)
	       img[i].resize(nColMax) ;

        for (i=0 ; i < nLineMax ; i++)
        {
            for (j=0 ; j < nColMax ; j++)
                {

                    k=0 ;
					while ((sSource[nIndex]=='\t') || (sSource[nIndex]==' ')) nIndex++ ;
                    for (img[i][j].r=0;(sSource[nIndex] >='0') && (sSource[nIndex] <= '9') && (k<70); nIndex++,k++)
                    {  img[i][j].r=img[i][j].r*10+(sSource[nIndex]-0x30) ;
                    } ;
					while ((sSource[nIndex]=='\t') || (sSource[nIndex]==' ')) nIndex++ ;
                    for (img[i][j].g=0;(sSource[nIndex] >='0') && (sSource[nIndex] <= '9') && (k<70); nIndex++,k++)
                    {  img[i][j].g=img[i][j].g*10+(sSource[nIndex]-0x30) ;
                    } ;
					while ((sSource[nIndex]=='\t') || (sSource[nIndex]==' ')) nIndex++ ;
                    for (img[i][j].b=0;(sSource[nIndex] >='0') && (sSource[nIndex] <= '9') && (k<70); nIndex++,k++)
                    {  img[i][j].b=img[i][j].b*10+(sSource[nIndex]-0x30) ;
                    } ;
                    if (k==70) {cout  <<"erreur ligne" << endl; } ;

                }
        }
    }

	return img;
}

void writePPM(string s,Image img){



	ofstream fImg (s.c_str(),ofstream::out | ofstream::trunc);
	string line;
	int nCol ;
	int nLine ;
	int i,j ;



	nLine= img.size();
	nCol=img[0].size();

	fImg << "3P" << endl ;
	fImg << nCol << " " << nLine <<endl ;
	fImg <<"255"  << endl ;
	for (i=0; i< nLine ; i++)
	{
      for (j=0; j < nCol ; j++)
        {
            fImg << (int) img[i][j].r << " " << (int) img[i][j].g << " " << (int) img[i][j].b << endl;
        }
    }

	fImg.close();
	return ;
}

void testReadWriteImage(){
	cout << "verifier que vos images \"IO_couleur\" sont semblables a celles donnees en correction" << endl;
	writePPM("IO_couleur_Baboon.512.ppm",readPPM("Baboon.512.ppm"));
    writePPM("IO_couleur_Billes.256.ppm",readPPM("Billes.256.ppm"));
	writePPM("IO_couleur_Embryos.512.ppm",readPPM("Embryos.512.ppm"));
	writePPM("IO_couleur_House.256.ppm",readPPM("House.256.ppm"));
	writePPM("IO_couleur_Lena.512.ppm",readPPM("Lena.512.ppm"));
}

typedef vector<vector<double> > ImageGris;

ImageGris readGreyPPM(string s){
	int nLine;
	int nCol ;
	int i, j ;
	const double GREY_KRED = 2126 ;
	const double GREY_KGREEN = 7152 ;
	const double GREY_KBLUE= 722 ;
	Image imgRGB;
	ImageGris image;

	imgRGB = readPPM(s) ;
	nLine=imgRGB.size() ;
	nCol = imgRGB[0].size() ;

	image.resize(nLine) ;
	for (i=0 ; i <nLine; i++)
	  { image[i].resize(nCol) ;} ;

	for (i=0 ; i < nLine ; i++)
		for (j=0 ; j < nCol ; j++)
		{ image[i][j]= GREY_KRED*imgRGB[i][j].r + GREY_KGREEN*imgRGB[i][j].g + GREY_KBLUE*imgRGB[i][j].b ;
	      image[i][j] /= (GREY_KRED+GREY_KGREEN+GREY_KBLUE) ; } ;


	cout << " " << endl;

	return image;
}

void writeGreyPPM(string s,ImageGris img){

	Image imgRGB;
	int	i ;
	int j;
	int nLine;
	int nCol;

	nLine=img.size() ;
	nCol = img[0].size() ;

    imgRGB.resize(nLine) ;
	for (i=0 ; i <nLine; i++)
	  { imgRGB[i].resize(nCol) ;} ;

	for (i=0 ; i < nLine ; i++)
	  for (j=0 ; j< nCol ; j++)
	  {
		  imgRGB[i][j].r=imgRGB[i][j].g=imgRGB[i][j].b=img[i][j] ;
	  } ;
	writePPM(s,imgRGB) ;
	cout << "fait writeGreyImage" << endl;
}

void testReadWriteGreyImage(){
	cout << "verifier que vos images \"IO\" sont semblables a celles donnees en correction" << endl;
	writeGreyPPM("IO_Baboon.512.ppm",readGreyPPM("Baboon.512.ppm"));
	writeGreyPPM("IO_Billes.256.ppm",readGreyPPM("Billes.256.ppm"));
	writeGreyPPM("IO_Embryos.512.ppm",readGreyPPM("Embryos.512.ppm"));
	writeGreyPPM("IO_House.256.ppm",readGreyPPM("House.256.ppm"));
	writeGreyPPM("IO_Lena.512.ppm",readGreyPPM("Lena.512.ppm"));
}

ImageGris intensiteH(ImageGris img){
	int i, j ;
	int nLine, nCol ;
	ImageGris sortie;
	int iMin, iMax ;
	int jMin, jMax ;

	sortie = img ;

	nLine = img.size() ;
	nCol = img[0].size() ;


	for (i=0 ; i < nLine ; i++)
	{
		iMin= (i>1)? i-1 : 0 ;
		iMax= ((i+1)==nLine)? i : i+1 ;
		for (j=0 ; j < nCol ; j++)
		{

			jMin= (j>1)? j-1 : 0 ;
			jMax= ((j+1)==nCol)? j : j+1 ;
			sortie[i][j] = (img[iMin][jMin]+2*img[iMin][j]+img[iMin][jMax]) ;
			sortie[i][j] -= (img[iMax][jMin]+2*img[iMax][j]+img[iMax][jMax]) ;
			sortie[i][j] = fabs (sortie[i][j]) ;
		}
	}
	cout << "raf- intensiteH fait" << endl;

	return sortie;
}

ImageGris intensiteV(ImageGris img){
	int i, j ;
	int nLine, nCol ;
	int iMin, iMax ;
	int jMin, jMax ;
	ImageGris sortie;

    nLine = img.size() ;
	nCol = img[0].size() ;

	sortie = img ;

	for (i=0 ; i < nLine ; i++)
	{
		iMin= (i>1)? i-1 : 0 ;
		iMax= ((i+1)==nLine)? i : i+1 ;
		for (j=0 ; j < nCol ; j++)
		{

			jMin= (j>1)? j-1 : 0 ;
			jMax= ((j+1)==nCol)? j : j+1 ;
			sortie[i][j] = (img[iMin][jMin]+2*img[i][jMin]+img[iMax][jMin]) ;
			sortie[i][j] -= (img[iMin][jMax]+2*img[i][jMax]+img[iMax][jMax]) ;
		    sortie[i][j] = fabs (sortie[i][j]) ;
		}
	}
	cout << "raf- intensiteV fait" << endl;



	return sortie;
}

ImageGris intensite(ImageGris img){
	ImageGris  imgH ;
	ImageGris  imgV ;
	ImageGris sortie;
	int nLine ;
	int nCol ;
	int i,j ;

	imgH=intensiteV(img) ;
	imgV=intensiteV(img) ;

	nLine = img.size() ;
	nCol = img[0].size() ;


	sortie = img ;

	for (i=0 ; i < nCol ; i++)
	{
	   for (j=0 ; j < nLine ; j++)
	   {
	     sortie[i][j] = pow(imgH[i][j],2)+pow(imgV[i][j],2) ;
		 sortie[i][j] = sqrt(sortie[i][j]) ;
	   }
	}


	return sortie;
}

ImageGris scaleUp(ImageGris img){

	int nLine ;
	int nCol ;
	int i,j ;
	ImageGris sortie;
	double rScaleUp ;
	double maxGrey ;

	sortie=img ;

	nLine = img.size() ;
	nCol = img[0].size() ;

	for (maxGrey=1,i=0 ; i < nLine ; i++)
	{
		for (j=0; j < nCol ; j++)
		{
			maxGrey = (img[i][j] > maxGrey) ? img[i][j] : maxGrey ;
		}
	} ;

	rScaleUp = (maxGrey > 255)? 255 : 255/maxGrey ;

	for (i=0 ; i < nLine ; i++)
	{
		for (j=0; j < nCol ; j++)
		{
			sortie[i][j] *= rScaleUp ;
		}
	}


	cout << "raf - scaleup fait" << endl;
	return sortie;
}

void sauvegardeIntensite(string a,string b){
	writeGreyPPM(b,scaleUp(intensite(readGreyPPM(a))));
}

void testSobel(){
	cout << "verifier que les images \"sobel\" sont semblables a celles donnees en correction" << endl;
	sauvegardeIntensite("Baboon.512.ppm","sobel_Baboon.512.ppm");
	sauvegardeIntensite("Billes.256.ppm","sobel_Billes.256.ppm");
	sauvegardeIntensite("Embryos.512.ppm","sobel_Embryos.512.ppm");
	sauvegardeIntensite("House.256.ppm","sobel_House.256.ppm");
	sauvegardeIntensite("Lena.512.ppm","sobel_Lena.512.ppm");
}

ImageGris seuillage(ImageGris img, int seuil){

	int nLine ;
	int nCol ;
	int i,j ;
	ImageGris sortie;
	const int PIXEL_BLANC = 255 ;
	const int PIXEL_NOIR = 0 ;

	sortie=img ;

	nLine = img.size() ;
	nCol = img[0].size() ;


	for (i=0 ; i < nLine ; i++)
	{
		for (j=0; j < nCol ; j++)
		{
			sortie[i][j] = (img[i][j] >= seuil) ? PIXEL_NOIR : PIXEL_BLANC ;
		}
	}


	cout << "Raf: seuillage fait" << endl;
	return sortie;
}

void sauvegardeSeuillage(string a, string b, int seuil){
	writeGreyPPM(b,seuillage(readGreyPPM(a),seuil));
}

void testSeuillage(){
	cout << "verifier que les images \"seuillage\" sont semblables a celles donnees en correction" << endl;
	sauvegardeSeuillage("Billes.256.ppm","seuillage_Billes.256.ppm",80);
	sauvegardeSeuillage("Lena.512.ppm","seuillage_Lena.512.ppm",80);

	cout << "\tproposer des seuils pour Embryos.512.ppm - House.256.ppm" << endl;
}

ImageGris seuillageBis(ImageGris imgIntensite, int seuil){


	cout << "raf: seuillage bis fait" << endl;

	return seuillage ( intensite(imgIntensite), seuil);
}



ImageGris doubleSeuillage(ImageGris imgIntensite, ImageGris imgContour, int seuil){

	int nLine ;
	int nCol ;
	int i,j ;
	int iContour, jContour ;
	int iContourMin, iContourMax ;
	int jContourMin, jContourMax;
	bool bProche ;

	ImageGris sortie;
	const int PIXEL_BLANC = 255 ;
	const int PIXEL_NOIR = 0 ;


	sortie=seuillageBis(imgIntensite,seuil) ;

	nLine = imgIntensite.size() ;
	nCol = imgIntensite[0].size() ;


	for (i=0 ; i < nLine ; i++)
	{
		iContourMin=(i>1)?i-1:0 ;
		iContourMax=((i+1)==nLine)? i : i+1 ;

		for (j=0; j < nCol ; j++)
		{   if (sortie[i][j] == PIXEL_NOIR)
            {
				jContourMax=(j>1)?j-1:0 ;
				jContourMin=((j+1)==nCol)? j: j+1 ;

		       for (bProche=false, iContour=iContourMin; iContour <= iContourMax ; iContour++)
		       {
			     for (jContour=jContourMin ; jContour <= jContourMax ; jContour )
			     {
					if (!imgContour[iContour][jContour]) bProche = true ;
			     }
			   }
			   sortie[i][j]=(bProche)? PIXEL_NOIR : PIXEL_BLANC ;
		   }

		}
	}


	cout << "Raf : doubleSeuillage - Contour fait" << endl;
	return sortie;
}



typedef struct { int i ; int j; } coord ;

bool aUnIntenseProche(int nbAmelioration, coord ptRef, ImageGris &imgIntensite)
{
    int nLine;
	int	nCol;
	int i, j ;
	int iMin, iMax ;
	int jMin, jMax ;
	coord pt ;
	const int PIXEL_BLANC = 255 ;
	const int PIXEL_NOIR = 0 ;

	if (imgIntensite[ptRef.i][ptRef.j] == PIXEL_NOIR)
		return (true) ;

	if (!nbAmelioration)
		return (false);
	--nbAmelioration;

	nLine = imgIntensite.size() ;
	nCol = imgIntensite[0].size() ;

	iMin=(ptRef.i>1)?ptRef.i-1:0 ;
	iMax=((ptRef.i+1)==nLine)? ptRef.i : ptRef.i+1 ;

	jMax=(ptRef.j>1)?ptRef.j-1:0 ;
	jMin=((ptRef.j+1)==nCol)? ptRef.j: ptRef.j+1 ;

	for (i=iMin ; i<= iMax ; i++)
		for (j=jMin ; j<= jMax ; j++)
		{
			pt.i = i ;
		    pt.j = j ;
			if (aUnIntenseProche(nbAmelioration, pt, imgIntensite) )
				return (true) ;
		}
  return (false) ;
}

ImageGris doubleSeuillage(ImageGris img, int seuilFort, int seuilFaible, int nbAmelioration){

    int nLine ;
	int nCol ;
	int i,j ;
	coord ptRef ;
	ImageGris sortie;
	ImageGris imgIntensite ;
	const int PIXEL_BLANC = 255 ;
	const int PIXEL_NOIR = 0 ;


	imgIntensite=seuillageBis(img,seuilFort) ;
	sortie=imgIntensite ;

	nLine = img.size() ;
	nCol = img[0].size() ;


	for (i=0 ; i < nLine ; i++)
	{
		for (j=0; j < nCol ; j++)
		{
			if ((imgIntensite[i][j] == PIXEL_BLANC) && (img[i][j] >= seuilFaible))
			{

			   ptRef.i= i ;
			   ptRef.j = j ;
			   sortie[i][j] = (aUnIntenseProche(nbAmelioration,ptRef,imgIntensite))? PIXEL_NOIR : PIXEL_BLANC ;
			}
		}
	}
	cout << "Raf : doubleSeuillage - reste à faire" << endl;

	return sortie;
}

void sauvegardeSeuillageDouble(string a, string b, int seuilFort, int seuilFaible, int nbAmelioration){
	writeGreyPPM(b,doubleSeuillage(readGreyPPM(a),seuilFort,seuilFaible,nbAmelioration));
}

void testDoubleSeuillage(){
	cout << "verifier que les images \"seuillage_double\" sont semblables a celles donnees en correction" << endl;
	sauvegardeSeuillageDouble("Billes.256.ppm","seuillage_double_Billes.256.ppm",500,80,100);
	sauvegardeSeuillageDouble("Lena.512.ppm","seuillage_double_Lena.512.ppm",500,80,100);

	cout << "\tproposer des seuils pour Embryos.512.ppm - House.256.ppm" << endl;
}

typedef vector<double> Point;

typedef vector<vector<double> > EnsPoint;

double euclidianDistance(Point p,Point c){

            int dimR ;
            int i ;
            double distance ;


            if (p.size() != c.size())
             return (0) ;

            dimR = p.size()  ;
            for (i=0, distance=0 ; i<dimR ;i++)
            {
              distance += pow((p[i]-c[i]),2);
            }
            distance = sqrt ( distance) ;

            cout << "à tester &  à compiler" << endl;
            return distance;
}


double distanceAEnsemble(Point p,EnsPoint C){
	double dMin =0 ;
	double distance ;


  if (!C.size())
    return 0 ;

  for (int n=0; n < (int) C.size() ; n++)
  {
    distance =  euclidianDistance(p,C[n]) ;
    dMin= (distance > dMin)? dMin : distance ;
  }
  cout << "à tester & à compiler" << endl;
  return dMin ;
}

int plusProcheVoisin(Point p,EnsPoint C){
	double dMin ;

  dMin = distanceAEnsemble (p, C) ;
  for (int n=0 ; n< (int) C.size() ; n++)
  {
     if ( euclidianDistance(p, C[n]) <= dMin )
       return (n) ;
  }


  return -1;
}

EnsPoint sousEnsemble(EnsPoint P,EnsPoint C,int k){
	EnsPoint sortie;


  for (int n=0 ; n <(int)  P.size() ; n++)
  {
     if (k==plusProcheVoisin(P[n],C))
        sortie.push_back(P[n]) ;
  }

  cout << "ce qui est fait n'est plus à faire mais il faut tout de même compiler & tester " << endl;

  return sortie;
}

Point barycentre(EnsPoint Q){
	Point sortie;



  sortie = Q[0] ;

  for (int n=1; n  < (int) Q.size() ; n++)
  {
    for (int i=0 ; i< (int) sortie.size() ; i++)
       sortie[i]+=Q[n][i] ;
  }
  for (int i=0; i< (int) sortie.size(); i++)
     sortie[i] /= Q.size() ;

  return sortie;
}

EnsPoint Kmoyenne(EnsPoint P,EnsPoint C, int nbAmeliorations){
	cout << "a faire" << endl;
	return C;
}

EnsPoint FAST_Kmoyenne(EnsPoint P,EnsPoint C, int nbAmeliorations)
{
  vector<int> label;

  label.resize(P.size());
  for(int n=0;n<nbAmeliorations;n++)
  {
     vector<int> clusterSize;
     clusterSize.resize(C.size(),0);
     for (int p=((int)P.size())-1; p>=0 ;p--)
     {  double di = 0;
        int nn=0;
        for(int d=((int)P[0].size())-1;d>=0;d--)
          di+=(P[p][d]-C[0][d])*(P[p][d]-C[0][d]);
        for(int c=((int)C.size())-1;c>=1;c--)
        {  double dt=0;
           for(int d=((int)P[0].size())-1;d>=0;d--)
             dt+=(P[p][d]-C[c][d])*(P[p][d]-C[c][d]);

           if(dt<di)
             {di=dt;nn=c;}
        }
        label[p]=nn;
        clusterSize[nn]++;
     }
     for (int p=((int)P.size())-1;p>=0;p--)
       for(int d=((int)P[0].size())-1;d>=0;d--)
         C[label[p]][d]+=P[p][d];
     for(int c=((int)C.size())-1;c>=0;c--)
       if(clusterSize[c]!=0)
         for(int d=((int)P[0].size())-1;d>=0;d--)
           C[c][d] = C[c][d]/(clusterSize[c]+1);
  }
  return C;
}

void testKmoyenne(){
	EnsPoint P;
	for(int i=0;i<100;i++){
		double xbiais=0.;
		double ybiais=0.;
		if(i%3==1){
			xbiais = 100.;
			ybiais = 0.;
		}
		if(i%3==2){
			xbiais = 0.;
			ybiais = 100.;
		}
		Point p;
		p.push_back(xbiais + 1.*(rand()%10)-20.);
		p.push_back(ybiais + 1.*(rand()%10)-20.);
		P.push_back(p);
	}

	EnsPoint C;
	if(true){
		Point c;
		c.push_back(75.);
		c.push_back(50.);
		C.push_back(c);
	}
	if(true){
		Point c;
		c.push_back(50.);
		c.push_back(75.);
		C.push_back(c);
	}
	if(true){
		Point c;

		c.push_back(25.);
		c.push_back(25.);
		C.push_back(c);
	}

	EnsPoint Cout = Kmoyenne(P,C,100);
	EnsPoint CRef = FAST_Kmoyenne(P,C,100);

	cout << "verifier que les valeurs sont comparables : " << endl;
	cout << "C[0] : \t\tKmoyenne \t\tFAST_Kmoyenne" << endl;
	for(int j=0;j<2;j++)
		cout << "\t\t" << Cout[0][j] << "\t\t" <<CRef[0][j] << endl;
	cout << endl;

	cout << "C[1] : \t\tKmoyenne \t\tFAST_Kmoyenne" << endl;
	for(int j=0;j<2;j++)
		cout << "\t\t" << Cout[1][j] << "\t\t" <<CRef[1][j] << endl;
	cout << endl;

	cout << "C[2] : \t\tKmoyenne \t\tFAST_Kmoyenne" << endl;
	for(int j=0;j<2;j++)
		cout << "\t\t" << Cout[2][j] << "\t\t" <<CRef[2][j] << endl;
	cout << endl;
}




EnsPoint pivotSuperPixel(ImageGris img,double lambda, int mu){
	cout << "a faire" << endl;
	EnsPoint sortie;
	return sortie;
}

EnsPoint superPixel(ImageGris img,double lambda, int mu, int nbAmeliorations){
	cout << "a faire et penser a utiliser FAST_Kmoyenne !" << endl;
	EnsPoint sortie;
	return sortie;
}

ImageGris appliqueSuperPixel(ImageGris img,double lambda, int mu, int nbAmeliorations){
	cout << "a faire" << endl;
	return img;
}

void savegardeSuperPixel(string a,string b,double lambda, int mu, int nbAmeliorations){
	writeGreyPPM(b,appliqueSuperPixel(readGreyPPM(a),lambda,mu,nbAmeliorations));
}

void testSuperPixel(){
	cout << "verifier que les images \"superpixel\" sont semblables a celles donnees en correction" << endl;
	savegardeSuperPixel("Billes.256.ppm","superpixel_Billes.256.ppm",1,30,25);
	savegardeSuperPixel("Baboon.512.ppm","superpixel_Baboon.512.ppm",2,60,15);
	savegardeSuperPixel("Lena.512.ppm","superpixel_Lena.512.ppm",2,60,15);

	cout << "\tproposer des parametres pour Embryos.512.ppm - House.256.ppm" << endl;
}


int RFUmain(int, char*[]) {

	cout << "test main" << endl;
	testReadWriteImage();

    testReadWriteGreyImage();

    testSobel();
    testSeuillage();
    testDoubleSeuillage();
    //testKmoyenne();
    //testSuperPixel();
    return 0;
}


int main(int,char*[])
{


	testReadWriteImage();

    testReadWriteGreyImage();

    testSobel();
    testSeuillage();
    testDoubleSeuillage();
    //testKmoyenne();
    //testSuperPixel();
	return 0;
}
