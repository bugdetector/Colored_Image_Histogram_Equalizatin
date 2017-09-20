#include<stdio.h>
#include<conio.h>
#include<stdlib.h>
#include<math.h>
#include<ctype.h>
#include<string.h>
#include <fcntl.h>
#include <malloc.h>
#include <math.h>
#define PI 3.1415926535897932384626433832795
struct ppm_header
{
	char pgmtype1;
	char pgmtype2;
	int pwidth;
	int pheight;
	int pmax;
};
struct ppm_file
{
	struct ppm_header *pheader;
	unsigned char *rdata,*gdata,*bdata;
};
double min(double x,double y, double z){
    if(x<=y && x<=z){
        return x;
    }else if(y<=x && y<=z){
        return y;
    }else{
        return z;
    }
}
void Histogram_Equilization(double *array,int size);

void RGB2HSI(struct ppm_file *image,double *Hue,double *Saturation,double *Intensity);//http://answers.opencv.org/question/62446/conversion-from-rgb-to-hsi/
void HSI2RGB(struct ppm_file *image,double *Hue,double *Saturation,double *Intensity);//https://gist.github.com/rzhukov/9129585

void RGB2YCbCr(struct ppm_file *image,double *Y,double *Cb,double *Cr);//http://www.mir.com/DMG/ycbcr.html
void YCbCr2RGB(struct ppm_file *image,double *Y,double *Cb,double *Cr);

void get_image_data(char *filename,struct ppm_file *image);
void write_image(char *filename,struct ppm_file *image);

int main()
{
	struct ppm_file resim;
	get_image_data("quant1.ppm",&resim);

    int totalPixel = resim.pheader->pwidth*resim.pheader->pheight;
    double *Y = (double*)malloc(totalPixel*sizeof(double));// Cb ve Cr değerleri kesirli değerler olduğu için veri tipi doubledır.
    double *Cb = (double*)malloc(totalPixel*sizeof(double));
    double *Cr = (double*)malloc(totalPixel*sizeof(double));

    double *H = (double*)malloc(totalPixel*sizeof(double));//H ve S kesirli değerler
    double *S = (double*)malloc(totalPixel*sizeof(double));
    double *I = (double*)malloc(totalPixel*sizeof(double));
    RGB2YCbCr(&resim,Y,Cb,Cr);
    RGB2HSI(&resim,H,S,I);
    Histogram_Equilization(Y,totalPixel);//Y ve I için histogram eşleme
    Histogram_Equilization(I,totalPixel);

    struct ppm_file yeniresim = resim;
    YCbCr2RGB(&yeniresim,Y,Cb,Cr);
    write_image("YCbCr.ppm",&yeniresim);
    HSI2RGB(&yeniresim,H,S,I);
    write_image("HSI.ppm",&yeniresim);
	/*printf("pgmtype...=%c%c\n",resim.pheader->pgmtype1,resim.pheader->pgmtype2);
	printf("width...=%d\n",resim.pheader->pwidth);
	printf("height...=%d\n",resim.pheader->pheight);        Buralar hocanýn yazdýðý kodlar.
	printf("max gray level...=%d\n",resim.pheader->pmax);*/


	return 0;
}         //https://www.youtube.com/watch?v=eNBZI-qYhpg&t=178s
void Histogram_Equilization(double *array,int size){//Verilen diziye histogram eşleme yöntemi uygular. Dizinin boyutu gereklidir.
    int *frequencyCounter = (int*)calloc(255,sizeof(int));// Her bir pixelden kaç adet bulunduðunu sayarak saklar.
    double *probabilities = (double*)calloc(255,sizeof(double));// Her bir pixel deðerinin bulunma olasılığını tutar.
    int i;//https://github.com/bugdetector/Image_Processing/blob/master/Histogram_Equilization.c buradaki fonksiyonu değiştirdim.
    for(i=0;i<size;i++){//Pixeller sayılıyor...
        frequencyCounter[(unsigned char)array[i]]++;
    }
    for(i=0;i<256;i++){//Pixellerin bulunma olasılıkları hesaplanıyor.
        probabilities[i] = ((double)frequencyCounter[i])/((double)size);
    }
    for(i=1;i<256;i++){//Olasilikların artirimli (kümülatif) toplamlari hesaplaniyor.
        probabilities[i] += probabilities[i-1];
    }
    for(i=0;i<size;i++){// Yeni pixel deðerleri eşleniyor.
        array[i] = probabilities[(unsigned char)array[i]]*255;
    }
}
void RGB2HSI(struct ppm_file *image,double *Hue,double *Saturation,double *Intensity){
    int i;
    double R,G,B,H,S,I;
    double total;
    int totalPixel = image->pheader->pheight*image->pheader->pwidth;//https://gist.github.com/rzhukov/9129585 buradan alıp değiştirdim.
    for(i=0;i<totalPixel;i++){
        R = image->rdata[i];
        G = image->gdata[i];
        B = image->bdata[i];
        total = R+G+B;
        I = total / 3.0;
        double rn = R / total;
        double gn = G / total;
        double bn = B / total;

        H = acos((0.5 * ((rn - gn) + (rn - bn))) / (sqrt((rn - gn) * (rn - gn) + (rn - bn) * (gn - bn))));
        if(B > G)
        {
            H = 2 * PI - H;
        }

        S = 1 - 3 * min(rn , gn, bn);

        Hue[i] = H;
        Saturation[i] = S;
        Intensity[i] = I;
    }
}
void HSI2RGB(struct ppm_file *image,double *Hue,double *Saturation,double *Intensity){
    int totalPixels = image->pheader->pheight*image->pheader->pwidth;//toplam pixel sayısı
    int i;
    double H,S,I,R,G,B;
    for(i=0;i<totalPixels;i++){ //https://gist.github.com/rzhukov/9129585 buradan alıp değiştirdiðim bir kod.
        H = Hue[i];
        S = Saturation[i];
        I = Intensity[i];
        double x = I * (1 - S);
        if(H < 2 * PI / 3)
        {
            double y = I * (1 + (S * cos(H)) / (cos(PI / 3 - H)));
            double z = 3 * I - (x + y);
            B = x; R = y; G = z;
        }
        else if(H < 4 * PI / 3)
        {
            double y = I * (1 + (S * cos(H - 2 * PI / 3)) / (cos(PI / 3 - (H - 2 * PI / 3))));
            double z = 3 * I - (x + y);
            R = x; G = y; B = z;
        }
        else
        {
            double y = I * (1 + (S * cos(H - 4 * PI / 3)) / (cos(PI / 3 - (H - 4 * PI / 3))));
            double z = 3 * I - (x + y);
            R = z; G = x; B = y;
        }
        image->rdata[i] = (unsigned char)(R<0 ? 0 : (R>255 ? 255 : R));
        image->gdata[i] = (unsigned char)(G<0 ? 0 : (G>255 ? 255 : G));
        image->bdata[i] = (unsigned char)(B<0 ? 0 : (B>255 ? 255 : B));
    }
    return;
}
void RGB2YCbCr(struct ppm_file *image,double *Y,double *Cb,double *Cr){
    int totalPixels = image->pheader->pheight*image->pheader->pwidth;
    int i;
    double R,G,B;
    for(i=0;i<totalPixels;i++){// http://answers.opencv.org/question/62446/conversion-from-rgb-to-hsi/ bu sitedeki formülleri uyguladım.
        R = image->rdata[i];
        G = image->gdata[i];
        B = image->bdata[i];
        Y[i] = 0.299 * R + 0.587 *G	+ 0.114 * B;
        Cb[i] = -0.168736 *R - 0.331264 * G + 0.500 * B;
        Cr[i] = 0.500 * R - 0.418688 * G - 0.081312 * B;
    }
    return;
}
void YCbCr2RGB(struct ppm_file *image,double *Y,double *Cb,double *Cr){
    int totalPixels = image->pheader->pheight*image->pheader->pwidth;
    int i;
    double R,G,B;
    for(i=0;i<totalPixels;i++){// http://answers.opencv.org/question/62446/conversion-from-rgb-to-hsi/ bu site.
        R = (double) (1.0 * Y[i]+ 0 * Cb[i]	+ 1.402 * Cr[i] );
        G = (double) (1.0 * Y[i]-0.344136 * Cb[i]-0.714136 * Cr[i]);
        B = (double) (1.0 * Y[i] + 1.772 * Cb[i] + 0 * Cr[i]);

        image->rdata[i] = (unsigned char)(R>255 ? 255 : (R<0 ? 0 : R));
        image->gdata[i] = (unsigned char)(G>255 ? 255 : (G<0 ? 0 : G));
        image->bdata[i] = (unsigned char)(B>255 ? 255 : (B<0 ? 0 : B));
    }
    return;
}
void write_image(char *filename,struct ppm_file *image)
{
	FILE *fp;
	int i,max=0;
	fp=fopen(filename,"wb");
	fputc(image->pheader->pgmtype1,fp);
	fputc(image->pheader->pgmtype2,fp);
	fputc('\n',fp);
	fprintf(fp,"%d %d\n",image->pheader->pwidth,image->pheader->pheight);
	fprintf(fp,"%d\n",255/*max*/);
	for(i=0;i<image->pheader->pwidth*image->pheader->pheight;i++)
	{
		fwrite(&image->rdata[i],1,1,fp);
		fwrite(&image->gdata[i],1,1,fp);
		fwrite(&image->bdata[i],1,1,fp);
	}
	fclose(fp);
}
void get_image_data(char *filename, struct ppm_file *image )
{
	FILE* fp;
	int i=0;
	char temp[256];
	image->pheader=(struct ppm_header *)malloc(sizeof(struct ppm_header));
	fp = fopen(filename, "rb" );
	if (fp==NULL)
	{
		printf("Dosya acilamadi: %s.\n\n", filename);
		exit(1);
	}
	printf ("Okunan PPM dosyasi : %s...\n", filename);
	fscanf (fp, "%s", temp);
	if (strcmp(temp, "P6") == 0)
	{
		image->pheader->pgmtype1=temp[0];
		image->pheader->pgmtype2=temp[1];
		fscanf (fp, "%s", temp);
		if (temp[0]=='#')
		{

			while(fgetc(fp)!='\n');
			fscanf (fp, "%d %d\n",&image->pheader->pwidth,&image->pheader->pheight);
			fscanf (fp, "%d\n", &image->pheader->pmax);

		}
		else
		{
			sscanf (temp, "%d", &image->pheader->pwidth);
			fscanf (fp, "%d", &image->pheader->pheight);
			fscanf (fp, "%d", &image->pheader->pmax);
		}
		image->rdata=(unsigned char *)malloc(image->pheader->pheight*image->pheader->pwidth*sizeof(unsigned char));
		image->gdata=(unsigned char *)malloc(image->pheader->pheight*image->pheader->pwidth*sizeof(unsigned char));
		image->bdata=(unsigned char *)malloc(image->pheader->pheight*image->pheader->pwidth*sizeof(unsigned char));
		if (image->rdata==NULL) printf("bellek problemi...\n");
		for(i=0;i<image->pheader->pwidth*image->pheader->pheight;i++)
		{
			fread(&image->rdata[i],1,1,fp);
			fread(&image->gdata[i],1,1,fp);
			fread(&image->bdata[i],1,1,fp);
		}
	}
	else
	{
		printf ("\nHata Resim dosyasi PGM P6 formatinda degil");
		exit(1);
	}
	fclose(fp);
}
