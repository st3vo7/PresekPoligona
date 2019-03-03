/*
 * Pri zadavanju temena ova dva poligona
 * njihov redosled treba da bude ciklican,
 * ali i dat u redosledu obilaska, dakle,
 * u pozitivnom matematickom smeru
 * (suprotno od kretanja kazaljke na casovniku)
 * */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define X 0
#define Y 1
typedef enum {FALSE, TRUE} bool;
typedef enum {Pun, Qun, Nepoznato} tUnFlag;

#define DIM 2
#define TMAX 1000
typedef int tTackai[DIM];
typedef double tTackad[DIM];

typedef tTackai tPoligoni[TMAX];

//globalne promenljive
int n, m;
tPoligoni P, Q;
FILE *ispis; //bice fajl za ispis preseka

//definicije funkcija
int UcitajPoligon(tPoligoni P);
void OdrediPresek(tPoligoni P, tPoligoni Q, int n, int m);
int Nacrtaj();
void Vektor(tTackai a, tTackai b, tTackai c);
int Znak( tTackai a, tTackai b, tTackai c );
char Seci(tTackai p, tTackai b, tTackai q, tTackai d, tTackad pr1, tTackad pr2);
int proizvodX(tTackai v, tTackai w);
int proizvodD(tTackai v, tTackai w);
tUnFlag PostaviFleg(tTackad p, tUnFlag unflag, int aHB, int bHA);
int Napreduj(int a, int *aa, int n, bool unutra, tTackai v);
double duzina_duzi(tTackai a, tTackai b);
int jednacina_prave(tTackai a, tTackai p, tTackai q);


int main(){
    
    
    
    clock_t t;
    
    n=UcitajPoligon(P);
    m=UcitajPoligon(Q);
    
    t=clock();
    
    OdrediPresek(P,Q,n,m);
    
    t=clock()-t;
    double vreme_izvrsavanja = ((double)t)/CLOCKS_PER_SEC;
    
    printf("Presek je odredjen u %f sekundi\n",vreme_izvrsavanja);
    
    Nacrtaj();
    
 
    exit(EXIT_SUCCESS);
}

int UcitajPoligon(tPoligoni P){
    
    /*ucitavam poligone a istovremeno upisujem
     * u fajl poligoni.txt kao parsirane za iscrtavanje pomocu gnuplota
     */
    
    FILE *fptr; //fajl za stampanje
    
    int n=0;
    int nin;
    
    if((fptr=fopen("poligoni.txt","a"))==NULL){
        printf("Greska pri otvaranju fajla\n");
        exit(EXIT_FAILURE);
    }
    
    
    tPoligoni Pom;  //pomocni, radi ispisa, zato sto gnuplot -,-'
    
    
    printf("Unesi broj temena:\n");
    scanf("%d", &nin);
    printf("%%Poligon:\n");
    printf("%% i    x    y\n");
    
    scanf("%d %d", &Pom[0][X], &Pom[0][Y]);
    fprintf(fptr, "%d %d\n",Pom[0][X],Pom[0][Y]);
    printf("%%%3d%4d%4d\n",n,Pom[0][X],Pom[0][Y]);
    P[0][X]=Pom[0][X];
    P[0][Y]=Pom[0][Y];
    
    ++n;
    while((n<nin)&&(scanf("%d %d", &P[n][X], &P[n][Y])!=EOF)){
        printf("%%%3d%4d%4d\n",n,P[n][X],P[n][Y]);
        fprintf(fptr, "%d %d\n",P[n][X],P[n][Y]);
        ++n;
    }
    fprintf(fptr, "%d %d\n\n\n",Pom[0][X],Pom[0][Y]);
    if(n>=TMAX)
        printf("Greska u UcitajPoligon. Previse temena. Maksimum je %d\n",TMAX);
    
    fclose(fptr);
    return n;
    
}

void OdrediPresek(tPoligoni P, tPoligoni Q, int n, int m){
    
    //funkcija dobija poligone P i Q sa po n i m temena, redom
    
    int a, b;       //indeksi temena poligona P, Q
    int a1, b1;     //a-1, b-1
    tTackai A, B;   //usmerene ivice na P, Q
    int presek;     //znak rezultujuceg vektora AxB
    int bHA, aHB;   //b pripada H(A) tj. a pripada H(B)
    tTackai Pocetak = {0,0};
    tTackad p;      //tacka preseka (koordinate tipa double)
    tTackad q;      //druga tacka preseka (ako postoji preklapanje)
    tUnFlag unflag; //{Pun, Qun, Nepoznato}; koji je unutra
    int aa, ba;     //broj napredaka po a i b indeksima (nakon prvog preseka)
    bool PrvaTacka; //kod inicijalizacije
    tTackad p0;     //prva tacka
    int rez;        //ono sto vrati f-ja Seci
    
    a = 0;
    b = 0;
    aa = 0;
    ba = 0;
    unflag = Nepoznato;
    PrvaTacka = TRUE;
    
    
    
    if((ispis=fopen("resenje.txt","w"))==NULL){
        printf("Greska pri otvaranju fajla\n");
        exit(EXIT_FAILURE);
    }
    
    do{
        printf("Pre napredovanja: a=%d aa=%d b=%d ba=%d unflag=%d\n",a,aa,b,ba,unflag);
        a1 = (a+n-1)%n;
        b1 = (b+m-1)%m;
        
       // printf("a1 je: %d\n",a1);
       // printf("b1 je: %d\n",b1);
        
       // printf("P[a]: {%d,%d}\n",P[a][X],P[a][Y]);
       // printf("P[a1]: {%d,%d}\n",P[a1][X],P[a1][Y]);
        
       // printf("Q[b]: {%d,%d}\n",Q[b][X],Q[b][Y]);
       // printf("Q[b1]: {%d,%d}\n",Q[b1][X],Q[b1][Y]);
        
        Vektor(P[a],P[a1],A);   //napravi usmerene vektore ivica
        Vektor(Q[b],Q[b1],B);
        
        printf("A[X]: %d A[Y]: %d\n",A[X],A[Y]);
        printf("B[X]: %d B[Y]: %d\n",B[X],B[Y]);
        
        presek = Znak(Pocetak, A, B);
        aHB = Znak(Q[b1],Q[b],P[a]);
        bHA = Znak(P[a1],P[a],Q[b]);
        
        //A i B su paralelni <=> presek = 0 linija do
        
       printf("presek: %d aHB: %d bHA: %d\n",presek,aHB,bHA);
        
        rez = Seci(P[a1],P[a],Q[b1],Q[b],p,q);
        printf("Rezultat Secija: %c\n",rez);
        
        if(rez=='1' || rez =='v'){
            //presek je jedna tacka
            if(unflag == Nepoznato && PrvaTacka){
                aa = ba = 0;
                PrvaTacka = FALSE;
                p0[X] = p[X];
                p0[Y] = p[Y];
                
                printf("------------{%lf,%lf}-----------\n",p0[X],p0[Y]);
                fprintf(ispis, "%lf %lf\n",p0[X],p0[Y]);
            }
            unflag = PostaviFleg(p,unflag,aHB,bHA);
            printf("PostaviFleg je postavio unflag na: %d\n",unflag);
        }
        
        
        
        
        
        //A i B se preklapaju i suprotno orijentisani
        //dakle presek poligona = presek tih duzi
        printf("proizvodD(A,B) = %d\n",proizvodD(A,B));
        if((rez=='e') && (proizvodD(A,B)<0)){
            printf("presek ova dva poligona je duz ");
            printf(" [%lf,%lf]-----------",p[X],p[Y]);
            printf(" [%lf,%lf]-----------\n",q[X],q[Y]);
            return;
            
        }
        
       
        
        //A i B su paralelni i razliciti
        else if((presek==0)&&(aHB<0)&&(bHA<0)){
            printf("Granice P i Q se ne seku\n");
            return;
        }
        
        //-----Pravila za napredovanje-----
        
        //A i B su kolinearni i razliciti *----A-----*   *-----B------*
        else if((presek==0)&&(aHB==0)&&(bHA==0)){
            if(unflag==Pun)
                b = Napreduj(b,&ba,m,unflag==Qun,Q[b]);
            else
                a = Napreduj(a,&aa,m,unflag==Pun,P[a]);
        }
        
        else if (presek>=0){
            if(bHA>0)
                a=Napreduj(a,&aa,n,unflag==Pun,P[a]);
            else
                b=Napreduj(b,&ba,m,unflag==Qun,Q[b]);
        }
        
        else{
            //presek<0
            if(aHB>0)
                b=Napreduj(b,&ba,m,unflag==Qun,Q[b]);
            else
                a=Napreduj(a,&aa,n,unflag==Pun,P[a]);
            
        }
        
        
        printf("Nakon napredovanja: a=%d aa=%d b=%d ba=%d unflag=%d\n",a,aa,b,ba,unflag);
        
        //obradi specijalne slucajeve :(
        
    
    } while ( ((aa<n) || (ba < m)) && (aa<2*n) && (ba<2*m));
    
    
    if(!PrvaTacka){
        printf("++++++++++++++++++++++++++++++++++%lf %lf linija do\n",p0[X], p0[Y]);
        fprintf(ispis, "%lf %lf\n",p0[X],p0[Y]);
        
    }
    
    if(unflag == Nepoznato){
        printf("Granice P i Q se ne seku\n");
        
    }
    
    fclose(ispis);
    
}

/*f-ja seci pronalazi tacku preseka izmedju dva zatvorena segmenta pb i qd
Vraca pr1 i karakter sa sledecim znacenjem:
'e': segmenti se kolinearno preklapaju
'v': kraj jednog segmenta (teme) pripada drugoj stranici, ali nisu kolinearne
'1': segmenti se "lepo" seku (dakle dele tacku i ne vaze ni 'v' ni 'e')
'0': segmenti se ne seku
Napomena: slucaj kolinearnih segmenata koji dele teme vraca rez 'e'*/

char Seci(tTackai p, tTackai b, tTackai q, tTackai d, tTackad pr1, tTackad pr2){
    
    // ovde ce mozda biti problema
    tTackai r = {b[X]-p[X],b[Y]-p[Y]};
    tTackai s = {d[X]-q[X],d[Y]-q[Y]};
    
    
    tTackai q_minus_p = {q[X]-p[X],q[Y]-p[Y]};
    
    double imenilac = proizvodX(r,s);
    double brojilac1 = proizvodX(q_minus_p,s);
    double brojilac2 = proizvodX(q_minus_p,r);
    
    double t = brojilac1/imenilac;
    double u = brojilac2/imenilac;
    
    printf("imenilac = %lf\n",imenilac);
    printf("brojilac1 = %lf\n",brojilac1);
    printf("brojilac2 = %lf\n",brojilac2);
    
    if(imenilac==0 && brojilac2==0){
        //kolinearni su
        printf("kolinearni su\n");
        
        double t0 = (double) proizvodD(q_minus_p,r) / (double) proizvodD(r,r);
        double t1 = t0 + (double) proizvodD(s,r) / (double) proizvodD(r,r);
        printf("t0: %lf, t1: %lf\n",t0,t1);
        
        double temin, temaks;
        
        if(t0<t1){
            temin = t0;
            temaks = t1;
        }
        else {
            temin = t1;
            temaks = t0;
        }
        
        if(temaks < 0 || temin > 1){
            printf("ne postoji presek\n");
            return '0';
        }
        else{
            printf("jedan je sadrzan u drugome, ili se dodiruju u tacki\n");
            
            //sad bi ovde trebalo da vidis koji je slucaj i da prosledis odgovarajuce pr1 i pr2
            if(t0==0){
                if(t1>0){
                    //duz
                    pr1[X]=p[X];
                    pr1[Y]=p[Y];
                    
                    pr2[X]=b[X]<d[X]?b[X]:d[X];
                    pr2[Y]=b[Y]<d[Y]?b[Y]:d[Y];
                }
                else if(t1<0){
                    //tacka
                    pr1[X]=p[X];
                    pr1[Y]=p[Y];
                    
                    pr2[X]=p[X];
                    pr2[Y]=p[Y];
                    
                }
            }
            if(t1==0){
                if(t0>0){
                    //duz
                    if(duzina_duzi(p,b)<duzina_duzi(q,d)){
                        printf("Duzina duzi (p,b): %lf\n",duzina_duzi(p,b));
                        printf("Duzina duzi (q,d): %lf\n",duzina_duzi(q,d));
                        pr1[X]=p[X];
                        pr1[Y]=p[Y];
                        pr2[X]=b[X];
                        pr2[Y]=b[Y];
                        
                    }
                    else{
                        pr1[X]=q[X];
                        pr1[Y]=q[Y];
                        pr2[X]=d[X];
                        pr2[Y]=d[Y];
                    }
                    
                }
                else if(t0<0){
                    //tacka
                    if((p[X]==d[X])&&(p[Y]==d[Y])){
                        pr1[X]=p[X];
                        pr1[Y]=p[Y];
                    
                        pr2[X]=p[X];
                        pr2[Y]=p[Y];
                    }
                    else{
                        pr1[X]=q[X];
                        pr1[Y]=q[Y];
                    
                        pr2[X]=b[X];
                        pr2[Y]=b[Y];
                        
                    }
                        
                }
            }
            if(t0==1){
                if(t1<1){
                    
                    //duz
                    if(duzina_duzi(p,b)<duzina_duzi(q,d)){
                        pr1[X]=p[X];
                        pr1[Y]=p[Y];
                        pr2[X]=b[X];
                        pr2[Y]=b[Y];
                        
                    }
                    else{
                        pr1[X]=q[X];
                        pr1[Y]=q[Y];
                        pr2[X]=d[X];
                        pr2[Y]=d[Y];
                    }
      
                }
                else if(t1>1){
                    //tacka
                    if((p[X]==d[X])&&(p[Y]==d[Y])){
                        pr1[X]=p[X];
                        pr1[Y]=p[Y];
                    
                        pr2[X]=p[X];
                        pr2[Y]=p[Y];
                    }
                    else{
                        pr1[X]=q[X];
                        pr1[Y]=q[Y];
                    
                        pr2[X]=b[X];
                        pr2[Y]=b[Y];
                        
                    }
                    
                }
            }
            if(t1==1){
                if(t0>1){
                    //tacka
                    printf("tacka:\n");
                        pr1[X]=b[X];
                        pr1[Y]=b[Y];
                    
                        pr2[X]=b[X];
                        pr2[Y]=b[Y];
                    
                }
                else if(t0<1){
                    //duz
                    if(duzina_duzi(p,b)<duzina_duzi(q,d)){
                        pr1[X]=p[X];
                        pr1[Y]=p[Y];
                        pr2[X]=b[X];
                        pr2[Y]=b[Y];
                        
                    }
                    else{
                        pr1[X]=q[X];
                        pr1[Y]=q[Y];
                        pr2[X]=d[X];
                        pr2[Y]=d[Y];
                    }
                    
                }
            }
            
            return 'e';
        }
        
    }
    
    else if(imenilac==0 && brojilac2!=0){
        //paralelni su i ne seku se
        printf("paralelni su i nemaju presek\n");
        return '0';
        
        
        
    }
    
    else if(imenilac!=0 && t>=0 && t<=1 && u>=0 && u<=1){
        
        if(t==0 ||t==1 || u==0 || u==1){
            //kapiram da jedno teme pripada drugoj liniji
            printf("jedno teme pripada drugoj liniji\n");
            //treba da zadovolji jednacinu prave kroz 2 tacke
            printf("t je : %lf\tu je :%lf\n",t,u);
            if(jednacina_prave(b,q,d)){
                pr1[X]=b[X];
                pr1[Y]=b[Y];
                
            }
            else if(jednacina_prave(p,q,d)){
                pr1[X]=p[X];
                pr1[Y]=p[Y];
                
            }
            else if(jednacina_prave(q,p,d)){
                pr1[X]=q[X];
                pr1[Y]=q[Y];
            }
            else{
                pr1[X]=d[X];
                pr1[Y]=d[Y];
                
            }
            
            return 'v';
        }
        
        else{
        //seku se u p+t*r = q+u*s
        pr1[X]=p[X]+t*r[X];
        pr1[Y]=p[Y]+t*r[Y];
        printf("seku se u {%lf, %lf}\n",pr1[X],pr1[Y]);
        return '1';
        }
    }
    
    //i za ovo nisam bas siguran
   /* else if(imenilac!=0 && (t==0 ||t==1 || u==0 || u==1)){
        //kapiram da jedno teme pripada drugoj liniji
        printf("jedno teme pripada drugoj liniji\n");
        //treba da zadovolji jednacinu prave kroz 2 tacke
        printf("t je : %lf\tu je :%lf\n",t,u);
        if(jednacina_prave(b,q,d)){
            pr1[X]=b[X];
            pr1[Y]=b[Y];
        }
        else if(jednacina_prave(p,q,d)){
            pr1[X]=p[X];
            pr1[Y]=p[Y];
        }
        else if(jednacina_prave(q,p,d)){
            pr1[X]=q[X];
            pr1[Y]=q[Y];
            
        }
        else{
            pr1[X]=d[X];
            pr1[Y]=d[Y];
        }
        return 'v';
        
    }*/
    
    else{
        //nisu paralelni ali se ne seku
        printf("nisu paralelni ali se ne seku\n");
        return '0';
    }
    
    return 'o';
    
    
}

tUnFlag PostaviFleg(tTackad p, tUnFlag unflag, int aHB, int bHA){
    printf("***************{%lf,%lf}*************\n",p[X], p[Y]);
    fprintf(ispis, "%lf %lf\n",p[X],p[Y]);
    
    printf("Pun: %d, Qun: %d, Nepoznato: %d\n",Pun,Qun,Nepoznato);
    
    if(aHB>0)
        return Pun;
    else if (bHA>0)
        return Qun;
    else
        return unflag;
    
}

int Napreduj(int a, int *aa, int n, bool unutra, tTackai v){
    
    if(unutra){
        printf("---------------------------------------------------------------{%d,%d} linija do\n",v[X],v[Y]);
        fprintf(ispis, "%d %d\n",v[X],v[Y]);
    }
    (*aa)++;
    
    return (a+1)%n;    
    
}

int jednacina_prave(tTackai a, tTackai p, tTackai q){
    
    if (((double)(a[X]-p[X])/(double)(q[X]-p[X]))==((double)(a[Y]-p[Y])/(double)(q[Y]-p[Y])))
        return 1;
    return 0;
    
    
}

double duzina_duzi(tTackai a, tTackai b){
    
    return abs((double)((b[X]-a[X])*(b[X]-a[X]) - (b[Y]-a[Y])*(b[Y]-a[Y])));
    
    
}




void Vektor(tTackai a, tTackai b, tTackai c){
    
    int i;
    for(i=0;i<DIM;i++)
        c[i]=a[i]-b[i];
    
    
}

int Znak( tTackai a, tTackai b, tTackai c ){
    
    double rez;
    
    rez = ( b[0] - a[0] ) * (double)( c[1] - a[1] ) - (double)( b[1] - a[1] ) * ( c[0] - a[0]);
   // printf("-----rez je: %lf\n",rez);
    
    if(rez>0.5) return 1;
    else if(rez<-0.5) return -1;
    else return 0;
    
}


int proizvodX(tTackai v, tTackai w){
    
    return v[X]*w[Y]-v[Y]*w[X];
    
}

int proizvodD(tTackai v, tTackai w){
    
    int i;
    int rez =0;
    
    for(i=0;i<2;i++)
        rez+=v[i]*w[i];
    
    return rez;
}


int Nacrtaj(){
    
    /*char *line = NULL;
    size_t len = 0;
    int x, y;*/
    
    FILE *gnuplot;
    
    
    
    
    
    if((gnuplot = popen("gnuplot", "w"))==NULL){
        printf("Greska pri otvaranju gnuplota\n");
        exit(EXIT_FAILURE);
    }
        
    
    fprintf(gnuplot, "plot 'poligoni.txt' with linespoints linestyle 2\n");
    fflush(gnuplot);
    
    
    /*if((f=fopen("preseci.txt","w"))==NULL){
        printf("Greska pri otvaranju fajla za upis\n")
    
    while((getline(&line, &len, fptr)) != -1){
        sscanf(line, "%d %d", &x, &y);
        fprintf("%d %d\n",x, y);
    }*/
    /*  for (i = 0; i < count; i++)
        fprintf(gnuplot, "%g %g\n", c[0], c[1]);
    fprintf(gnuplot, "e\n");*/
    
    
    
    while(1){
        
    }
    
    pclose(gnuplot);
    
    
    exit(EXIT_SUCCESS);
    
}
